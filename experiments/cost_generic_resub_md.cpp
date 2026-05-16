/* Multi-die-aware arrival time resub on i2c.
 *
 *   Md_AT(n) = max over fanin v of { Md_AT(v) + delay(v, n) }
 *   delay(v, n) = 1                if die(v) == die(n)
 *               = LARGE (= 10000)  otherwise
 *
 *  Forest-internal new nodes inherit the pivot's die (paper's policy). This
 *  driver demonstrates the post-Phase-1/2 cost-generic API:
 *    - cost function is a plain instance (no statics)
 *    - pivot context is delivered via the on_pivot hook
 *    - new source nodes are auto-labelled via the on_add hook
 */

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/experimental/cost_generic_resub.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/views/cost_view.hpp>

using namespace mockturtle;
using namespace mockturtle::experimental;

/* ---------------------------------------------------------------------- */
/* multi-die-aware arrival-time cost function                              */
/* ---------------------------------------------------------------------- */

template<class Ntk>
struct md_arrival_cost_function
{
  using node = typename Ntk::node;

  struct context_t
  {
    uint32_t md_at{ 0 };
    uint8_t die{ 0 };
  };

  static constexpr uint32_t LARGE = 10000u;

  /* User wires these once on the instance; copies (made by cost_view /
     cost_resyn) keep the same captures, so the underlying die_map is
     shared. */
  std::function<std::optional<uint8_t>( node )> source_die_get;
  std::function<void( node, uint8_t )>          source_die_set;

  /* Updated by cost_resyn at the start of each resub problem (Phase 2 hook). */
  mutable uint8_t pivot_die{ 0 };

  /* Phase 2 hook — receives the pivot's full context once per resub problem. */
  void on_pivot( context_t const& pctx ) const
  {
    pivot_die = pctx.die;
  }

  /* Phase 2.5 hook — fires for every node added to the MAIN cost_view's
     underlying network. cost_view auto-skips this for disposable solution
     forests (is_main == false), so the user code is now forest-agnostic. */
  void on_add( Ntk const&, node const& n, context_t const& ctx ) const
  {
    if ( source_die_set ) source_die_set( n, ctx.die );
  }

  /* Propagation. */
  context_t operator()( Ntk const& ntk, node const& n,
                        std::vector<context_t> const& fc = {} ) const
  {
    if ( ntk.is_constant( n ) ) return context_t{};

    if ( ntk.is_pi( n ) )
    {
      /* In the source ntk this is a true PI (known die from partitioning).
         In the forest, divisor PIs have their context written directly by
         cost_resyn via set_context, so this branch does not fire for them. */
      uint8_t d = source_die_get
                    ? source_die_get( n ).value_or( pivot_die )
                    : pivot_die;
      return { 0u, d };
    }

    /* Internal node: source nodes have a known die from the partitioner;
       forest-internal new nodes inherit the pivot's die. */
    auto opt = source_die_get ? source_die_get( n ) : std::nullopt;
    uint8_t my_die = opt.value_or( pivot_die );

    uint32_t arr = 0u;
    for ( auto const& f : fc )
    {
      uint32_t d = ( f.die == my_die ) ? 1u : LARGE;
      arr = std::max( arr, f.md_at + d );
    }
    return { arr, my_die };
  }

  /* Contribution. */
  void operator()( Ntk const&, node const&,
                   uint32_t& total, context_t const ctx ) const
  {
    total = std::max( total, ctx.md_at );
  }
};

/* ---------------------------------------------------------------------- */
/* mini BLIF parser (each .names is a 2-input AND in AIG form)             */
/* ---------------------------------------------------------------------- */

struct blif_result
{
  aig_network aig;
  std::unordered_map<std::string, aig_network::signal> name_to_signal;
  std::vector<std::string> po_order;
};

blif_result read_aig_blif( std::string const& path )
{
  blif_result r;
  auto& aig = r.aig;
  auto& nts = r.name_to_signal;

  std::ifstream fin( path );
  if ( !fin ) { std::cerr << "cannot open " << path << "\n"; std::exit( 1 ); }

  std::string line;
  while ( std::getline( fin, line ) )
  {
    if ( line.empty() || line[0] == '#' ) continue;

    std::istringstream iss( line );
    std::string tok;
    iss >> tok;

    if ( tok == ".model" ) continue;
    if ( tok == ".end" ) break;

    if ( tok == ".inputs" )
    {
      std::string n;
      while ( iss >> n ) nts[n] = aig.create_pi();
      continue;
    }

    if ( tok == ".outputs" )
    {
      std::string n;
      while ( iss >> n ) r.po_order.push_back( n );
      continue;
    }

    if ( tok == ".names" )
    {
      std::vector<std::string> nodes;
      std::string n;
      while ( iss >> n ) nodes.push_back( n );

      std::vector<std::string> cover;
      std::streampos prev;
      while ( true )
      {
        prev = fin.tellg();
        if ( !std::getline( fin, line ) ) break;
        if ( line.empty() ) continue;
        if ( line[0] == '.' ) { fin.seekg( prev ); break; }
        cover.push_back( line );
      }

      std::string const out = nodes.back();
      size_t const nin = nodes.size() - 1;

      if ( nin == 0 )
      {
        bool val = false;
        if ( cover.size() == 1 && cover[0].size() && cover[0][0] == '1' ) val = true;
        nts[out] = aig.get_constant( val );
        continue;
      }

      if ( nin == 1 )
      {
        assert( cover.size() == 1 );
        auto s = nts.at( nodes[0] );
        bool const inv = ( cover[0][0] == '0' );
        nts[out] = inv ? !s : s;
        continue;
      }

      assert( nin == 2 && cover.size() == 1 );
      auto a = nts.at( nodes[0] );
      auto b = nts.at( nodes[1] );
      bool const ia = ( cover[0][0] == '0' );
      bool const ib = ( cover[0][1] == '0' );
      nts[out] = aig.create_and( ia ? !a : a, ib ? !b : b );
    }
  }

  for ( auto const& po : r.po_order )
  {
    auto it = nts.find( po );
    assert( it != nts.end() );
    aig.create_po( it->second );
  }

  return r;
}

/* ---------------------------------------------------------------------- */
std::unordered_map<std::string, uint8_t> read_die_csv( std::string const& path )
{
  std::unordered_map<std::string, uint8_t> r;
  std::ifstream fin( path );
  if ( !fin ) { std::cerr << "cannot open " << path << "\n"; std::exit( 1 ); }
  std::string line;
  std::getline( fin, line ); /* header */
  while ( std::getline( fin, line ) )
  {
    auto c = line.find( ',' );
    if ( c == std::string::npos ) continue;
    r[line.substr( 0, c )] = (uint8_t)std::stoi( line.substr( c + 1 ) );
  }
  return r;
}

/* ---------------------------------------------------------------------- */
uint32_t count_sll( aig_network const& aig,
                    std::function<std::optional<uint8_t>( aig_network::node )> const& get )
{
  uint32_t sll = 0;
  aig.foreach_gate( [&]( auto const& n ) {
    auto dn = get( n );
    if ( !dn ) return;
    aig.foreach_fanin( n, [&]( auto const& f ) {
      auto fn = aig.get_node( f );
      if ( aig.is_constant( fn ) ) return;
      auto df = get( fn );
      if ( df && *df != *dn ) ++sll;
    } );
  } );
  return sll;
}

/* ---------------------------------------------------------------------- */
int main( int argc, char** argv )
{
  std::string const blif = ( argc > 1 ) ? argv[1] : "/tmp/i2c_data/i2c_aig.blif";
  std::string const csv  = ( argc > 2 ) ? argv[2] : "/tmp/i2c_data/i2c_die_labels.csv";

  /* 1. read BLIF */
  auto br = read_aig_blif( blif );
  auto& aig = br.aig;
  std::cout << "[i] AIG: " << aig.num_pis() << " PIs, " << aig.num_pos()
            << " POs, " << aig.num_gates() << " gates\n";

  /* 2. read CSV → name→die map */
  auto n2d = read_die_csv( csv );
  std::cout << "[i] CSV: " << n2d.size() << " entries\n";

  /* 3. build die_map.
     We use std::unordered_map so that reads for unmapped node ids return a
     default-constructed std::optional (= nullopt). With node_map's vector
     backing, indexing past the current size is UB, and propagation reads
     before our on_add hook gets a chance to resize. */
  auto die_map = std::make_shared<std::unordered_map<aig_network::node, uint8_t>>();
  uint32_t mapped = 0;
  for ( auto const& [name, sig] : br.name_to_signal )
  {
    auto it = n2d.find( name );
    if ( it != n2d.end() )
    {
      ( *die_map )[aig.get_node( sig )] = it->second;
      ++mapped;
    }
  }
  std::cout << "[i] die_map: " << mapped << " mapped\n";

  /* 4. set up cost function instance (no statics; lambdas hold die_map) */
  using cost_fn_t = md_arrival_cost_function<aig_network>;
  cost_fn_t costfn;
  costfn.source_die_get = [die_map]( aig_network::node n ) -> std::optional<uint8_t> {
    auto it = die_map->find( n );
    if ( it == die_map->end() ) return std::nullopt;
    return it->second;
  };
  costfn.source_die_set = [die_map]( aig_network::node n, uint8_t d ) {
    ( *die_map )[n] = d;
  };

  auto get_die = costfn.source_die_get;

  /* 5. evaluate cost before */
  auto eval_md_at = [&]() {
    cost_view<aig_network, cost_fn_t> v( aig, costfn );
    return v.get_cost();
  };

  uint32_t const md_at_before = eval_md_at();
  uint32_t const sll_before = count_sll( aig, get_die );
  uint32_t const size_before = aig.num_gates();
  std::cout << "[i] before: Md_AT = " << md_at_before
            << "  SLL = " << sll_before
            << "  size = " << size_before << "\n";

  /* 6. run cost_generic_resub */
  cost_generic_resub_params ps;
  cost_generic_resub_stats st;
  ps.verbose = true;
  ps.wps.max_pis = 8;
  ps.wps.max_divisors = 150;
  ps.wps.max_inserts = 2;

  stopwatch<>::duration time_resub{ 0 };
  {
    stopwatch t( time_resub );
    cost_generic_resub( aig, costfn, ps, &st );
  }

  /* 7. evaluate after (die_map kept up-to-date by on_add hook) */
  uint32_t const md_at_after = eval_md_at();
  uint32_t const sll_after = count_sll( aig, get_die );
  uint32_t const size_after = aig.num_gates();
  auto cleaned = cleanup_dangling( aig );

  std::cout << "[i] after:  Md_AT = " << md_at_after
            << "  SLL = " << sll_after
            << "  size = " << size_after
            << "  (post-cleanup size = " << cleaned.num_gates() << ")\n";
  std::cout << "[i] gains:  dMd_AT = " << (int64_t)md_at_after - (int64_t)md_at_before
            << "  dSLL = " << (int64_t)sll_after - (int64_t)sll_before
            << "  dSize = " << (int64_t)size_after - (int64_t)size_before
            << "  runtime = " << to_seconds( time_resub ) << " s\n";

  return 0;
}
