/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file contest.hpp
  \brief klut to xag method for IWLS contest

  \author Hanyu WANG
*/

#pragma once

#include "../../algorithms/balancing/sop_balancing.hpp"
#include "../../algorithms/cleanup.hpp"
#include "../../algorithms/cut_rewriting.hpp"
#include "../../algorithms/detail/resub_utils.hpp"
#include "../../algorithms/node_resynthesis/xag_npn.hpp"
#include "../../algorithms/permute_outputs.hpp"
#include "../../algorithms/simulation.hpp"
#include "../../networks/xag.hpp"
#include "../../networks/klut.hpp"

namespace mockturtle
{
struct tts_to_xags_params
{
  bool verbose{ false };
  bool permute_output{ false };
  uint32_t population{ 3u };
  uint32_t max_unate_pairs{ 1000 };
};

namespace detail
{

class tts_to_xags_impl
{
public:
  using TT = kitty::dynamic_truth_table;
  using signal = xag_network::signal;
  using node = xag_network::node;
  using resyn_t = std::function<void( xag_network& )>;

private:

  inline const double gibbs( TT const& on, TT const& off )
  {
    auto p = []( double x ) { return x == 0 ? 0 : x * x * log2( x ); }; /* this is good for XOR */
    double x = kitty::count_ones( on );
    double y = kitty::count_ones( off );
    return p( x + y ) - p( x ) - p( y );
  }

  /* sort output based on the distance to const */
  void sort_output()
  {
    std::iota( std::begin( outputs_order ), std::end( outputs_order ), 0 );
    if ( ps.permute_output )
    {
      std::sort( std::begin( outputs_order ), std::end( outputs_order ), [&]( auto x, auto y ) {
        double score_x = gibbs( targets[x], ~targets[x] );
        double score_y = gibbs( targets[y], ~targets[y] );
        return score_x < score_y;
      } );
    }
  }

  /* generate a initial network */
  template<class ContainerType>
  signal synthesis_rec( xag_network& xag, TT const& on, TT const& off, ContainerType& ntts )
  {
    if ( kitty::count_ones( on ) == 0 )
      return xag.get_constant( 0 ); /* prefer 0 for next AND gate */
    if ( kitty::count_ones( off ) == 0 )
      return !xag.get_constant( 0 );
    double best = gibbs( on, off );
    node n_;

    // find the best candidate to split
    xag.foreach_node( [&]( node n, uint32_t i ) {
      if ( xag.is_dead( n ) || xag.is_constant( n ) )
        return;
      if ( xag.value( n ) == tabu_marker )
        return;
      double curr = gibbs( on & ntts[n], off & ntts[n] ) + gibbs( on & ~ntts[n], off & ~ntts[n] );
      if ( curr < best )
      {
        best = curr;
        n_ = n;
      }
    } );
    signal pos_cofactor = synthesis_rec( xag, on & ntts[n_], off & ntts[n_], ntts );   /* f=DC if ntt=0 */
    signal neg_cofactor = synthesis_rec( xag, on & ~ntts[n_], off & ~ntts[n_], ntts ); /* f=DC if ntt=1 */
    return xag.create_or(
        xag.create_and( xag.make_signal( n_ ), pos_cofactor ),
        xag.create_and( !xag.make_signal( n_ ), neg_cofactor ) );
  }

  signal on_synthesis( xag_network& xag, TT const& on, TT const& off )
  {
    default_simulator<TT> nsim( xag.num_pis() );
    unordered_node_map<TT, xag_network> ntts( xag );
    simulate_nodes<TT>( xag, ntts, nsim );
    auto add_event = xag.events().register_add_event( [&]( const auto& n ) {
      ntts.resize();
      std::vector<TT> fanin_values( xag.fanin_size( n ) );
      xag.foreach_fanin( n, [&]( auto const& f, auto i ) {
        fanin_values[i] = ntts[xag.get_node( f )]; /* compute will take care of the complement */
      } );
      ntts[n] = xag.compute( n, std::begin( fanin_values ), std::end( fanin_values ) );
    } );

    for ( auto i = 0u; i < xag.num_pis(); ++i )
    {
      signal si = xag.make_signal( xag.pi_at( i ) );

      for ( auto j = i + 1; j < xag.num_pis(); ++j )
      {
        signal sj = xag.make_signal( xag.pi_at( j ) );

        xag.create_and( si, sj );
        xag.create_and( si, !sj );
        xag.create_and( !si, sj );
        xag.create_and( !si, !sj );
        xag.create_xor( si, sj );
      }
    }

    signal s = synthesis_rec( xag, on, off, ntts );
    if ( add_event )
      xag.events().release_add_event( add_event );
    return s;
  }

  signal on_resynthesis( xag_network& xag, node n ) /* resynthesis an exisiting node */
  {
    default_simulator<TT> nsim( xag.num_pis() );
    unordered_node_map<TT, xag_network> ntts( xag );
    simulate_nodes<TT>( xag, ntts, nsim );
    auto add_event = xag.events().register_add_event( [&]( const auto& n_ ) {
      ntts.resize();
      std::vector<TT> fanin_values( xag.fanin_size( n_ ) );
      xag.foreach_fanin( n_, [&]( auto const& f, auto i ) {
        fanin_values[i] = ntts[xag.get_node( f )]; /* compute will take care of the complement */
      } );
      ntts[n_] = xag.compute( n_, std::begin( fanin_values ), std::end( fanin_values ) );
    } );

    tabu_marker++;
    mark_fanout_cone( xag, n );

    signal s = synthesis_rec( xag, ntts[n], ~ntts[n], ntts );
    xag.substitute_node( n, s );
    if ( add_event )
      xag.events().release_add_event( add_event );
    return s;
  }

  uint32_t mark_fanout_cone( xag_network& xag, node n ) /* this will modify the trav id */
  {
    fanout_view fxag{ xag }; 
    typename mockturtle::detail::node_mffc_inside<xag_network> mffc_mgr( fxag );
    std::queue<node> q;
    return mffc_mgr.call_on_mffc_and_count( n, {}, [&]( node const& _n ) {
      q.push( n );
    } );
    fxag.incr_trav_id();
    while ( !q.empty() )
    {
      fxag.set_visited( q.front(), fxag.trav_id() );
      fxag.foreach_fanout( q.front(), [&]( const auto& f ) { if ( fxag.visited( f ) != fxag.trav_id() ) q.push( f ); } );
      fxag.set_value( q.front(), tabu_marker );
      q.pop();
    }
  }

public:
  template<class iterator_type>
  tts_to_xags_impl( iterator_type begin, iterator_type end, tts_to_xags_params const& ps, resyn_t const& resyn_fn )
      : ps( ps ), resyn_fn( resyn_fn )
  {
    while ( begin != end ) targets.emplace_back( *begin++ );
    pis.resize( targets[0].num_vars() );
    outputs_order.resize( targets.size() );
  }

  std::vector<xag_network> run()
  {
    std::vector<xag_network> xags;
    sort_output(); /* initial output */
    for ( uint32_t j = 0u; j < ps.population; j++ )
    {
      xag_network xag;
      std::generate( std::begin( pis ), std::end( pis ), [&]() { return xag.create_pi(); } );
      for ( auto i = 0u; i < targets.size(); i++ )
      {
        TT tt = targets[outputs_order[i]];
        xag.create_po( on_synthesis( xag, tt, ~tt ) );
        xag = cleanup_dangling( xag );
        resyn_fn( xag );
      }
      xag = permute_outputs_back( xag, outputs_order );
      std::shuffle( std::begin( outputs_order ), std::end( outputs_order ), std::mt19937() );
      xags.emplace_back( xag );
    }
    // for ( auto i = 0u; i < ntk.num_pos(); i++ )
    // {
    //   on_resynthesis( xag.get_node( xag.po_at( i ) ) );
    //   xag = cleanup_dangling( xag );
    //   resyn_fn( xag );
    // }
    return xags;
  }

private:
  std::vector<TT> targets;
  std::vector<uint64_t> outputs_order;
  tts_to_xags_params const& ps;

  resyn_t resyn_fn;

  std::vector<signal> pis;
  uint32_t tabu_marker{ 1u };
};
} // namespace detail

void null_resynthesis_xag( xag_network& xag )
{
  (void)xag;
  xag = cleanup_dangling( xag );
}

template<class Ntk>
std::vector<xag_network> klut_to_xags( Ntk const& ntk, tts_to_xags_params const& ps = {}, std::function<void( xag_network& )> const& resyn_fn = null_resynthesis_xag )
{
  using TT = kitty::dynamic_truth_table;
  std::vector<TT> tts;
  default_simulator<TT> sim( ntk.num_pis() );
  unordered_node_map<TT, Ntk> targets( ntk );
  simulate_nodes<TT>( ntk, targets, sim );
  for ( uint32_t i = 0; i < ntk.num_pos(); i++ )
  {
    TT tt = targets[ntk.get_node( ntk.po_at( i ) )];
    tts.emplace_back( ntk.is_complemented( ntk.po_at( i ) ) ? ~tt : tt );
  }
  detail::tts_to_xags_impl p( std::begin( tts ), std::end( tts ), ps, resyn_fn );
  const auto ret = p.run();
  return ret;
}

template<class Ntk>
xag_network klut_to_xag( Ntk const& ntk, tts_to_xags_params const& ps = {}, std::function<void( xag_network& )> const& resyn_fn = null_resynthesis_xag )
{
  using TT = kitty::dynamic_truth_table;
  std::vector<TT> tts;
  default_simulator<TT> sim( ntk.num_pis() );
  unordered_node_map<TT, Ntk> targets( ntk );
  simulate_nodes<TT>( ntk, targets, sim );
  for ( uint32_t i = 0; i < ntk.num_pos(); i++ )
  {
    TT tt = targets[ntk.get_node( ntk.po_at( i ) )];
    tts.emplace_back( ntk.is_complemented( ntk.po_at( i ) ) ? ~tt : tt );
  }
  detail::tts_to_xags_impl p( std::begin( tts ), std::end( tts ), ps, resyn_fn );
  const auto ret = p.run();
  return ret[0];
}

template<class iterator_type>
std::vector<xag_network> tt_to_xags( iterator_type begin, iterator_type end, tts_to_xags_params const& ps = {}, std::function<void( xag_network& )> const& resyn_fn = null_resynthesis_xag )
{
  detail::tts_to_xags_impl p( begin, end, ps, resyn_fn );
  const auto ret = p.run();
  return ret;
}

} // namespace mockturtle