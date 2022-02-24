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
  \file xag_resyn.hpp
  \brief Resynthesis by recursive decomposition for AIGs or XAGs.
  (based on ABC's implementation in `giaResub.c` by Alan Mishchenko)

  \author Siang-Yun Lee
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/node_map.hpp"

#include <kitty/kitty.hpp>
#include <fmt/format.h>
#include <abcresub/abcresub.hpp>

#include <vector>
#include <algorithm>
#include <type_traits>
#include <optional>
#include <queue>
#include <unordered_set>

namespace mockturtle
{

struct xag_resyn_static_params
{
  using base_type = xag_resyn_static_params;

  /*! \brief Maximum number of binate divisors to be considered. */
  static constexpr uint32_t max_binates{50u};

  /*! \brief Reserved capacity for divisor truth tables (number of divisors). */
  static constexpr uint32_t reserve{200u};

  /*! \brief Whether to consider single XOR gates (i.e., using XAGs instead of AIGs). */
  static constexpr bool use_xor{true};

  /*! \brief Whether to copy truth tables. */
  static constexpr bool copy_tts{false};

  /*! \brief Whether to preserve depth. */
  static constexpr bool preserve_depth{false};

  /*! \brief Whether the divisors have uniform costs (size and depth, whenever relevant). */
  static constexpr bool uniform_div_cost{true};

  /*! \brief Size cost of each AND gate. */
  static constexpr uint32_t size_cost_of_and{1u};

  /*! \brief Size cost of each XOR gate (only relevant when `use_xor = true`). */
  static constexpr uint32_t size_cost_of_xor{1u};

  /*! \brief Depth cost of each AND gate (only relevant when `preserve_depth = true`). */
  static constexpr uint32_t depth_cost_of_and{1u};

  /*! \brief Depth cost of each XOR gate (only relevant when `preserve_depth = true` and `use_xor = true`). */
  static constexpr uint32_t depth_cost_of_xor{1u};

  static constexpr uint32_t max_enqueue{1000u};

  static constexpr uint32_t max_xor{1u};

  static constexpr uint32_t max_neighbors{10u};

  using truth_table_storage_type = void;
  using node_type = void;
};

template<class TT>
struct xag_resyn_static_params_default : public xag_resyn_static_params
{
  using truth_table_storage_type = std::vector<TT>;
  using node_type = uint32_t;
};

template<class TT>
struct aig_resyn_static_params_default : public xag_resyn_static_params_default<TT>
{
  static constexpr bool use_xor = false;
};

template<class TT>
struct xag_resyn_static_params_preserve_depth : public xag_resyn_static_params_default<TT>
{
  static constexpr bool preserve_depth = true;
  static constexpr bool uniform_div_cost = false;
};

template<class TT>
struct aig_resyn_static_params_preserve_depth : public xag_resyn_static_params_preserve_depth<TT>
{
  static constexpr bool use_xor = false;
};

template<class Ntk>
struct xag_resyn_static_params_for_sim_resub : public xag_resyn_static_params
{
  using truth_table_storage_type = incomplete_node_map<kitty::partial_truth_table, Ntk>;
  using node_type = typename Ntk::node;
};

template<class Ntk>
struct aig_resyn_static_params_for_sim_resub : public xag_resyn_static_params_for_sim_resub<Ntk>
{
  static constexpr bool use_xor = false;
};

struct xag_resyn_stats
{
  /*! \brief Time for finding 0-resub and collecting unate literals. */
  stopwatch<>::duration time_unate{0};

  /*! \brief Time for finding 1-resub. */
  stopwatch<>::duration time_resub1{0};

  /*! \brief Time for finding 2-resub. */
  stopwatch<>::duration time_resub2{0};

  /*! \brief Time for finding 3-resub. */
  stopwatch<>::duration time_resub3{0};

  /*! \brief Time for sorting unate literals and unate pairs. */
  stopwatch<>::duration time_sort{0};

  /*! \brief Time for collecting unate pairs. */
  stopwatch<>::duration time_collect_pairs{0};

  /*! \brief Time for dividing the target and recursive call. */
  stopwatch<>::duration time_divide{0};

  /*! \brief Time for checking unateness of literals. */
  stopwatch<>::duration time_check_unateness{0};

  stopwatch<>::duration time_enqueue{0};
  stopwatch<>::duration time_tt_calculation{0};
  stopwatch<>::duration time_check_unate{0};
  stopwatch<>::duration time_move_tt{0};

  void report() const
  {
    fmt::print( "[i]         <xag_resyn_decompose>\n" );
    fmt::print( "[i]             0-resub      : {:>5.2f} secs\n", to_seconds( time_unate ) );
    fmt::print( "[i]             1-resub      : {:>5.2f} secs\n", to_seconds( time_resub1 ) );
    fmt::print( "[i]             2-resub      : {:>5.2f} secs\n", to_seconds( time_resub2 ) );
    fmt::print( "[i]             3-resub      : {:>5.2f} secs\n", to_seconds( time_resub3 ) );
    fmt::print( "[i]             sort         : {:>5.2f} secs\n", to_seconds( time_sort ) );
    fmt::print( "[i]             collect pairs: {:>5.2f} secs\n", to_seconds( time_collect_pairs ) );
    fmt::print( "[i]             dividing     : {:>5.2f} secs\n", to_seconds( time_divide ) );
    fmt::print( "[i]             enqueue      : {:>5.2f} secs\n", to_seconds( time_enqueue ) );
    fmt::print( "[i]             tt calc      : {:>5.2f} secs\n", to_seconds( time_tt_calculation ) );
    fmt::print( "[i]             check unate  : {:>5.2f} secs\n", to_seconds( time_check_unate ) );
    fmt::print( "[i]             tt move      : {:>5.2f} secs\n", to_seconds( time_move_tt ) );
  }
};

/*! \brief Logic resynthesis engine for AIGs or XAGs.
 *
 * The algorithm is based on ABC's implementation in `giaResub.c` by Alan Mishchenko.
 *
 * Divisors are classified as positive unate (not overlapping with target offset),
 * negative unate (not overlapping with target onset), or binate (overlapping with
 * both onset and offset). Furthermore, pairs of binate divisors are combined with
 * an AND operation and considering all possible input polarities and again classified
 * as positive unate, negative unate or binate. Simple solutions of zero cost 
 * (one unate divisor), one node (two unate divisors), two nodes (one unate divisor + 
 * one unate pair), and three nodes (two unate pairs) are exhaustively examined.
 * When no simple solutions can be found, the algorithm heuristically chooses an unate
 * divisor or an unate pair to divide the target function with and recursively calls
 * itself to decompose the remainder function.
   \verbatim embed:rst

   Example

   .. code-block:: c++

      using TT = kitty::static_truth_table<6>;
      const std::vector<aig_network::node> divisors = ...;
      const node_map<TT, aig_network> tts = ...;
      const TT target = ..., care = ...;
      xag_resyn_stats st;
      xag_resyn_decompose<TT, node_map<TT, aig_network>, false, false, aig_network::node> resyn( st );
      auto result = resyn( target, care, divisors.begin(), divisors.end(), tts );
   \endverbatim
 */
template<class TT, class static_params = xag_resyn_static_params_default<TT>>
class xag_resyn_decompose
{
public:
  using stats = xag_resyn_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;

private:
  struct unate_lit
  {
    unate_lit( uint32_t l )
      : lit( l )
    { }

    bool operator==( unate_lit const& other ) const
    {
      return lit == other.lit;
    }

    uint32_t lit;
    uint32_t score{0};
  };

  struct fanin_pair
  {
    fanin_pair( uint32_t l1, uint32_t l2 )
      : lit1( l1 < l2 ? l1 : l2 ), lit2( l1 < l2 ? l2 : l1 )
    { }

    fanin_pair( uint32_t l1, uint32_t l2, bool is_xor )
      : lit1( l1 > l2 ? l1 : l2 ), lit2( l1 > l2 ? l2 : l1 )
    { (void)is_xor; }

    bool operator==( fanin_pair const& other ) const
    {
      return lit1 == other.lit1 && lit2 == other.lit2;
    }

    uint32_t lit1, lit2;
    uint32_t score{0};
  };

  enum gate_type {AND, OR, XOR, NONE};
  enum lit_type {EQUAL, EQUAL_INV, POS_UNATE, NEG_UNATE, POS_UNATE_INV, NEG_UNATE_INV, BINATE, DONT_CARE};

  std::unordered_map<TT, uint32_t, kitty::hash<TT>> tt_to_id;
  std::vector<TT> id_to_tt;
  std::vector<uint32_t> id_to_num;
  
  uint32_t to_id ( const TT & tt )
  {
    if ( tt_to_id.find( tt ) != tt_to_id.end() ) return tt_to_id[tt];
    tt_to_id[tt] = id_to_tt.size();
    id_to_tt.emplace_back( tt );
    id_to_num.emplace_back( kitty::count_ones( tt ) );
    return tt_to_id[tt];
  }

  const auto & to_tt ( uint32_t id )
  {
    assert( id < id_to_tt.size() );
    return id_to_tt[id];
  }

  const uint32_t to_num ( uint32_t id )
  {
    assert( id < id_to_num.size() );
    return id_to_num[id];
  }

  std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> best_cost;

  bool check_cost ( uint32_t x, uint32_t y, uint32_t cost )
  {
    uint32_t _x = std::min( x, y );
    uint32_t _y = std::max( x, y );
    bool ret = false; /* if best is updated */
    if ( best_cost.find( _x ) == best_cost.end() )
    {
      ret = true;
      best_cost[ _x ] = std::unordered_map<uint32_t, uint32_t>();
    }
    if ( best_cost[ _x ].find( _y ) == best_cost[ _x ].end() )
    {
      ret = true;
      best_cost[ _x ][ _y ] = cost;
    }
    if ( best_cost[ _x ][ _y ] > cost )
    {
      ret = true;
      best_cost[ _x ][ _y ] = cost;
    }
    return ret;
  }

  struct task 
  {
    std::array<uint32_t, 2> sets; /* the on-off set of each task (could be optimized) */
    std::pair<uint32_t, uint32_t> cost; /* the lower bound of the cost */
    uint32_t score;
    std::size_t prev;
    bool done;
    gate_type ntype;
    uint32_t lit;
    uint32_t num_xor;
    const bool operator > ( const task & t ) const
    {
      if ( cost.first != t.cost.first ) /* area */
      {
        return cost.first > t.cost.first;
      }
      /* the most likely first */
      return score > t.score;
    }

    task( auto _done, auto _prev, auto _lit, auto _ntype, auto _cost ): done(_done), prev(_prev), lit(_lit), ntype(_ntype), cost(_cost), score(0), num_xor(0) {}
  };

  struct deq_task
  {
    std::pair<uint32_t, uint32_t> cost; /* the lower bound of the cost */
    std::size_t prev;
    gate_type ntype;
    uint32_t lit;
    deq_task( const task & t ): prev( t.prev ), ntype( t.ntype ), lit( t.lit ), cost( t.cost ) {}
  };

public:
  explicit xag_resyn_decompose( stats& st ) noexcept
    : st( st )
  {
    static_assert( std::is_same_v<typename static_params::base_type, xag_resyn_static_params>, "Invalid static_params type" );
    static_assert( !( static_params::uniform_div_cost && static_params::preserve_depth ), "If depth is to be preserved, divisor depth cost must be provided (usually not uniform)" );
    divisors.reserve( static_params::reserve );
  }

  /*! \brief Perform XAG resynthesis.
   *
   * `tts[*begin]` must be of type `TT`.
   * Moreover, if `static_params::copy_tts = false`, `*begin` must be of type `static_params::node_type`.
   *
   * \param target Truth table of the target function.
   * \param care Truth table of the care set.
   * \param begin Begin iterator to divisor nodes.
   * \param end End iterator to divisor nodes.
   * \param tts A data structure (e.g. std::vector<TT>) that stores the truth tables of the divisor functions.
   * \param max_size Maximum number of nodes allowed in the dependency circuit.
   */
  template<class iterator_type, 
           bool enabled = static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );

    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( (*ptts)[*begin] );
      }
      else
      {
        divisors.emplace_back( *begin );
      }
      ++begin;
    }

    return compute_function( max_size );
  }

  template<class iterator_type, class Fn, 
           bool enabled = !static_params::uniform_div_cost && !static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max() )
  {}

  template<class iterator_type, class Fn, 
           bool enabled = !static_params::uniform_div_cost && static_params::preserve_depth, typename = std::enable_if_t<enabled>>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, typename static_params::truth_table_storage_type const& tts, Fn&& size_cost, Fn&& depth_cost, uint32_t max_size = std::numeric_limits<uint32_t>::max(), uint32_t max_depth = std::numeric_limits<uint32_t>::max() )
  {
    static_assert( static_params::copy_tts || std::is_same_v<typename std::iterator_traits<iterator_type>::value_type, typename static_params::node_type>, "iterator_type does not dereference to static_params::node_type" );

    ptts = &tts;

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    while ( begin != end )
    {
      if constexpr ( static_params::copy_tts )
      {
        divisors.emplace_back( (*ptts)[*begin] );
      }
      else
      {
        divisors.emplace_back( *begin );
      }
      ++begin;
    }

    index_list.clear();
    index_list.add_inputs( divisors.size() - 1 );

    /* check trivial solution */
    if ( kitty::count_ones( ~target & care ) == 0 )
    {
      index_list.add_output( 1 );
      return index_list;
    }
    if ( kitty::count_ones( target & care ) == 0 )
    {
      index_list.add_output( 0 );
      return index_list;
    }

    depth_fn = depth_cost;
    upper_bound = max_size;

    mem.clear();
    id_to_num.clear();
    id_to_tt.clear();
    tt_to_id.clear();

    best_cost.clear();

    task init_task( false, 0, 0, NONE, std::pair(0,0) );
    init_task.sets[0] = to_id( ~target & care );
    init_task.sets[1] = to_id( target & care );

    std::priority_queue<task, std::vector<task>, std::greater<>> q;
    /* prepare the initial task */
    call_with_stopwatch( st.time_enqueue, [&]() {
      q.push( init_task );
    } );

    while ( !q.empty() ) 
    {
      /* get the current lower bound */
      auto t = q.top(); q.pop(); 
      mem.emplace_back( deq_task( t ) );
      /* back trace succeed tasks */
      if ( t.done == true )
      {
        auto output = back_trace( mem.size() - 1 );
        index_list.add_output( output.second );
        return index_list;
      }
      if ( t.cost.first >= upper_bound ) break;
      if ( q.size() >= static_params::max_enqueue ) break;

      add_neighbors ( t, q );
    }
    return std::nullopt;
  }

private:
  /* the list of temp nodes */
  std::vector<deq_task> mem;

  /* the depth cost function */
  std::function<uint32_t(uint32_t)> depth_fn;

  /* cost upper bound */
  uint32_t upper_bound;

  template<class Q>
  auto add_neighbors ( const task & t, Q & q )
  {
    for ( auto v = 1u; v < divisors.size(); ++v )
    {
      auto _t = call_with_stopwatch( st.time_tt_calculation, [&]() {
        return find_unate_subtask( t, v );
      } );
      if ( _t ) /* prune dont cares */
      {
        if ( _t->done == true )
        {
          upper_bound = (_t->cost).first;
        }
        call_with_stopwatch( st.time_enqueue, [&]() {
          q.push( *_t );
        } );
      }
    }
  }

  /* */
  auto tt_move ( uint32_t off, uint32_t on, uint32_t lit, gate_type ntype )
  {
    const auto & tt = lit & 0x1? ~get_div( lit>>1 ) : get_div( lit>>1 ) ;
    uint32_t _off = 0u;
    uint32_t _on  = 0u;
    switch ( ntype )
    {
    case OR:
      _off = off;
      _on = to_id( ~tt & to_tt( on ) );
      break;
    case AND:
      _off = to_id( tt & to_tt( off ) );
      _on = on;
      break;
    case XOR:
      _off = to_id( ( ~tt & to_tt( off ) ) | ( tt & to_tt( on ) ) );
      _on = to_id( ( ~tt & to_tt( on ) ) | ( tt & to_tt( off ) ) );
      break;
    }
    return std::make_tuple( _off, _on );
  }

  using cand_t = std::pair<uint32_t, uint32_t>;
  cand_t back_trace( size_t pos )
  {
    size_t p = pos;
    std::priority_queue<cand_t, std::vector<cand_t>, std::greater<>> cand_q;
    cand_q.push( std::pair( depth_fn( mem[p].lit >> 1 ), mem[p].lit ) );
    while ( mem[p].prev != 0 )
    {
      for ( p=mem[p].prev; ; p=mem[p].prev )
      {
        cand_q.push( std::pair( depth_fn( mem[p].lit >> 1 ), mem[p].lit ) );
        if ( mem[p].ntype != mem[mem[p].prev].ntype ) break;
      }
      /* add the nodes */
      while ( cand_q.size() > 1 )
      {
        auto fanin1 = cand_q.top(); cand_q.pop();
        auto fanin2 = cand_q.top(); cand_q.pop();
        uint32_t new_lit = 0u;
        if ( mem[p].ntype == AND ) new_lit = index_list.add_and( fanin1.second, fanin2.second );
        else if ( mem[p].ntype == OR ) new_lit = index_list.add_and( fanin1.second ^ 0x1, fanin2.second ^ 0x1) ^ 0x1;
        else if ( mem[p].ntype == XOR ) new_lit = index_list.add_xor( fanin1.second, fanin2.second );
        auto new_cost = fanin2.first + 1; // TODO: change this "1" to cost
        cand_q.push( std::pair( new_cost, new_lit ) );
      }
    }
    return cand_q.top();
  }

  std::optional<std::pair<uint32_t, uint32_t>> get_cost( size_t pos, uint32_t lit, gate_type _ntype, bool balancing = false ) const 
  {
    uint32_t size_cost, depth_cost;
    size_cost = mem[pos].cost.first; 
    if ( mem[pos].ntype != NONE) // TODO: size cost is not uniform
    {
      size_cost ++;
    }
    if ( balancing ) /* a better estimation of depth cost */
    {
      std::priority_queue<uint32_t, std::vector<uint32_t>, std::greater<>> cost_q;
      auto p = pos;
      cost_q.push( depth_fn( mem[p].lit>>1 ) );
      while ( mem[p].prev != 0 )
      {
        for ( p = mem[p].prev; ; p = mem[p].prev ) /* get all the same node type */
        {
          cost_q.push( depth_fn( mem[p].lit >> 1 ) );
          if ( mem[p].ntype != mem[mem[p].prev].ntype ) break;
        }
        while ( cost_q.size() > 1 ) /* add node while maintaining the depth order */
        {
          cost_q.pop();
          cost_q.push( cost_q.top() + 1 );
          cost_q.pop();
        }
      }
      depth_cost = cost_q.top();
    }
    else
    {
      depth_cost = std::max( mem[pos].cost.second, (uint32_t)depth_fn( lit>>1 ) ) + 1;
    }
    return std::pair( size_cost, depth_cost );
  };

  lit_type check_unateness ( const TT & off_set, const TT & on_set, const TT & tt )
  {
    bool unateness[4] = {
      kitty::intersection_is_empty<TT, 1, 1>( tt, off_set ),
      kitty::intersection_is_empty<TT, 0, 1>( tt, off_set ),
      kitty::intersection_is_empty<TT, 1, 1>( tt, on_set ),
      kitty::intersection_is_empty<TT, 0, 1>( tt, on_set ),
    };
    if ( ( unateness[0] && unateness[2] ) || ( unateness[1] && unateness[3] ) ) return DONT_CARE;
    if ( unateness[0] && unateness[3] ) return EQUAL;
    if ( unateness[1] && unateness[2] ) return EQUAL_INV;
    if ( unateness[0] ) return POS_UNATE;
    if ( unateness[1] ) return POS_UNATE_INV;
    if ( unateness[2] ) return NEG_UNATE_INV;
    if ( unateness[3] ) return NEG_UNATE;
    return BINATE;
  }
  std::optional<task> find_unate_subtask ( const task & _t, uint32_t v )
  {
    auto const & tt = get_div(v);
    auto off = _t.sets[0];
    auto on  = _t.sets[1];
    auto ltype = call_with_stopwatch( st.time_check_unate, [&] () {
      return check_unateness( to_tt( off ), to_tt( on ), tt );
    } );
    gate_type _ntype = NONE;
    bool done = false;
    uint32_t lit = v << 1;
    switch ( ltype )
    {
    case DONT_CARE: return std::nullopt;
    case EQUAL:     done = true; _ntype = NONE;         break;
    case EQUAL_INV: done = true; _ntype = NONE; lit++ ; break;
    case POS_UNATE:              _ntype = OR;           break;
    case POS_UNATE_INV:          _ntype = OR;   lit++;  break;
    case NEG_UNATE:              _ntype = AND;          break;
    case NEG_UNATE_INV:          _ntype = AND;  lit++;  break;
    case BINATE:                 _ntype = XOR;          break;
    }

    if constexpr ( static_params::use_xor == false )
    {
      if ( _ntype == XOR ) return std::nullopt;
    }

    if ( _ntype != NONE && _ntype == _t.ntype && ( lit >> 1 ) <= ( _t.lit >> 1 ) ) 
    {
      return std::nullopt; /* commutativity */
    } 
    auto cost = get_cost( mem.size() - 1, lit, _ntype, false ); // TODO: assume prev task always at back()
    if ( !cost || (*cost).first >= upper_bound )
    {
      return std::nullopt; /* task is pruned */
    }
    task t( done, mem.size() - 1, lit, _ntype, *cost );
    if ( _ntype == XOR ) 
    {
      if ( _t.num_xor >= static_params::max_xor ) return std::nullopt;
      t.num_xor = _t.num_xor + 1;
    }
    if ( done == false )
    {
      auto [ _off, _on ] = call_with_stopwatch( st.time_move_tt, [&] () {
        return tt_move( off, on, lit, _ntype );
      } );
      if ( check_cost( _off, _on, (*cost).first ) == false )
      {
        return std::nullopt;
      }
      t.sets[0] = _off;
      t.sets[1] = _on;
      t.score = to_num( _off ) + to_num( _on );
    }
    return t;
  }

  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    index_list.clear();
    index_list.add_inputs( divisors.size() - 1 );
    auto const lit = compute_function_rec( num_inserts );
    if ( lit )
    {
      assert( index_list.num_gates() <= num_inserts );
      index_list.add_output( *lit );
      return index_list;
    }
    return std::nullopt;
  }

  std::optional<uint32_t> compute_function_rec( uint32_t num_inserts )
  {
    pos_unate_lits.clear();
    neg_unate_lits.clear();
    binate_divs.clear();
    pos_unate_pairs.clear();
    neg_unate_pairs.clear();

    /* try 0-resub and collect unate literals */
    auto const res0 = call_with_stopwatch( st.time_unate, [&]() {
      return find_one_unate();
    });
    if ( res0 )
    {
      return *res0;
    }
    if ( num_inserts == 0u )
    {
      return std::nullopt;
    }

    /* sort unate literals and try 1-resub */
    call_with_stopwatch( st.time_sort, [&]() {
      sort_unate_lits( pos_unate_lits, 1 );
      sort_unate_lits( neg_unate_lits, 0 );
    });
    auto const res1or = call_with_stopwatch( st.time_resub1, [&]() {
      return find_div_div( pos_unate_lits, 1 );
    });
    if ( res1or )
    {
      return *res1or;
    }
    auto const res1and = call_with_stopwatch( st.time_resub1, [&]() {
      return find_div_div( neg_unate_lits, 0 );
    });
    if ( res1and )
    {
      return *res1and;
    }

    if ( binate_divs.size() > static_params::max_binates )
    {
      binate_divs.resize( static_params::max_binates );
    }

    if constexpr ( static_params::use_xor )
    {
      /* collect XOR-type unate pairs and try 1-resub with XOR */
      auto const res1xor = find_xor();
      if ( res1xor )
      {
        return *res1xor;
      }
    }
    if ( num_inserts == 1u )
    {
      return std::nullopt;
    }

    /* collect AND-type unate pairs and sort (both types), then try 2- and 3-resub */
    call_with_stopwatch( st.time_collect_pairs, [&]() {
      collect_unate_pairs();
    });
    call_with_stopwatch( st.time_sort, [&]() {
      sort_unate_pairs( pos_unate_pairs, 1 );
      sort_unate_pairs( neg_unate_pairs, 0 );
    });
    auto const res2or = call_with_stopwatch( st.time_resub2, [&]() {
      return find_div_pair( pos_unate_lits, pos_unate_pairs, 1 );
    });
    if ( res2or )
    {
      return *res2or;
    }
    auto const res2and = call_with_stopwatch( st.time_resub2, [&]() {
      return find_div_pair( neg_unate_lits, neg_unate_pairs, 0 );
    });
    if ( res2and )
    {
      return *res2and;
    }

    if ( num_inserts >= 3u )
    {
      auto const res3or = call_with_stopwatch( st.time_resub3, [&]() {
        return find_pair_pair( pos_unate_pairs, 1 );
      });
      if ( res3or )
      {
        return *res3or;
      }
      auto const res3and = call_with_stopwatch( st.time_resub3, [&]() {
        return find_pair_pair( neg_unate_pairs, 0 );
      });
      if ( res3and )
      {
        return *res3and;
      }
    }

    /* choose something to divide and recursive call on the remainder */
    /* Note: dividing = AND the on-set (if using positive unate) or the off-set (if using negative unate)
                        with the *negation* of the divisor/pair (subtracting) */
    uint32_t on_off_div, on_off_pair;
    uint32_t score_div = 0, score_pair = 0;

    call_with_stopwatch( st.time_divide, [&]() {
      if ( pos_unate_lits.size() > 0 )
      {
        on_off_div = 1; /* use pos_lit */
        score_div = pos_unate_lits[0].score;
        if ( neg_unate_lits.size() > 0 && neg_unate_lits[0].score > pos_unate_lits[0].score )
        {
          on_off_div = 0; /* use neg_lit */
          score_div = neg_unate_lits[0].score;
        }
      }
      else if ( neg_unate_lits.size() > 0 )
      {
        on_off_div = 0; /* use neg_lit */
        score_div = neg_unate_lits[0].score;
      }

      if ( num_inserts > 3u )
      {
        if ( pos_unate_pairs.size() > 0 )
        {
          on_off_pair = 1; /* use pos_pair */
          score_pair = pos_unate_pairs[0].score;
          if ( neg_unate_pairs.size() > 0 && neg_unate_pairs[0].score > pos_unate_pairs[0].score )
          {
            on_off_pair = 0; /* use neg_pair */
            score_pair = neg_unate_pairs[0].score;
          }
        }
        else if ( neg_unate_pairs.size() > 0 )
        {
          on_off_pair = 0; /* use neg_pair */
          score_pair = neg_unate_pairs[0].score;
        }
      }
    });

    if ( score_div > score_pair / 2 ) /* divide with a divisor */
    {
      /* if using pos_lit (on_off_div = 1), modify on-set and use an OR gate on top;
         if using neg_lit (on_off_div = 0), modify off-set and use an AND gate on top
       */
      uint32_t const lit = on_off_div ? pos_unate_lits[0].lit : neg_unate_lits[0].lit;
      call_with_stopwatch( st.time_divide, [&]() {
        on_off_sets[on_off_div] &= lit & 0x1 ? get_div( lit >> 1 ) : ~get_div( lit >> 1 );
      });

      auto const res_remain_div = compute_function_rec( num_inserts - 1 );
      if ( res_remain_div )
      {
        auto const new_lit = index_list.add_and( ( lit ^ 0x1 ), *res_remain_div ^ on_off_div );
        return new_lit + on_off_div;
      }
    }
    else if ( score_pair > 0 ) /* divide with a pair */
    {
      fanin_pair const pair = on_off_pair ? pos_unate_pairs[0] : neg_unate_pairs[0];
      call_with_stopwatch( st.time_divide, [&]() {
        if constexpr ( static_params::use_xor )
        {
          if ( pair.lit1 > pair.lit2 ) /* XOR pair: ~(lit1 ^ lit2) = ~lit1 ^ lit2 */
          {
            on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) )
                                  ^ ( pair.lit2 & 0x1 ? ~get_div( pair.lit2 >> 1 ) : get_div( pair.lit2 >> 1 ) );
          }
          else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
          {
            on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) )
                                  | ( pair.lit2 & 0x1 ? get_div( pair.lit2 >> 1 ) : ~get_div( pair.lit2 >> 1 ) );
          }
        }
        else /* AND pair: ~(lit1 & lit2) = ~lit1 | ~lit2 */
        {
          on_off_sets[on_off_pair] &= ( pair.lit1 & 0x1 ? get_div( pair.lit1 >> 1 ) : ~get_div( pair.lit1 >> 1 ) )
                                | ( pair.lit2 & 0x1 ? get_div( pair.lit2 >> 1 ) : ~get_div( pair.lit2 >> 1 ) );
        }
      });

      auto const res_remain_pair = compute_function_rec( num_inserts - 2 );
      if ( res_remain_pair )
      {
        uint32_t new_lit1;
        if constexpr ( static_params::use_xor )
        {
          new_lit1 = ( pair.lit1 > pair.lit2 ) ? index_list.add_xor( pair.lit1, pair.lit2 ) : index_list.add_and( pair.lit1, pair.lit2 );
        }
        else
        {
          new_lit1 = index_list.add_and( pair.lit1, pair.lit2 );
        }
        auto const new_lit2 = index_list.add_and( new_lit1 ^ 0x1, *res_remain_pair ^ on_off_pair );
        return new_lit2 + on_off_pair;
      }
    }

    return std::nullopt;
  }

  /* See if there is a constant or divisor covering all on-set bits or all off-set bits.
     1. Check constant-resub
     2. Collect unate literals
     3. Find 0-resub (both positive unate and negative unate) and collect binate (neither pos nor neg unate) divisors
   */
  std::optional<uint32_t> find_one_unate()
  {
    num_bits[0] = kitty::count_ones( on_off_sets[0] ); /* off-set */
    num_bits[1] = kitty::count_ones( on_off_sets[1] ); /* on-set */
    if ( num_bits[0] == 0 )
    {
      return 1;
    }
    if ( num_bits[1] == 0 )
    {
      return 0;
    }

    for ( auto v = 1u; v < divisors.size(); ++v )
    {
      bool unateness[4] = {false, false, false, false};
      /* check intersection with off-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div(v), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 );
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div(v), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[1] = true;
      }

      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div(v), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 );
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div(v), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 | 0x1 );
        unateness[3] = true;
      }

      /* 0-resub */
      if ( unateness[0] && unateness[3] )
      {
        return ( v << 1 );
      }
      if ( unateness[1] && unateness[2] )
      {
        return ( v << 1 ) + 1;
      }
      /* useless unate literal */
      if ( ( unateness[0] && unateness[2] ) || ( unateness[1] && unateness[3] ) )
      {
        pos_unate_lits.pop_back();
        neg_unate_lits.pop_back();
      }
      /* binate divisor */
      else if ( !unateness[0] && !unateness[1] && !unateness[2] && !unateness[3] )
      {
        binate_divs.emplace_back( v );
      }
    }
    return std::nullopt;
  }

  /* Sort the unate literals by the number of minterms in the intersection.
     - For `pos_unate_lits`, `on_off` = 1, sort by intersection with on-set;
     - For `neg_unate_lits`, `on_off` = 0, sort by intersection with off-set
   */
  void sort_unate_lits( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto& l : unate_lits )
    {
      l.score = kitty::count_ones( ( l.lit & 0x1 ? ~get_div( l.lit >> 1 ) : get_div( l.lit >> 1 ) ) & on_off_sets[on_off] );
    }
    std::sort( unate_lits.begin(), unate_lits.end(), [&]( unate_lit const& l1, unate_lit const& l2 ) {
        return l1.score > l2.score; // descending order
    });
  }

  void sort_unate_pairs( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto& p : unate_pairs )
    {
      if constexpr ( static_params::use_xor )
      {
        p.score = ( p.lit1 > p.lit2 ) ?
                    kitty::count_ones( ( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) )
                                     ^ ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) )
                                     & on_off_sets[on_off] )
                  : kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) )
                                     & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) )
                                     & on_off_sets[on_off] );
      }
      else
      {
        p.score = kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) )
                                   & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) )
                                   & on_off_sets[on_off] );
      }
    }
    std::sort( unate_pairs.begin(), unate_pairs.end(), [&]( fanin_pair const& p1, fanin_pair const& p2 ) {
        return p1.score > p2.score; // descending order
    });
  }

  /* See if there are two unate divisors covering all on-set bits or all off-set bits.
     - For `pos_unate_lits`, `on_off` = 1, try covering all on-set bits by combining two with an OR gate;
     - For `neg_unate_lits`, `on_off` = 0, try covering all off-set bits by combining two with an AND gate
   */
  std::optional<uint32_t> find_div_div( std::vector<unate_lit>& unate_lits, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      if ( unate_lits[i].score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_lits.size(); ++j )
      {
        uint32_t const& lit2 = unate_lits[j].lit;
        if ( unate_lits[i].score + unate_lits[j].score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        auto const ntt2 = lit2 & 0x1 ? get_div( lit2 >> 1 ) : ~get_div( lit2 >> 1 );
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          auto const new_lit = index_list.add_and( ( lit1 ^ 0x1 ), ( lit2 ^ 0x1 ) );
          return new_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_div_pair( std::vector<unate_lit>& unate_lits, std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_lits.size(); ++i )
    {
      uint32_t const& lit1 = unate_lits[i].lit;
      for ( auto j = 0u; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( unate_lits[i].score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        auto const ntt1 = lit1 & 0x1 ? get_div( lit1 >> 1 ) : ~get_div( lit1 >> 1 );
        TT ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                 ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                 | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
               | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }
        
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t new_lit1;
          if constexpr ( static_params::use_xor )
          {
            if ( pair2.lit1 > pair2.lit2 )
            {
              new_lit1 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            new_lit1 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          auto const new_lit2 = index_list.add_and( ( lit1 ^ 0x1 ), new_lit1 ^ 0x1 );
          return new_lit2 + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_pair_pair( std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_pairs.size(); ++i )
    {
      fanin_pair const& pair1 = unate_pairs[i];
      if ( pair1.score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_pairs.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs[j];
        if ( pair1.score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        TT ntt1, ntt2;
        if constexpr ( static_params::use_xor )
        {
          if ( pair1.lit1 > pair1.lit2 )
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) )
                 ^ ( pair1.lit2 & 0x1 ? ~get_div( pair1.lit2 >> 1 ) : get_div( pair1.lit2 >> 1 ) );
          }
          else
          {
            ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) )
                 | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          }
          if ( pair2.lit1 > pair2.lit2 )
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                 ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
          }
          else
          {
            ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                 | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
          }
        }
        else
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) )
               | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
               | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }

        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          uint32_t fanin_lit1, fanin_lit2;
          if constexpr ( static_params::use_xor )
          {
            if ( pair1.lit1 > pair1.lit2 )
            {
              fanin_lit1 = index_list.add_xor( pair1.lit1, pair1.lit2 );
            }
            else
            {
              fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            }
            if ( pair2.lit1 > pair2.lit2 )
            {
              fanin_lit2 = index_list.add_xor( pair2.lit1, pair2.lit2 );
            }
            else
            {
              fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
            }
          }
          else
          {
            fanin_lit1 = index_list.add_and( pair1.lit1, pair1.lit2 );
            fanin_lit2 = index_list.add_and( pair2.lit1, pair2.lit2 );
          }
          uint32_t const output_lit = index_list.add_and( fanin_lit1 ^ 0x1, fanin_lit2 ^ 0x1 );
          return output_lit + on_off;
        }
      }
    }
    return std::nullopt;
  }

  std::optional<uint32_t> find_xor()
  {
    /* collect XOR-type pairs (d1 ^ d2) & off = 0 or ~(d1 ^ d2) & on = 0, selecting d1, d2 from binate_divs */
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        auto const tt_xor = get_div( binate_divs[i] ) ^ get_div( binate_divs[j] );
        bool unateness[4] = {false, false, false, false};
        /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[0] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[1] = true;
        }

        /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
          unateness[2] = true;
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
          unateness[3] = true;
        }

        if ( unateness[0] && unateness[2] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ), ( binate_divs[j] << 1 ) );
        }
        if ( unateness[1] && unateness[3] )
        {
          return index_list.add_xor( ( binate_divs[i] << 1 ) + 1, ( binate_divs[j] << 1 ) );
        }
      }
    }

    return std::nullopt;
  }

  /* collect AND-type pairs (d1 & d2) & off = 0 or ~(d1 & d2) & on = 0, selecting d1, d2 from binate_divs */
  void collect_unate_pairs()
  {
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        collect_unate_pairs_detail<1, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<1, 0>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 0>( binate_divs[i], binate_divs[j] );
      }
    }
  }

  template<bool pol1, bool pol2>
  void collect_unate_pairs_detail( uint32_t div1, uint32_t div2 )
  {
    /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
    if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) )
    {
      pos_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)(!pol1), ( div2 << 1 ) + (uint32_t)(!pol2) );
    }
    /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
    else if ( kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[1] ) && !kitty::intersection_is_empty<TT, pol1, pol2>( get_div( div1 ), get_div( div2 ), on_off_sets[0] ) )
    {
      neg_unate_pairs.emplace_back( ( div1 << 1 ) + (uint32_t)(!pol1), ( div2 << 1 ) + (uint32_t)(!pol2) );
    }
  }

  inline TT const& get_div( uint32_t idx ) const
  {
    if constexpr ( static_params::copy_tts )
    {
      return divisors[idx];
    }
    else
    {
      return (*ptts)[divisors[idx]];
    }
  }

private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */

  const typename static_params::truth_table_storage_type* ptts;
  std::vector<std::conditional_t<static_params::copy_tts, TT, typename static_params::node_type>> divisors;

  index_list_t index_list;

  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<unate_lit> pos_unate_lits, neg_unate_lits;
  std::vector<uint32_t> binate_divs;
  std::vector<fanin_pair> pos_unate_pairs, neg_unate_pairs;

  stats& st;
}; /* xag_resyn_decompose */


struct xag_resyn_abc_stats
{};

template<class TT, class static_params = xag_resyn_static_params_default<TT>>
class xag_resyn_abc
{
public:
  using stats = xag_resyn_abc_stats;
  using index_list_t = large_xag_index_list;
  using truth_table_t = TT;

  explicit xag_resyn_abc( stats& st) noexcept
    : st( st ), counter( 0 )
  {
    static_assert( std::is_same_v<typename static_params::base_type, xag_resyn_static_params>, "Invalid static_params type" );
    static_assert( !static_params::preserve_depth && static_params::uniform_div_cost, "Advanced resynthesis is not implemented for this solver" );
  }

  virtual ~xag_resyn_abc()
  {
    abcresub::Abc_ResubPrepareManager( 0 );
    release();
  }

  template<class iterator_type, class truth_table_storage_type>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, iterator_type begin, iterator_type end, truth_table_storage_type const& tts, uint32_t max_size = std::numeric_limits<uint32_t>::max(), uint32_t max_level = std::numeric_limits<uint32_t>::max() )
  {
    (void)max_level;
    num_divisors = std::distance( begin, end ) + 2;
    num_blocks_per_truth_table = target.num_blocks();
    abcresub::Abc_ResubPrepareManager( num_blocks_per_truth_table );
    alloc();

    add_divisor( ~target & care ); /* off-set */
    add_divisor( target & care ); /* on-set */

    while ( begin != end )
    {
      add_divisor( tts[*begin] );
      ++begin;
    }

    return compute_function( max_size );
  }

protected:
  void add_divisor( TT const& tt )
  {
    assert( tt.num_blocks() == num_blocks_per_truth_table );
    for ( uint64_t i = 0ul; i < num_blocks_per_truth_table; ++i )
    {
      if constexpr ( std::is_same_v<TT, kitty::partial_truth_table> || std::is_same_v<TT, kitty::dynamic_truth_table> )
        Vec_WrdPush( abc_tts, tt._bits[i] );
      else // static_truth_table
        Vec_WrdPush( abc_tts, tt._bits );
    }
    Vec_PtrPush( abc_divs, Vec_WrdEntryP( abc_tts, counter * num_blocks_per_truth_table ) );
    ++counter;
  }

  std::optional<index_list_t> compute_function( uint32_t num_inserts )
  {
    int nLimit = num_inserts > std::numeric_limits<int>::max() ? std::numeric_limits<int>::max() : num_inserts;
    int * raw_list;
    int size = abcresub::Abc_ResubComputeFunction( 
               /* ppDivs */(void **)Vec_PtrArray( abc_divs ), 
               /* nDivs */Vec_PtrSize( abc_divs ), 
               /* nWords */num_blocks_per_truth_table, 
               /* nLimit */nLimit, 
               /* nDivsMax */static_params::max_binates, 
               /* iChoice */0, /* fUseXor */int(static_params::use_xor), /* fDebug */0, /* fVerbose */0, 
               /* ppArray */&raw_list );

    if ( size )
    {
      index_list_t xag_list;
      xag_list.add_inputs( num_divisors - 2 );
      for ( int i = 0; i < size - 1; i += 2 )
      {
        if ( raw_list[i] < raw_list[i+1] )
          xag_list.add_and( raw_list[i] - 2, raw_list[i+1] - 2 );
        else
          xag_list.add_xor( raw_list[i] - 2, raw_list[i+1] - 2 );
      }
      xag_list.add_output( raw_list[size - 1] < 2 ? raw_list[size - 1] : raw_list[size - 1] - 2 );
      return xag_list;
    }

    return std::nullopt;
  }

  void dump( std::string const file = "dump.txt" ) const
  {
    abcresub::Abc_ResubDumpProblem( file.c_str(), (void **)Vec_PtrArray( abc_divs ),  Vec_PtrSize( abc_divs ), num_blocks_per_truth_table );
  }

  void alloc()
  {
    assert( abc_tts == nullptr );
    assert( abc_divs == nullptr );
    abc_tts = abcresub::Vec_WrdAlloc( num_divisors * num_blocks_per_truth_table );
    abc_divs = abcresub::Vec_PtrAlloc( num_divisors );
  }

  void release()
  {
    assert( abc_divs != nullptr );
    assert( abc_tts != nullptr );
    Vec_PtrFree( abc_divs );
    Vec_WrdFree( abc_tts );
    abc_divs = nullptr;
    abc_tts = nullptr;
  }

protected:
  uint64_t num_divisors;
  uint64_t num_blocks_per_truth_table;
  uint64_t counter;

  abcresub::Vec_Wrd_t * abc_tts{nullptr};
  abcresub::Vec_Ptr_t * abc_divs{nullptr};

  stats& st;
}; /* xag_resyn_abc */

} /* namespace mockturtle */
