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
  \file search_core.hpp
  \brief the lower level searching core

  \author Hanyu Wang
*/

#pragma once

#include "../../utils/index_list.hpp"
#include "../../algorithms/cleanup.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>
#include <optional>
#include <queue>
#include <unordered_set>

namespace mockturtle::experimental
{

struct search_core_stats
{
  /* time for cost view the solution network */
  stopwatch<>::duration time_eval{0};

  /* time for searching the equivalent network */
  stopwatch<>::duration time_search{0};

  /* data */
  void report() const
  { 
    fmt::print( "[i]         <xag_resyn_decompose>\n" );
    fmt::print( "[i]             Evalutation      : {:>5.2f} secs\n", to_seconds( time_eval ) );
    fmt::print( "[i]             Searching        : {:>5.2f} secs\n", to_seconds( time_search ) );
  }
};

template<class Ntk, class TT>
class search_core
{
public:
  using stats = search_core_stats;
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;
  using cost_t = cost<typename Ntk::costfn_t>;
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
  inline TT const& get_div( uint32_t idx ) const
  {
    return (*ptts)[divisors[idx]];
  }

/**
 * Cost Evaluation Functions
 */
  uint32_t eval_result( Ntk & forest, index_list_t const& il )
  {
    Ntk tmp;
    decode( tmp, il );
    tmp.incr_trav_id();
    tmp.foreach_pi( [&]( auto const& n, auto i ){
      tmp.set_cost_val( n, div_costs[i] );
    } );
    uint32_t eval_1 = tmp.get_tmp_cost(); 

    uint32_t eval_2 = 0u;
    // insert il to forest
    insert( forest, std::begin( forest_leaves ), std::end( forest_leaves ), il, [&]( signal const& g ){
      forest.incr_trav_id();
      eval_2 = forest.get_cost( forest.get_node( g ), forest_leaves );
    } );

    // assert( forest.num_pos() == 1 );

    if( eval_1 != eval_2 )
    {
      std::cout << eval_1 << "->" << eval_2 << std::endl;
      std::cout << to_index_list_string( il ) << std::endl;
      assert(false);
    }

    return eval_1; // the cost of the whole network
  }
  bool update_result( Ntk & forest, std::optional<index_list_t> il )
  {
    if ( !il ) return false;
    uint32_t curr_cost = eval_result( forest, *il );
    if ( curr_cost < best_cost )
    {
      best_cost = curr_cost;
      index_list = *il;
      return true;
    }
    return false;
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
      p.score = ( p.lit1 > p.lit2 ) ?
                  kitty::count_ones( ( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) )
                                    ^ ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) ) )
                                    & on_off_sets[on_off] )
                : kitty::count_ones( ( p.lit1 & 0x1 ? ~get_div( p.lit1 >> 1 ) : get_div( p.lit1 >> 1 ) )
                                    & ( p.lit2 & 0x1 ? ~get_div( p.lit2 >> 1 ) : get_div( p.lit2 >> 1 ) )
                                    & on_off_sets[on_off] );
    }
    std::sort( unate_pairs.begin(), unate_pairs.end(), [&]( fanin_pair const& p1, fanin_pair const& p2 ) {
        return p1.score > p2.score; // descending order
    });
  }

  std::optional<index_list_t> find_and_detail( std::vector<unate_lit> const& unate_lits, uint32_t on_off )
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
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          auto const new_lit = il.add_and( ( lit1 ^ 0x1 ), ( lit2 ^ 0x1 ) );
          il.add_output( new_lit + on_off );
          return il;
        }
      }
    }
    return std::nullopt;
  }
  template<bool is_xor=false>
  std::optional<index_list_t> find_and_and_helper( std::vector<unate_lit>& unate_lits, std::vector<fanin_pair>& unate_pairs, uint32_t on_off )
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
        if constexpr ( is_xor )
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          uint32_t new_lit1;
          if constexpr ( is_xor )
          {
            new_lit1 = il.add_xor( pair2.lit1, pair2.lit2 );
          }
          else
          {
            new_lit1 = il.add_and( pair2.lit1, pair2.lit2 );
          }
          auto const new_lit2 = il.add_and( ( lit1 ^ 0x1 ), new_lit1 ^ 0x1 );
          il.add_output( new_lit2 + on_off );
          return il;
        }
      }
    }
    return std::nullopt;
  }
  template<bool left_xor=false, bool right_xor=false>
  std::optional<index_list_t> find_and_and_and_helper( std::vector<fanin_pair>& unate_pairs_1, std::vector<fanin_pair>& unate_pairs_2, uint32_t on_off )
  {
    for ( auto i = 0u; i < unate_pairs_1.size(); ++i )
    {
      fanin_pair const& pair1 = unate_pairs_1[i];
      if ( pair1.score * 2 < num_bits[on_off] )
      {
        break;
      }
      for ( auto j = i + 1; j < unate_pairs_2.size(); ++j )
      {
        fanin_pair const& pair2 = unate_pairs_2[j];
        if ( pair1.score + pair2.score < num_bits[on_off] )
        {
          break;
        }
        TT ntt1, ntt2;
        if constexpr ( left_xor )
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) )
               ^ ( pair1.lit2 & 0x1 ? ~get_div( pair1.lit2 >> 1 ) : get_div( pair1.lit2 >> 1 ) );
        }
        else
        {
          ntt1 = ( pair1.lit1 & 0x1 ? get_div( pair1.lit1 >> 1 ) : ~get_div( pair1.lit1 >> 1 ) )
               | ( pair1.lit2 & 0x1 ? get_div( pair1.lit2 >> 1 ) : ~get_div( pair1.lit2 >> 1 ) );
        }
        if constexpr ( right_xor )
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                ^ ( pair2.lit2 & 0x1 ? ~get_div( pair2.lit2 >> 1 ) : get_div( pair2.lit2 >> 1 ) );
        }
        else
        {
          ntt2 = ( pair2.lit1 & 0x1 ? get_div( pair2.lit1 >> 1 ) : ~get_div( pair2.lit1 >> 1 ) )
                | ( pair2.lit2 & 0x1 ? get_div( pair2.lit2 >> 1 ) : ~get_div( pair2.lit2 >> 1 ) );
        }
        if ( kitty::intersection_is_empty( ntt1, ntt2, on_off_sets[on_off] ) )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          uint32_t fanin_lit1, fanin_lit2;
          if constexpr ( left_xor )
          {
            fanin_lit1 = il.add_xor( pair1.lit1, pair1.lit2 );
          }
          else
          {
            fanin_lit1 = il.add_and( pair1.lit1, pair1.lit2 );
          }
          if constexpr ( right_xor )
          {
            fanin_lit2 = il.add_xor( pair2.lit1, pair2.lit2 );
          }
          else
          {
            fanin_lit2 = il.add_and( pair2.lit1, pair2.lit2 );
          }
          uint32_t const output_lit = il.add_and( fanin_lit1 ^ 0x1, fanin_lit2 ^ 0x1 );
          il.add_output( output_lit + on_off );
          return il;
        }
      }
    }
    return std::nullopt;
  }
  void prepare_clear()
  {
    pos_unate_lits.clear();
    neg_unate_lits.clear();
    binate_divs.clear();
    pos_unate_pairs.clear();
    neg_unate_pairs.clear();
    pos_unate_xor_pairs.clear();
    neg_unate_xor_pairs.clear();
    mem_xor.clear();
    mem_xor_xor.clear();
    mem_xor_and.clear();
    has_xor_xor = false;
    has_xor = false;
    has_xor_and = false;
    has_unateness = false;
    has_and_pairs = false;
    has_xor_pairs = false;
    has_lit_xor = false;
    has_init = false;

    index_list = std::nullopt;
  }
  void prepare_task()
  {
    assert( has_init == false );
    num_bits[0] = kitty::count_ones( on_off_sets[0] ); /* off-set */
    num_bits[1] = kitty::count_ones( on_off_sets[1] ); /* on-set */
    has_init = true;
  }
  void prepare_unateness()
  {
    assert( has_unateness == false && "already have unateness" );
    if ( has_init == false )
    {
      prepare_task();
    }
    for ( auto v = 1u; v < divisors.size(); ++v )
    {
      bool unateness[4] = {false, false, false, false};
      /* check intersection with off-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( v << 1 );
        unateness[0] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[0] ) )
      {
        pos_unate_lits.emplace_back( ( v << 1 ) | 0x1 );
        unateness[1] = true;
      }
      /* check intersection with on-set */
      if ( kitty::intersection_is_empty<TT, 1, 1>( get_div( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( v << 1 );
        unateness[2] = true;
      }
      else if ( kitty::intersection_is_empty<TT, 0, 1>( get_div( v ), on_off_sets[1] ) )
      {
        neg_unate_lits.emplace_back( ( v << 1 ) | 0x1 );
        unateness[3] = true;
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
    sort_unate_lits( pos_unate_lits, 1 );
    sort_unate_lits( neg_unate_lits, 0 );
    has_unateness = true;
  }
  void prepare_xor()
  {
    assert( has_xor == false );
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      mem_xor[ on_off_sets[1] ^ get_div( i ) ] = i;
    }
    has_xor = true;
  }
  void prepare_xor_xor()
  {
    assert( has_xor_xor == false );
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        mem_xor_xor[ get_div( i ) ^ get_div( j ) ^ on_off_sets[1] ] = i * divisors.size() + j;
      }
    }
    has_xor_xor = true;
  }
  void prepare_xor_and()
  {
    assert( has_xor_and == false );
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        for ( auto on_off_1 = 0u; on_off_1 < 2; on_off_1++ )
        {
          for ( auto on_off_2 = 0u; on_off_2 < 2; on_off_2++ )
          {
            auto const tt = ( on_off_1? ~get_div( i ) : get_div( i ) ) & ( on_off_2? ~get_div( j ) : get_div( j ) );
            mem_xor_and[ tt ^ on_off_sets[1] ] = ( ( i << 1 ) + on_off_1 ) * 2 * divisors.size() + ( ( j << 1 ) + on_off_2 );
          }
        }
      }
    }
    has_xor_and = true;
  }
  void prepare_and_pairs()
  {
    if ( has_unateness == false )
    {
      prepare_unateness();
    }
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        collect_unate_pairs_detail<1, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 1>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<1, 0>( binate_divs[i], binate_divs[j] );
        collect_unate_pairs_detail<0, 0>( binate_divs[i], binate_divs[j] );
      }
    };
    sort_unate_pairs( pos_unate_pairs, 1 );
    sort_unate_pairs( neg_unate_pairs, 0 );
    has_and_pairs = true;
  }
  void prepare_xor_pairs()
  {
    if ( has_unateness == false )
    {
      prepare_unateness();
    }
    for ( auto i = 0u; i < binate_divs.size(); ++i )
    {
      for ( auto j = i + 1; j < binate_divs.size(); ++j )
      {
        auto const tt_xor = get_div( binate_divs[i] ) ^ get_div( binate_divs[j] );
        /* check intersection with off-set; additionally check intersection with on-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_xor_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) )
        {
          pos_unate_xor_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
        }
        /* check intersection with on-set; additionally check intersection with off-set is not empty (otherwise it's useless) */
        if ( kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 1, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_xor_pairs.emplace_back( binate_divs[i] << 1, binate_divs[j] << 1, true );
        }
        if ( kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[1] ) && !kitty::intersection_is_empty<TT, 0, 1>( tt_xor, on_off_sets[0] ) )
        {
          neg_unate_xor_pairs.emplace_back( ( binate_divs[i] << 1 ) + 1, binate_divs[j] << 1, true );
        }
      }
    }
    sort_unate_pairs( pos_unate_xor_pairs, 1 );
    sort_unate_pairs( neg_unate_xor_pairs, 0 );
    has_xor_pairs = true;
  }

  std::optional<index_list_t> find_wire()
  {
    if ( has_init == false )
    {
      prepare_task();
    }
    if ( num_bits[0] == 0 )
    {
      index_list_t il;
      il.clear();
      il.add_inputs( divisors.size() - 1 );
      il.add_output( 1 );
      return il;
    }
    if ( num_bits[1] == 0 )
    {
      index_list_t il;
      il.clear();
      il.add_inputs( divisors.size() - 1 );
      il.add_output( 0 );
      return il;
    }
    for ( auto v = 1u; v < divisors.size(); ++v )
    {
      if ( get_div( v ) == on_off_sets[1] )
      {
        index_list_t il;
        il.clear();
        il.add_inputs( divisors.size() - 1 );
        il.add_output( v << 1 );
        return il;
      }
      if ( get_div( v ) == on_off_sets[0] )
      {
        index_list_t il;
        il.clear();
        il.add_inputs( divisors.size() - 1 );
        il.add_output( ( v << 1 ) + 1 );
        return il;
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_and()
  {
    if ( has_unateness == false )
    {
      prepare_unateness();
    }
    return find_and_detail( neg_unate_lits, 0 );
  }
  std::optional<index_list_t> find_or()
  {
    if ( has_unateness == false )
    {
      prepare_unateness();
    }
    return find_and_detail( pos_unate_lits, 1 );
  }
  std::optional<index_list_t> find_xor()
  {
    if ( has_xor == false )
    {
      prepare_xor();
    }
    for ( auto i = 1u; i < divisors.size(); ++i )
    {
      auto const tt = get_div( i );
      if ( mem_xor.find( tt ) != mem_xor.end() )
      {
        index_list_t il;
        il.clear();
        il.add_inputs( divisors.size() -1 );
        il.add_output( il.add_xor( ( i << 1 ), mem_xor[ tt ] << 1 ) );
        return il;
      }
      if ( mem_xor.find( ~tt ) != mem_xor.end() )
      {
        index_list_t il;
        il.clear();
        il.add_inputs( divisors.size() -1 );
        il.add_output( il.add_xor( ( i << 1 ) + 1, mem_xor[ ~tt ] << 1 ) );
        return il;
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_or_and()
  {
    if ( has_and_pairs == false )
    {
      prepare_and_pairs();
    }
    return find_and_and_helper( pos_unate_lits, pos_unate_pairs, 1 );
  }
  std::optional<index_list_t> find_and_and()
  {
    if ( has_and_pairs == false )
    {
      prepare_and_pairs();
    }
    return find_and_and_helper( neg_unate_lits, neg_unate_pairs, 0 );
  }
  std::optional<index_list_t> find_and_xor()
  {
    if ( has_xor_pairs == false )
    {
      prepare_xor_pairs();
    }
    auto ret = find_and_and_helper<true>( pos_unate_lits, pos_unate_xor_pairs, 1 );
    if ( !ret ) ret = find_and_and_helper<true>( neg_unate_lits, neg_unate_xor_pairs, 0 );
    return ret;
  }
  std::optional<index_list_t> find_xor_xor()
  {
    if ( has_xor == false )
    {
      prepare_xor();
    }
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        auto const tt = get_div( i ) ^ get_div( j );
        if ( mem_xor.find( tt ) != mem_xor.end() )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          auto new_lit1 = il.add_xor( (i << 1), (j << 1) );
          auto new_lit2 = il.add_xor( new_lit1, mem_xor[ tt ] << 1 );
          il.add_output( new_lit2 );
          return il;
        }
        if ( mem_xor.find( ~tt ) != mem_xor.end() )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          auto new_lit1 = il.add_xor( (i << 1), (j << 1) );
          auto new_lit2 = il.add_xor( new_lit1 ^ 0x1, mem_xor[ ~tt ] << 1 );
          il.add_output( new_lit2 );
          return il;
        }
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_xor_xor_xor()
  {
    if ( has_xor_xor == false )
    {
      prepare_xor_xor();
    }
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        auto const tt = get_div( i ) ^ get_div( j );
        if ( mem_xor_xor.find( tt ) != mem_xor_xor.end() )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          auto new_lit1 = il.add_xor( (i << 1), (j << 1) );
          auto new_lit2 = il.add_xor( ( mem_xor_xor[ tt ] % divisors.size() ) << 1, ( mem_xor_xor[ tt ] / divisors.size() ) << 1 );
          auto new_lit3 = il.add_xor( new_lit2, new_lit1 );
          il.add_output( new_lit3 );
          return il;
        }
        if ( mem_xor_xor.find( ~tt ) != mem_xor_xor.end() )
        {
          index_list_t il;
          il.clear();
          il.add_inputs( divisors.size() - 1 );
          auto new_lit1 = il.add_xor( (i << 1), (j << 1) );
          auto new_lit2 = il.add_xor( ( mem_xor_xor[ ~tt ] % divisors.size() ) << 1, ( mem_xor_xor[ ~tt ] / divisors.size() ) << 1 );
          auto new_lit3 = il.add_xor( new_lit2 ^ 0x1, new_lit1 );
          il.add_output( new_lit3 );
          return il;
        }
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_xor_xor_and()
  {
    if ( has_xor_xor == false )
    {
      prepare_xor_xor();
    }
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        for ( auto on_off_1 = 0u; on_off_1 < 2; on_off_1++ )
        {
          for ( auto on_off_2 = 0u; on_off_2 < 2; on_off_2++ )
          {
            auto const tt = ( on_off_1? ~get_div( i ) : get_div( i ) ) & ( on_off_2? ~get_div( j ) : get_div( j ) );
            if ( mem_xor_xor.find( tt ) != mem_xor_xor.end() )
            {
              index_list_t il;
              il.clear();
              il.add_inputs( divisors.size() - 1 );
              auto new_lit1 = il.add_and( (i << 1) + on_off_1, (j << 1) + on_off_2);
              auto new_lit2 = il.add_xor( ( mem_xor_xor[ tt ] % divisors.size() ) << 1, ( mem_xor_xor[ tt ] / divisors.size() ) << 1 );
              auto new_lit3 = il.add_xor( new_lit2, new_lit1 );
              il.add_output( new_lit3 );
              return il;
            }
            if ( mem_xor_xor.find( ~tt ) != mem_xor_xor.end() )
            {
              index_list_t il;
              il.clear();
              il.add_inputs( divisors.size() - 1 );
              auto new_lit1 = il.add_and( (i << 1) + on_off_1, (j << 1) + on_off_2);
              auto new_lit2 = il.add_xor( ( mem_xor_xor[ ~tt ] % divisors.size() ) << 1, ( mem_xor_xor[ ~tt ] / divisors.size() ) << 1 );
              auto new_lit3 = il.add_xor( new_lit2 ^ 0x1, new_lit1 );
              il.add_output( new_lit3 );
              return il;
            }
          }
        }
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_xor_and_and()
  {
    if ( has_xor_and == false )
    {
      prepare_xor_and();
    }
    for ( auto i = 1u; i < divisors.size(); i++ )
    {
      for ( auto j = i+1; j < divisors.size(); j++ )
      {
        for ( auto on_off_1 = 0u; on_off_1 < 2; on_off_1++ )
        {
          for ( auto on_off_2 = 0u; on_off_2 < 2; on_off_2++ )
          {
            auto const tt = ( on_off_1? ~get_div( i ) : get_div( i ) ) & ( on_off_2? ~get_div( j ) : get_div( j ) );
            if ( mem_xor_and.find( tt ) != mem_xor_and.end() )
            {
              index_list_t il;
              il.clear();
              il.add_inputs( divisors.size() - 1 );
              auto new_lit1 = il.add_and( (i << 1) + on_off_1, (j << 1) + on_off_2 );
              auto new_lit2 = il.add_and( mem_xor_and[ tt ] % ( 2 * divisors.size() ), mem_xor_and[ tt ] / ( 2 * divisors.size() ) );
              auto new_lit3 = il.add_xor( new_lit1, new_lit2 );
              il.add_output( new_lit3 );
              return il;
            }
            if ( mem_xor_and.find( ~tt ) != mem_xor_and.end() )
            {
              index_list_t il;
              il.clear();
              il.add_inputs( divisors.size() - 1 );
              auto new_lit1 = il.add_and( (i << 1) + on_off_1, (j << 1) + on_off_2 );
              auto new_lit2 = il.add_and( mem_xor_and[ ~tt ] % ( 2 * divisors.size() ), mem_xor_and[ ~tt ] / ( 2 * divisors.size() ) );
              auto new_lit3 = il.add_xor( new_lit1 ^ 0x1, new_lit2 );
              il.add_output( new_lit3 );
              return il;
            }
          }
        }
      }
    }
    return std::nullopt;
  }
  std::optional<index_list_t> find_and_and_and()
  {
    if ( has_and_pairs == false )
    {
      prepare_and_pairs();
    }
    auto ret = find_and_and_and_helper( pos_unate_pairs, pos_unate_pairs, 1 );
    if ( !ret ) ret = find_and_and_and_helper( neg_unate_pairs, neg_unate_pairs, 0 );
    return ret;
  }
  std::optional<index_list_t> find_and_and_xor()
  {
    if ( has_and_pairs == false )
    {
      prepare_and_pairs();
    }
    if ( has_xor_pairs == false )
    {
      prepare_xor_pairs();
    }
    auto ret = find_and_and_and_helper<true>( pos_unate_xor_pairs, pos_unate_pairs,  1 );
    if ( !ret ) ret = find_and_and_and_helper<true>( neg_unate_xor_pairs, neg_unate_pairs,  0 );
    return ret;
  }
  std::optional<index_list_t> find_and_xor_xor()
  {
    if ( has_xor_pairs == false )
    {
      prepare_xor_pairs();
    }
    auto ret = find_and_and_and_helper<true, true>( pos_unate_xor_pairs, pos_unate_xor_pairs, 1 );
    if ( !ret ) ret = find_and_and_and_helper<true, true>( neg_unate_xor_pairs, neg_unate_xor_pairs, 0 );
    return ret;
  }
public:
  explicit search_core( Ntk const& ntk, stats& st ) noexcept
    : ntk( ntk ), st( st )
  {
    divisors.reserve( 200u );
  }
  template<class iterator_type, class truth_table_storage_type>
  std::optional<index_list_t> operator()( TT const& target, TT const& care, std::vector<signal> const& divs, iterator_type begin, iterator_type end, truth_table_storage_type const& tts, uint32_t max_cost = std::numeric_limits<uint32_t>::max() )
  {
    ptts = &tts;
    on_off_sets[0] = ~target & care;
    on_off_sets[1] = target & care;

    divisors.resize( 1 ); /* clear previous data and reserve 1 dummy node for constant */
    while ( begin != end )
    {
      divisors.emplace_back( *begin );
      ++begin;
    }

    // prepare solution forest
    Ntk forest; // create an empty network
    forest_leaves.clear();
    forest_root = std::nullopt;
    div_costs.clear();
    
    for( signal div : divs )
    {
      signal const& s = forest.create_pi(); 
      node n = forest.get_node( s );
      forest_leaves.emplace_back( s );
      cost_t div_cost = ntk.get_cost_val( ntk.get_node( div ) );
      forest.set_cost_val( n, div_cost );
      div_costs.emplace_back( div_cost );
    }

    best_cost = max_cost; // upper bound
    prepare_clear();

    bool done = false;

    std::optional<index_list_t> il = find_wire(); if ( il ) return( *il ); /* searching constant is useless and cause error */
    il = call_with_stopwatch( st.time_search, [&](){ return find_or()          ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and()         ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_xor()         ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and_and()     ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_or_and()      ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and_xor()     ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_xor_xor()     ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_xor_xor_xor() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_xor_xor_and() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_xor_and_and() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and_and_and() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and_and_xor() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;
    il = call_with_stopwatch( st.time_search, [&](){ return find_and_xor_xor() ;}   ); done = call_with_stopwatch( st.time_eval, [&](){ return update_result( forest, il ); } );// if (done) return index_list;

    return index_list;
  }
private:
  std::array<TT, 2> on_off_sets;
  std::array<uint32_t, 2> num_bits; /* number of bits in on-set and off-set */

  const std::vector<TT>* ptts;
  std::vector<uint32_t> divisors;
  std::vector<cost_t> div_costs;
  std::array<TT, 0x100u> tts_xors;
  std::unordered_map<TT, uint32_t, kitty::hash<TT>> mem_xor;
  std::unordered_map<TT, uint32_t, kitty::hash<TT>> mem_xor_xor;
  std::unordered_map<TT, uint32_t, kitty::hash<TT>> mem_xor_and;
  bool has_xor;
  bool has_xor_xor;
  bool has_xor_and;
  bool has_unateness;
  bool has_and_pairs;
  bool has_xor_pairs;
  bool has_init;
  bool has_lit_xor;
  /* positive unate: not overlapping with off-set
     negative unate: not overlapping with on-set */
  std::vector<unate_lit> pos_unate_lits, neg_unate_lits;
  std::vector<uint32_t> binate_divs;
  std::vector<fanin_pair> pos_unate_pairs, neg_unate_pairs;
  std::vector<fanin_pair> pos_unate_xor_pairs, neg_unate_xor_pairs;

  Ntk const& ntk;
  std::vector<signal> forest_leaves;
  std::optional<signal> forest_root;
  stats& st;

  std::optional<index_list_t> index_list;
  uint32_t best_cost;
};

}