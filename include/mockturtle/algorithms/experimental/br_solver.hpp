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
  \file br_solver.hpp
  \brief Boolean relation solver

  \author Hanyu Wang
*/

#pragma once

#include "../../algorithms/simulation.hpp"
#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include "../../traits.hpp"
#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../resyn_engines/xag_resyn.hpp"
#include <kitty/kitty.hpp>

#include <functional>
#include <optional>
#include <vector>

namespace mockturtle::experimental
{

struct br_solver_params
{
};

struct br_solver_stats
{

  void report() const
  {
    // clang-format off
    fmt::print( "[i] Boolean solver report\n" );
    // clang-format on
  }
};

namespace detail
{

} /* namespace detail */

class br_solver
{
public:
  using index_list_t = large_xag_index_list;
  using TT = kitty::partial_truth_table;
  using stats = br_solver_stats;
  using Ntk = aig_network;
  using node = Ntk::node;
  using signal = Ntk::signal;

private:
  enum lit_type
  {
    EQUAL,
    EQUAL_INV,
    POS_UNATE,
    NEG_UNATE,
    POS_UNATE_INV,
    NEG_UNATE_INV,
    BINATE,
    DONT_CARE
  };

private:
  void print_target()
  {
    for ( uint32_t idx = 0; idx < targets.size(); ++idx )
    {
      kitty::print_binary( targets[idx] & ~mask[idx] );
      std::cout << "\n";
    }
  }

  void update_projection( uint32_t tid )
  {
    /** TODO: get the projection
     * for each assignment of target, for each bit with 1. 
    */
    fmt::print( "[i] projection of output {}\n", tid );
    for ( uint64_t offset = 0; offset < tt_size; ++offset )
    {
      std::optional<bool> temp_val;
      for ( uint32_t i = 0; i < targets.size(); ++i )
      {
        if ( kitty::get_bit( targets[i] & ~mask[i], offset ) == 1 )
        {
          bool val = ( i >> tid ) & 0b01;
          if ( temp_val && *temp_val != val )
          {
            // mark as dont care
            kitty::clear_bit( care_out[tid], offset );
            temp_val = std::nullopt;
            break;
          }
          temp_val = val;
        }
      }
      if ( temp_val )
      {
        kitty::set_bit( care_out[tid], offset );
        if ( *temp_val == 1 )
        {
          kitty::set_bit( tt_out[tid], offset );
        }
        if ( *temp_val == 0 )
        {
          kitty::clear_bit( tt_out[tid], offset );
        }
      }
    }
    std::cout << "TT = ";
    kitty::print_binary( tt_out[tid] );
    std::cout << std::endl;
    std::cout << "Care = ";
    kitty::print_binary( care_out[tid] );
    std::cout << std::endl;
  }

  void propagate_and_mask( uint32_t oid )
  {
    incomplete_node_map<TT, Ntk> tts( ntk );
    simulate_nodes<Ntk>( ntk, tts, sim, true );
    signal s = ntk.po_at( oid );
    node n = ntk.get_node( s );
    const TT& tt = tts[n];
    fmt::print( "[i] propagate PO {} = ", oid );
    kitty::print_binary( tt );
    std::cout << std::endl;
    for ( uint32_t offset = 0; offset < tt_size; ++offset )
    {
      uint32_t val = kitty::get_bit( tt, offset ) & 0b01;
      for ( uint32_t tid = 0; tid < targets.size(); ++tid )
      {
        uint32_t ref_val = ( tid >> oid ) & 0b01;
        if ( ref_val != val ) // mask if is different / conflict
        {
          kitty::set_bit( mask[tid], offset );
        }
      }
    }
  }

  std::optional<uint32_t> propagate_and_verify()
  {
    incomplete_node_map<TT, Ntk> tts( ntk );
    simulate_nodes<Ntk>( ntk, tts, sim, true );

    for ( uint32_t offset = 0; offset < tt_size; ++offset )
    {
      uint32_t tid = 0;
      ntk.foreach_po( [&]( auto const& n, auto i ) {
        uint32_t val = kitty::get_bit( tts[n], offset ) & 0b01;
        tid &= ( val << i );
      } );
      // check if edge in the adjacent matrix exists
      assert( tid < targets.size() );
      if ( kitty::get_bit( targets[tid] & ~mask[tid], offset ) == 0 )
      {
        return offset;
      }
    }
    // no violation found
    return std::nullopt;
  }

  void solve_single_output( uint32_t tid )
  {
    xag_resyn_stats st;
    xag_resyn_decompose<TT, xag_resyn_static_params_for_sim_resub<Ntk>> engine( st );

    incomplete_node_map<TT, Ntk> tts( ntk )
        simulate_nodes<Ntk>( ntk, tts, sim, false );
    // TODO: collect more divisors

    fmt::print( "[i] solving output {} \n", tid );
    for ( auto div : divs )
    {
      kitty::print_binary( tts[div] );
      std::cout << std::endl;
    }
    const std::optional<index_list_t> res = engine( tt_out[tid], care_out[tid], std::begin( divs ), std::end( divs ), tts );
    assert( res );
    insert<false>( ntk, divs.begin(), divs.end(), *res, [&]( signal const& g ) {
      ntk.create_po( g );
    } );
    fmt::print( "[i] network has {} gates.\n", ntk.num_gates() );
  }

  std::optional<index_list_t> br_naive()
  {
    index_list_t res;
    for ( uint32_t tid = 0; tid < num_target; ++tid )
    {
      print_target();
      update_projection( tid );
      solve_single_output( tid );
      propagate_and_mask( tid );
    }
    encode( res, ntk );
    return res;
  }

public:
  br_solver()
  {
  }
  std::optional<index_list_t> operator()( std::vector<TT> const& _divs, std::vector<TT> const& _targets, uint32_t _max_cost = std::numeric_limits<uint32_t>::max() )
  {
    tt_size = _divs[0].num_bits();
    max_cost = _max_cost;
    num_target = std::log2( _targets.size() );
    num_divs = _divs.size();

    for ( uint32_t i = 0; i < _divs.size(); ++i )
    {
      tt_divs.emplace_back( _divs[i] );
      signal s = ntk.create_pi();
      divs.emplace_back( ntk.get_node( s ) );
    }
    for ( uint32_t i = 0; i < _targets.size(); ++i )
    {
      targets.emplace_back( _targets[i] );
      mask.emplace_back( TT( tt_size ) );
      assert( targets[i].num_bits() == mask[i].num_bits() );
    }
    for ( uint32_t tid = 0; tid < num_target; ++tid )
    {
      tt_out.emplace_back( TT( tt_size ) );
      care_out.emplace_back( TT( tt_size ) );
    }

    // initialize simulator
    sim = partial_simulator( tt_divs );

    return br_naive();
  }

private:
  uint32_t max_cost;
  uint32_t num_target;
  uint32_t num_divs;
  uint32_t tt_size;
  std::vector<TT> targets;
  std::vector<TT> tt_out;
  std::vector<TT> care_out;
  std::vector<TT> tt_divs;
  std::vector<TT> mask;
  std::vector<node> divs;
  partial_simulator sim; // to simulate pattern
  Ntk ntk;
};

} // namespace mockturtle::experimental