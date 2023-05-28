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
  \file evolutionary_algorithm.hpp

  \author Hanyu WANG
*/

#pragma once

#include "../balancing/sop_balancing.hpp"
#include "../cleanup.hpp"
#include "../cut_rewriting.hpp"
#include "../node_resynthesis/xag_npn.hpp"
#include "../../utils/network_utils.hpp"

#include <algorithm>
#include <random>

namespace mockturtle
{
struct evolutionary_algorithm_params
{
  bool verbose{ true };

  uint32_t size_limit{ 200u };
  uint32_t num_parents{ 3 };
  uint32_t num_offsprings{ 5 };
  uint32_t num_generations{ 10 };

  double mutation_rate{ 0.05 };
};

namespace detail
{
template<class Ntk, class ResynFn, class MutFn>
class evolutionary_algorithm_impl
{
public:
  using TT = kitty::dynamic_truth_table;
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;

public:
  evolutionary_algorithm_impl( Ntk const& ntk, ResynFn const& fn, MutFn const& _fn, evolutionary_algorithm_params const& ps = {} )
      : ps( ps ), fn( fn ), _fn( _fn )
  {
    ntks.emplace_back( ntk.clone() );
    sort_chromosomes();
    while ( ntks.size() < ps.num_parents + ps.num_offsprings )
    {
      ntks.emplace_back( Ntk() ); /* place-holders */
    }
  }

  template<typename NtkIterator>
  evolutionary_algorithm_impl( NtkIterator begin, NtkIterator end, ResynFn const& fn, MutFn const& _fn, evolutionary_algorithm_params const& ps = {} )
      : ps( ps ), fn( fn ), _fn( _fn )
  {
    while ( begin != end && ntks.size() < ps.num_parents + ps.num_offsprings )
    {
      Ntk ntk = ( begin++ )->clone();
      fn( ntk );
      ntks.emplace_back( ntk );
    }
    sort_chromosomes();
    while ( ntks.size() < ps.num_parents + ps.num_offsprings )
    {
      ntks.emplace_back( Ntk() ); /* place-holders */
    }
  }

  Ntk run()
  {
    if ( ntks[0].num_pos() == 1 ) return ntks[0]; // skip single output problem
    if ( ntks[0].num_gates() >= ps.size_limit ) return ntks[0]; // skip large problem
    for ( int i = 0; i < ps.num_generations; i++ )
    {
      evolve( i );
      // fmt::print( "{},{}\n", i, ntks[0].num_gates() );
    }
    return ntks[0];
  }

  void sort_chromosomes()
  {
    // the fitness is defined as the number of gates in the network
    std::sort( ntks.begin(), ntks.end(), [&]( Ntk const& ntk1, Ntk const& ntk2 ) { 
      return ntk1.num_gates() < ntk2.num_gates(); } );
  }

  Ntk crossover( std::vector<uint32_t> const& parents_order )
  {
    Ntk ntk;
    std::vector<signal> pis( ntks[0].num_pis() );
    std::generate( std::begin( pis ), std::end( pis ), [&]() { return ntk.create_pi(); } );
    std::vector<std::vector<signal>> pos;
    for ( uint32_t i = 0; i < ps.num_parents; i++ )
    {
      pos.emplace_back( cleanup_dangling( ntks[i], ntk, std::begin( pis ), std::end( pis ) ) );
    }
    for ( uint32_t i = 0; i < ntks[0].num_pos(); i++ )
    {
      ntk.create_po( pos[parents_order[i]][i] );
    }
    return cleanup_dangling( ntk );
  }

  void evolve( uint32_t age )
  {
    std::vector<uint32_t> parents_order( ntks[0].num_pos() );
    for ( int i = ps.num_offsprings + ps.num_parents - 1; i >= 0; i-- )
    {
      if ( i >= ps.num_parents ) // crossover
      {
        std::generate( std::begin( parents_order ), std::end( parents_order ), [&]() { return rand() % ps.num_parents; } );
        ntks[i] = crossover( parents_order );
        fn( ntks[i] );
      }
      else // mutation
      {
        if ( (double)rand() / RAND_MAX < ps.mutation_rate )
        {
          _fn( ntks[i] );
          fn( ntks[i] );
        }
      }
    }
    sort_chromosomes();
  }

private:
  std::vector<Ntk> ntks;
  ResynFn const& fn;
  MutFn const& _fn;
  evolutionary_algorithm_params const& ps;
};
} // namespace detail

/**
 * @brief Use evolutionary algorithm to optimize a network
 * 
 */
template<class Ntk, class ResynFn, class MutFn>
Ntk evolutionary_algorithm( Ntk const& ntk, ResynFn const& fn, MutFn const& _fn, evolutionary_algorithm_params const& ps = {} )
{
  detail::evolutionary_algorithm_impl<Ntk, ResynFn, MutFn> p( ntk, fn, _fn, ps );
  const auto ret = p.run();
  return ret;
}

/**
 * \brief Use evolutionary algorithm on first generation
 * 
 */
template<typename Ntk, typename NtkIterator, class ResynFn, class MutFn>
Ntk evolutionary_algorithm( NtkIterator begin, NtkIterator end, ResynFn const& fn, MutFn const& _fn, evolutionary_algorithm_params const& ps = {} )
{
  detail::evolutionary_algorithm_impl<Ntk, ResynFn, MutFn> p( begin, end, fn, _fn, ps );
  const auto ret = p.run();
  return ret;
}

} // namespace mockturtle