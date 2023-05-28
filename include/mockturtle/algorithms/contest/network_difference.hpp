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
  \file network difference.hpp
  \brief function to calculate the distance between two networks

  \author Hanyu WANG
*/

#pragma once

#include "../../io/write_patterns.hpp"
#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include "../../traits.hpp"
#include "../../utils/index_list.hpp"
#include "../../views/depth_view.hpp"
#include "../../views/fanout_view.hpp"
#include "../circuit_validator.hpp"
#include "../detail/resub_utils.hpp"
#include "../dont_cares.hpp"
#include "../pattern_generation.hpp"
#include "../resyn_engines/aig_enumerative.hpp"
#include "../resyn_engines/mig_enumerative.hpp"
#include "../resyn_engines/mig_resyn.hpp"
#include "../resyn_engines/xag_resyn.hpp"
#include "../simulation.hpp"
#include <kitty/kitty.hpp>

#include <functional>
#include <optional>
#include <vector>

namespace mockturtle::experimental
{
struct network_difference_params
{
  double decay{ 0 };

  /* use validator to verify if the function is exactly the same */
  bool use_validator{ false };

  /* with same simulation pattern but structurally different */
  bool structural_penalty{ false };

  /* maximum number of patterns */
  uint32_t max_patterns{ 1024 };

  /* verbose */
  bool verbose;
};

namespace detail
{

template<class Ntk, class TT>
class network_difference_impl
{

public:
  using signal = typename Ntk::signal;
  using node = typename Ntk::node;

public:
  network_difference_impl( Ntk const& basis, network_difference_params const& ps )
      : basis( basis ), ps( ps ), sim( basis.num_pis() ), tts( basis )
  {
    /* precomputations: prepare the simulation patterns */
    simulate_nodes<TT>( basis, tts, sim );
    basis.foreach_gate( [&]( auto n ) {
      basis_tts[tts[n]] = n;
      basis_tts[~tts[n]] = n;
    } );
  }

  /* return the difference */
  double run( Ntk const& ntk )
  {
    default_simulator<TT> _sim( ntk.num_pis() );
    unordered_node_map<TT, Ntk> _tts( ntk );
    simulate_nodes<TT>( ntk, _tts, _sim );

    double difference = 0.0;
    ntk.foreach_gate( [&]( auto n ) {
      if ( basis_tts.find( _tts[n] ) == basis_tts.end() )
        difference += 1.0;
    } );
    return difference;
  }

private:
  Ntk const& basis;
  default_simulator<TT> sim;
  unordered_node_map<TT, Ntk> tts;
  network_difference_params const& ps;
  std::unordered_map<TT, node, kitty::hash<TT>> basis_tts;
};
} // namespace detail

/**
 * @brief Measure the difference between two logic networks, or more specifically, the distance from 
 * ntk to basis
 * 
 */
template<class Ntk, class NtkSrc>
double network_difference( Ntk const& ntk, NtkSrc const& basis, network_difference_params const& ps = {} )
{
  static_assert( std::is_same_v<typename Ntk::base_type, typename NtkSrc::base_type> && "network and basis should of the same type" );
  assert( ntk.num_pis() == basis.num_pis() && "network inputs not matched" );
  using TT = kitty::dynamic_truth_table;
  detail::network_difference_impl<Ntk, TT> p( basis, ps );
  const auto ret = p.run( ntk );
  assert( ret <= ntk.num_gates() - ntk.num_pos() );
  return ret;
}

} // namespace mockturtle::experimental