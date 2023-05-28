/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file permute_outputs.hpp
  \brief Randomizes the topological ordering of a network

  \author Alessandro Tempia Calvino
  \author Andrea Costamagna
*/

#pragma once

#include <vector>
#include <random>
#include <algorithm>

#include "../traits.hpp"
#include "../utils/node_map.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class permute_outputs_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  permute_outputs_impl( Ntk& ntk, std::vector<uint64_t> const& order )
      : ntk( ntk ), order( order )
  {
  }

  Ntk run()
  {
    Ntk dest;
    node_map<signal, Ntk> old2new( ntk );

    ntk.incr_trav_id();
    ntk.incr_trav_id();

    old2new[ntk.get_constant( false )] = dest.get_constant( false );
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );

    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = dest.get_constant( true );
      ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
    }

    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = dest.create_pi();
      ntk.set_visited( n, ntk.trav_id() );
    } );

    std::vector<signal> outputs;
    outputs.reserve( ntk.num_pos() );

    ntk.foreach_po( [&]( auto const& f ) {
      outputs.push_back( f );
    } );

    for( uint32_t i = 0; i < order.size(); ++i )
    {
      signal q = outputs[order[i]];

      topo_rec( dest, old2new, ntk.get_node( q ));

      if ( ntk.is_complemented( q ) )
        dest.create_po( dest.create_not( old2new[ q ] ) );
      else
        dest.create_po( old2new[ q ] );
    }

    return dest;
  }

private:
  void topo_rec( Ntk& dest, node_map<signal, Ntk>& old2new, node const& n )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    assert( ntk.visited( n ) != ntk.trav_id() - 1 );

    ntk.set_visited( n, ntk.trav_id() - 1 );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      topo_rec( dest, old2new, ntk.get_node( f ) );
    } );

    std::vector<signal> children;
    children.reserve( ntk.fanin_size( n ) );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
        children.push_back( dest.create_not( old2new[f] ) );
      else
        children.push_back( old2new[f] );
    } );

    old2new[n] = dest.clone_node( ntk, n, children );

    ntk.set_visited( n, ntk.trav_id() );
  }

private:
  Ntk& ntk;
  std::vector<uint64_t> const& order;
};

} /* namespace detail */

/*! \brief Outputs permutation.
 *
 * This method performs the outputs permutation for a given outputs order
 *
   \verbatim embed:rst

   .. note::

      This method returns a network with the desired outputs order.  It does
      *not* modify the input network.
   \endverbatim
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_node`
 * - `foreach_pi`
 * - `foreach_po`
 * - `clone_node`
 * - `is_pi`
 * - `is_constant`
 */
template<class Ntk>
Ntk permute_outputs( Ntk& ntk, std::vector<uint64_t> const& order )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_clone_node_v<Ntk>, "Ntk does not implement the clone_node method" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );
  static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  assert( order.size() == ntk.num_pos() );

  detail::permute_outputs_impl p( ntk, order );
  return p.run();
}

/*! \brief Reverse outputs permutation.
 *
 * This method reverts the outputs permutation performed with permute_outputs
 *
   \verbatim embed:rst

   .. note::

    This method restores the network outputs order assuming a previous application of permute_outputs.  It does
    *not* modify the input network.
   \endverbatim
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_node`
 * - `foreach_pi`
 * - `foreach_po`
 * - `clone_node`
 * - `is_pi`
 * - `is_constant`
 */
template<class Ntk>
Ntk permute_outputs_back( Ntk& ntk, std::vector<uint64_t> const& order )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_clone_node_v<Ntk>, "Ntk does not implement the clone_node method" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );
  static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  assert( order.size() == ntk.num_pos() );

  std::vector<uint64_t> rev_order = order;
  for( uint32_t i = 0; i < order.size(); ++i )
    rev_order[order[i]] = i;

  detail::permute_outputs_impl p( ntk, rev_order );
  return p.run();
}

} // namespace mockturtle
