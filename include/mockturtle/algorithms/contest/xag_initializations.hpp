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
  \brief klut to ntk method for IWLS contest

  \author Hanyu WANG
*/

#pragma once

/* network */
#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include "../../networks/klut.hpp"


#include "../klut_to_graph.hpp"

/* for output permutation */
#include "../permute_outputs.hpp"

namespace mockturtle
{

std::vector<xag_network> xag_initializations(klut_network const& _klut)
{

  std::vector<xag_network> xags;

  // clone the klut network
  klut_network klut = _klut;

  /* random number generator */
  std::mt19937 g( 888 );
  std::vector<uint64_t> outputs_order( klut.num_pos() );
  std::iota( std::begin( outputs_order ), std::end( outputs_order ), 0 );

  uint32_t n_inits = 6u;
  for ( uint32_t i = 0; i < n_inits; ++i )
  {
    std::shuffle( outputs_order.begin(), outputs_order.end(), g );

    klut_network output_permuted_klut = permute_outputs( klut, outputs_order );
    xag_network xag = convert_klut_to_graph<xag_network>( output_permuted_klut );
    xag = permute_outputs_back( xag, outputs_order );

    xag = cleanup_dangling( xag );

    xags.push_back( xag );
  }

  return xags;
}

} // namespace mockturtle