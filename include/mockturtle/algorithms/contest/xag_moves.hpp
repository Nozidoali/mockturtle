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

/* for balancing */
#include "../balancing.hpp"
#include "../balancing/esop_balancing.hpp"
#include "../balancing/sop_balancing.hpp"

/* for cut rewriting*/
#include "../cut_rewriting.hpp"
#include "../node_resynthesis/xag_npn.hpp"

/* for optimization */
#include "../sim_resub.hpp"

/* for LUT mapping */
#include "../../views/mapping_view.hpp"
#include "../collapse_mapped.hpp"
#include "../lut_mapping.hpp"
#include "../klut_to_graph.hpp"


#include <algorithm>
#include <random>

namespace mockturtle
{

void xag_compression(xag_network& xag)
{

  xag_network xag_opt = xag.clone();
    
  // optimization
  resubstitution_params resubstitution_parameters;
  resubstitution_parameters.max_inserts = 3;
  resubstitution_parameters.max_divisors = 1000;
  resubstitution_parameters.max_pis = 20;
  resubstitution_parameters.odc_levels = 3;
  resubstitution_parameters.conflict_limit = 1000000;
  resubstitution_parameters.max_clauses = 100000;

  // rewriting
  xag_npn_resynthesis<xag_network> resyn;
  cut_rewriting_params rewriting_parameters;
  rewriting_parameters.cut_enumeration_ps.cut_size = 4;

  while ( 1 )
  {
    uint32_t size_before = xag_opt.num_gates();
    sim_resubstitution( xag_opt, resubstitution_parameters );
    xag_opt = cleanup_dangling( xag_opt );
    cut_rewriting( xag_opt, resyn, rewriting_parameters );
    xag_opt = cleanup_dangling( xag_opt );
    if ( xag_opt.num_gates() >= size_before )
      break;
  }

  xag = xag_opt.clone();
}

void xag_decompression(xag_network& xag)
{
  /* restructure, perturbation */
  uint32_t cut_size = random() % 4 + 3;
  lut_mapping_params ps;
  mapping_view<xag_network, true> mapped_xag{ xag };
  lut_mapping<decltype( mapped_xag ), true>( mapped_xag, ps );
  ps.cut_enumeration_ps.cut_size = cut_size;
  const auto klut = *collapse_mapped_network<klut_network>( mapped_xag );
  xag = convert_klut_to_graph<xag_network>( klut );

  xag = balancing( xag, { sop_rebalancing<xag_network>{} } );
  xag = cleanup_dangling( xag );
}

} // namespace mockturtle