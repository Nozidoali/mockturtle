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
  \file klut_to_graph.hpp
  \brief Convert a k-LUT network into AIG, XAG, MIG or XMG.
  \author Andrea Costamagna
*/

#pragma once

#include "../../io/truth_reader.hpp"
#include "../../networks/aig.hpp"
#include "../../networks/xag.hpp"
#include "../klut_to_graph.hpp"
#include "../sim_resub.hpp"

/* for balancing */
#include "../balancing.hpp"
#include "../balancing/esop_balancing.hpp"
#include "../balancing/sop_balancing.hpp"

/* for LUT mapping */
#include "../../views/mapping_view.hpp"
#include "../collapse_mapped.hpp"
#include "../lut_mapping.hpp"

/* for cut rewriting*/
#include "../cut_rewriting.hpp"
#include "../node_resynthesis/xag_npn.hpp"

#include "cec.hpp"

namespace mockturtle
{
namespace contest
{

class contest_method_xag_params
{
public:
  uint32_t search_length = 10;
  uint32_t num_search = 1;
};
class contest_method_xag
{
public:
  contest_method_xag(contest_method_xag_params const& ps = {})
    : ps(ps)
  {

  }

  xag_network run(klut_network const& klut)
  {
    // get the initial network
    xag_network xag;
    convert_klut_to_graph( xag, klut );
    
    // optimization
    resubstitution_params resubstitution_parameters;
    resubstitution_parameters.max_inserts = 3;
    resubstitution_parameters.max_divisors = 1000;
    resubstitution_parameters.max_pis = 20;
    resubstitution_parameters.odc_levels = 3;
    resubstitution_parameters.conflict_limit = 1000000;
    resubstitution_parameters.max_clauses = 100000;

    xag_network best_xag = xag.clone();
    uint32_t prev_size = xag.num_gates();
    xag_npn_resynthesis<xag_network> resyn;
    cut_rewriting_params rewriting_parameters;
    rewriting_parameters.cut_enumeration_ps.cut_size = 4;

    /* every search starts from this initial xag */
    xag_network start_xag = xag.clone();


    for ( uint32_t search_idx = 0; search_idx < ps.num_search; search_idx++ )
    {
      if ( xag.num_gates() > start_xag.num_gates() )
        xag = start_xag.clone();
      for ( uint32_t op_idx = 0; op_idx < ps.search_length; op_idx ++ )
      {
        /* high-effort logic minimization */
        while ( 1 )
        {
          uint32_t size_before = xag.num_gates();
          sim_resubstitution( xag, resubstitution_parameters );
          xag = cleanup_dangling( xag );
          cut_rewriting( xag, resyn, rewriting_parameters );
          xag = cleanup_dangling( xag );
          if ( xag.num_gates() >= size_before )
            break;
        }
        if ( xag.num_gates() < prev_size )
          best_xag = xag.clone();
          prev_size = xag.num_gates();

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
    }
    return best_xag;
  }

  std::string name() const
  {
    return "convert_klut_to_graph + high-effort sim-resub + balancing";
  }
private:

  contest_method_xag_params ps;
};

} // namespace contest
} // namespace mockturtle
