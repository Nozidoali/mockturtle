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
#include <experiments.hpp>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/esop_balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/experimental/cost_generic_resub.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/utils/cost_functions.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_minmc2.hpp>
#include <mockturtle/algorithms/experimental/window_resub.hpp>


int main(int argc, char ** argv)
{
    
    using namespace mockturtle;
    using namespace mockturtle::experimental;

    std::string benchmark = argv[1];

    xag_network xag;

    if ( std::string( benchmark ).find( ".v" ) != std::string::npos )
    {
        lorina::read_verilog( benchmark, verilog_reader( xag ) );
    }
    else if ( std::string( benchmark ).find( ".aig" ) != std::string::npos )
    {
        lorina::read_aiger( benchmark, aiger_reader( xag ) );
    }

    // xag = convert_klut_to_graph<aig_network>( klut );
    assert( result == lorina::return_code::success );

    /* the cost functions for evalution */
    auto costfn = xag_multiplicative_complexity_cost_function<xag_network>();

    int c3 = cost_view( xag, costfn ).get_cost();

    cost_generic_resub_params ps;
    cost_generic_resub_stats st;
    ps.verbose = false;
    ps.rps.use_esop = true;
    ps.rps.max_solutions = 0; /* = 1: collect one, =0: collect all */

    resubstitution_params sps;

    sps.max_inserts = 20;
    sps.max_pis = 8;
    sps.max_divisors = std::numeric_limits<uint32_t>::max();

    window_resub_params wps;
    wps.wps.max_inserts = 3;

    future::xag_minmc_resynthesis resyn;
    cut_rewriting_params cps;
    cps.cut_enumeration_ps.cut_size = 4;

    xag_network xag_opt = xag.clone();
    while ( true ) {
        int cost_before = cost_view( xag, costfn ).get_cost();
        
        // window_xag_heuristic_resub( xag, wps );

        sim_resubstitution( xag, sps );
        xag = cleanup_dangling( xag );
        
        // xag = balancing( xag, { esop_rebalancing<xag_network>{} } );
        // xag = cleanup_dangling( xag );

        while ( true )
        {
            int cost_before_rewriting = cost_view( xag, costfn ).get_cost();
            xag = cut_rewriting<xag_network, decltype( resyn ), mc_cost<xag_network>>( xag, resyn, cps, nullptr );
            xag = cleanup_dangling( xag );

            cost_generic_resub( xag, costfn, ps, &st );
            xag = cleanup_dangling( xag );
            
            int cost_after_rewriting = cost_view( xag, costfn ).get_cost();
            if ( cost_after_rewriting >= cost_before_rewriting )
            {
                break;
            }
        }

        
        int cost_after = cost_view( xag, costfn ).get_cost();
        if ( cost_after >= cost_before ) {
            break;
        }
        xag_opt = xag.clone();
    }

    int _c3 = cost_view( xag_opt, costfn ).get_cost();

    std::string output_filename = argv[2];

    write_verilog( xag_opt, output_filename );

    /* the costs after optimization */

    fmt::print( "num_ands = {}\n", _c3 );

    return 0;
}