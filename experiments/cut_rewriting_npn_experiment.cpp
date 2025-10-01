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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>

#include <algorithms/node_resynthesis/heuristic_npn_resynthesis.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, bool> exp( 
    "cut_rewriting_npn_experiment", "benchmark", "size_before", "size_after", "runtime", "equivalent" );

  xag_npn_resynthesis_params resyn_params;
  resyn_params.verbose = true;

  heuristic_npn_params npn_params;
  npn_params.method = npn_canon_method::flip_swap;

  heuristic_npn_resynthesis<aig_network> resyn( resyn_params, nullptr, 6, npn_params );

  cut_rewriting_params ps;
  ps.cut_enumeration_ps.cut_size = 6;
  ps.cut_enumeration_ps.cut_limit = 10;
  ps.progress = true;
  ps.allow_zero_gain = false;
  ps.use_dont_cares = false;

  fmt::print( "[i] ========================================\n" );
  fmt::print( "[i] 6-INPUT HEURISTIC NPN\n" );
  fmt::print( "[i] ========================================\n" );
  fmt::print( "[i] Cut size: {}, Cut limit: {}\n", 
              ps.cut_enumeration_ps.cut_size, ps.cut_enumeration_ps.cut_limit );
  fmt::print( "[i] Database: {} entries, {} circuits\n",
              resyn.get_database_size(), resyn.get_total_circuits() );
  fmt::print( "[i] ========================================\n\n" );

  auto benchmarks = epfl_benchmarks();
  fmt::print( "[i] Found {} EPFL benchmarks\n", benchmarks.size() );

  std::vector<std::string> filtered_benchmarks;
  for ( auto const& benchmark : benchmarks )
  {
    aig_network temp_aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( temp_aig ) ) != lorina::return_code::success )
      continue;
    
    if ( temp_aig.num_gates() < 200000 )
    {
      filtered_benchmarks.push_back( benchmark );
      fmt::print( "[i]   {} ({} nodes) - included\n", benchmark, temp_aig.num_gates() );
    }
    else
      fmt::print( "[i]   {} ({} nodes) - skipped\n", benchmark, temp_aig.num_gates() );
  }

  fmt::print( "[i] Processing {} benchmarks (< 2K nodes)\n\n", filtered_benchmarks.size() );

  uint32_t processed = 0, successful = 0;

  for ( auto const& benchmark : filtered_benchmarks )
  {
    fmt::print( "[i] Processing: {}\n", benchmark );
    
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
      continue;

    uint32_t size_before = aig.num_gates();
    auto start = std::chrono::high_resolution_clock::now();
    cut_rewriting_with_compatibility_graph( aig, resyn, ps );
    auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - start).count() / 1000.0;

    aig = cleanup_dangling( aig );
    uint32_t size_after = aig.num_gates();
    bool equivalent = abc_cec( aig, benchmark );

    fmt::print( "    {} → {} gates ({:.1f}%) {} {:.2f}s\n",
                size_before, size_after,
                100.0 * (size_before - size_after) / size_before,
                equivalent ? "PASS" : "FAIL", runtime );

    exp( benchmark, size_before, size_after, runtime, equivalent );
    processed++;
    if ( equivalent ) successful++;
  }

  exp.save();
  fmt::print( "\n[i] Completed: {}/{} PASS\n", successful, processed );
  exp.table();

  return 0;
}
