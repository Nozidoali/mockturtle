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
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/utils/cost_functions.hpp>
#include <mockturtle/utils/stopwatch.hpp>

using namespace mockturtle;


template<typename Ntk>
Ntk abc_opto( Ntk const& ntk, std::string str_code, std::string abc_script = "compress2rs" )
{
  write_aiger( ntk, "/tmp/" + str_code + ".aig" );
  std::string command = "abc -q \"r /tmp/" + str_code + ".aig; " + abc_script + "; write_aiger /tmp/" + str_code + ".aig\"";

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  Ntk res;
  std::string string_path = ( "/tmp/" + str_code + ".aig" );
  if ( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

int main()
{

  using namespace mockturtle::experimental;
  using namespace experiments;

  /* run the actual experiments */
  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool> exp( "cost_generic_resub", "benchmark", "c1", "c2", "c3", "_c1", "_c2", "_c3", "runtime", "impr %", "cec" );

  std::vector<std::string> bugs{ "bug" };
  for ( auto const& benchmark : iwls_benchmarks() )
  // for ( auto const& benchmark : bugs )
  {
    float run_time = 0;

    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    auto const result = lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
    assert( result == lorina::return_code::success );
    (void)result;

    if ( aig.num_gates() >= 120000 ) {
    // if ( benchmark != "vga_lcd" ) {
      fmt::print( "[w] skipping {} due to size\n", benchmark );
      continue;
    }

    /* the cost functions for evalution */
    auto costfn_1 = xag_size_cost_function<aig_network>();
    auto costfn_2 = xag_depth_cost_function<aig_network>();
    auto costfn_3 = aig_fflc_cost_function<aig_network>();
    auto costfn_4 = aig_fflc_cost_function_eval<aig_network>();

    /* the cost function for optimization */
    auto costfn = costfn_1;


    cost_generic_resub_params ps;
    cost_generic_resub_stats st;
    ps.verbose = false;
    ps.rps.use_esop = false;
    ps.rps.max_solutions = 0; /* = 1: collect one, =0: collect all */

    stopwatch<>::duration time_tot{ 0 };

    aig = abc_opto( aig, benchmark, "compress2rs" );
    aig = abc_opto( aig, benchmark, "compress2rs" );

    /* the costs before optimization */
    int c1 = cost_view( aig, costfn_1 ).get_cost();
    int c2 = cost_view( aig, costfn_2 ).get_cost();
    int c3 = cost_view( aig, costfn_4 ).get_cost() + aig.num_pos() - 1;

    call_with_stopwatch( time_tot, [&]() {
      for ( auto i = 0; i < 2; ++i ) {
        aig = abc_opto( aig, benchmark, "b -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 6" );
        ps.wps.max_pis = 6;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "rw -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 6 -N 2 -l" );
        ps.wps.max_pis = 6;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "rf -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 8 -l" );
        ps.wps.max_pis = 8;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "b -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 8 -N 2 -l" );
        ps.wps.max_pis = 8;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "rw -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 10 -l" );
        ps.wps.max_pis = 10;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );
        
        aig = abc_opto( aig, benchmark, "rwz -l" );
        // aig = abc_opto( aig, benchmark, "rs -K 10 -N 2 -l" );
        ps.wps.max_pis = 10;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "b -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 12 -l" );
        ps.wps.max_pis = 12;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );

        aig = abc_opto( aig, benchmark, "rfz -l" );

        // aig = abc_opto( aig, benchmark, "rs -K 12 -N 2 -l" );
        ps.wps.max_pis = 12;
        cost_generic_resub( aig, costfn, ps, &st );
        aig = cleanup_dangling( aig );
        
        aig = abc_opto( aig, benchmark, "rwz -l" );
        aig = abc_opto( aig, benchmark, "b -l" );
        aig = cleanup_dangling( aig );
      }
    } );

    run_time = to_seconds( time_tot );

    /* the costs after optimization */
    int _c1 = cost_view( aig, costfn_1 ).get_cost();
    int _c2 = cost_view( aig, costfn_2 ).get_cost();
    int _c3 = cost_view( aig, costfn_4 ).get_cost() + aig.num_pos() - 1;

    fmt::print( "[i] c1: {} -> {}, c2: {} -> {}, c3: {} -> {}\n", c1, _c1, c2, _c2, c3, _c3 );

    float impr = (float)( c3 - _c3 ) / (float)c3 * 100.0f;

    auto cec = true;
    // cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    exp( benchmark, c1, c2, c3, _c1, _c2, _c3, run_time, impr, cec );
  }
  exp.save();
  exp.table();
  return 0;
}