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

#include "experiments.hpp"


#include <fmt/format.h>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/experimental/sequential_mapping.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/sequential.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/algorithms/collapse_mapped_sequential.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/views/names_view.hpp>
#include <mockturtle/algorithms/retiming_network.hpp>



int main()
{

  std::string benchmark_dir = "../experiments/elastic_circuit/cutloopback";
  std::vector<std::string> benchmark_suite = {
    "covariance_float",
    "gaussian",
    "gemver",
    "gsum",
    "gsumif",
    "insertion_sort",
    "kmp",
    "matching_2",
    "matrix",
    "mvt_float",
    "stencil_2d"
  };

  using namespace mockturtle;
  using namespace mockturtle::experimental;
  using namespace experiments;

  auto abc_read = [&] ( auto const & benchmark ) 
  {
    std::string command = fmt::format( "abc -q \"r {}/{}.blif; write_blif /tmp/{}.blif\";", benchmark_dir, benchmark, benchmark );
    std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
    std::array<char, 128> buffer;
    std::string result;
    while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
    {
      result += buffer.data();
    }
    fmt::print( result );
  // names_view<sequential<klut_network>> sequential_klut;
    sequential<klut_network> sequential_klut;
    if ( lorina::read_blif( fmt::format( "/tmp/{}.blif", benchmark ), blif_reader( sequential_klut ) ) != lorina::return_code::success )
    {
      fmt::print( "FATAL: read blif failed\n" );
    }
    return sequential_klut;
  };


  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t> exp( "sequential_mapping", "benchmark", "#LUTs", "#FF", "CP", "#LUTs'", "#FF'", "CP'" );

  for ( auto const& benchmark : benchmark_suite )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    
    auto sequential_klut = abc_read( benchmark );

    uint32_t n_luts = sequential_klut.num_gates();
    uint32_t n_ffs = sequential_klut.num_registers();
    uint32_t cp = 0;

    sequential_klut = cleanup_dangling( sequential_klut );

    mapping_view<decltype(sequential_klut), true> viewed{sequential_klut};
    
    sequential_mapping_params ps;
    ps.cut_enumeration_ps.cut_size = 6;
    sequential_mapping<decltype(viewed), true>( viewed, ps );
    
    sequential_klut = *collapse_mapped_sequential_network<decltype(sequential_klut)>( viewed );

    retiming_network_params rps;
    rps.clock_period = 1;
    retiming_network( sequential_klut, rps );
    sequential_klut = cleanup_dangling( sequential_klut );

    exp( benchmark, n_luts, n_ffs, cp,
    sequential_klut.num_gates(), sequential_klut.num_registers(), rps.clock_period );
  }

  exp.save();
  exp.table();

  return 0;
}