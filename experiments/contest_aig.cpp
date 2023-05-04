#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/io/truth_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <lorina/truth.hpp>
#include "experiments.hpp"
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <mockturtle/algorithms/contest/cec.hpp>
#include <mockturtle/algorithms/contest/contest_aig.hpp>

#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
#include <random>

#include <iostream>

using namespace mockturtle;

using experiment_t = experiments::experiment<std::string, uint32_t, uint32_t, std::string>;

experiment_t exp_res( "contest_aig", "benchmark", "#gates", "depth", "method" );

struct contest_parameters
{
  bool verbose{false};
};

#pragma region mutex
std::atomic<uint32_t> exp_id{0};
std::mutex exp_mutex;
#pragma endregion

void thread_run( contest_parameters const& ps, std::string const& run_only_one )
{
  const std::string benchmark_path = "../experiments/contest_benchmarks/";
  const std::string output_path = "../experiments/contest_results/aigs/";

  uint32_t id = exp_id++;

  while ( id < 100 )
  {
    /* Step 1: Read benchmarks */
    std::string benchmark = fmt::format( "ex{:02}", id );
    if ( run_only_one != "" && benchmark != run_only_one )
    {
      id = exp_id++;
      continue;
    }
    auto current_best = *exp_res.get_entry<uint32_t>( benchmark, "#gates", "best" );
    std::cout << "[i] processing " << benchmark << " curr best = " << current_best << "\n";
    klut_network klut;
    auto res = lorina::read_truth( benchmark_path + benchmark + ".truth", truth_reader( klut ) );
    if ( res != lorina::return_code::success )
    {
      std::cout << "read failed\n";
      id = exp_id++;
      continue;
    }
    
    aig_network aig;
    auto start = std::chrono::high_resolution_clock::now();

    /* Step 2: define method here*/
    auto method = contest::contest_method_aig();
    aig = method.run( klut );

    if ( !abc_cec_truth( aig, klut, benchmark ) )
    {
      std::cout << "[w] cec = false!\n";
      id = exp_id++;
      continue;
    }

    /* Step 3: Evaluation */
    auto daig = depth_view{aig};
    if ( aig.num_gates() < current_best )
    {
      fmt::print( "[i] obtained better result on {}: {} < {}\n", benchmark, aig.num_gates(), current_best );
      exp_mutex.lock();
      exp_res( benchmark, aig.num_gates(), daig.depth(), method.name() );
      exp_mutex.unlock();
      aig = cleanup_dangling( aig );
      write_aiger( aig, output_path + benchmark + ".aig" );
    }
    else
    {
      fmt::print( "[i] obtained worse result on {}: {} >= {}\n", benchmark, aig.num_gates(), current_best );
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    if ( ps.verbose )
    {
      std::cout << ".b " << benchmark << "\n";
      std::cout << ".g " << aig.num_gates() << "\n";
      std::cout << ".d " << daig.depth() << "\n";
      std::cout << ".t " << duration.count() << std::endl;
    }

    id = exp_id++;
  }
}

int main( int argc, char* argv[] )
{
  using namespace experiments;

  contest_parameters ps_contest;
  ps_contest.verbose = false;

  std::string run_only_one = "";

  if ( argc == 2 )
    run_only_one = std::string( argv[1] );

  const auto processor_count = run_only_one != "" ? 1 : std::thread::hardware_concurrency();

  /* starting benchmark id */
  exp_id.store( 0 );
  std::vector<std::thread> threads;

  /* generate threads */
  fmt::print( "[i] Running on {} threads\n", processor_count );

  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads.emplace_back( thread_run, ps_contest, run_only_one );
  }

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads[i].join();
  }
  exp_res.update( "best" );
  return 0;
}
