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
#include <mockturtle/algorithms/contest/contest_xag.hpp>

#include <thread>
#include <mutex>
#include <algorithm>
#include <set>
#include <random>

#include <iostream>

using namespace mockturtle;

using experiment_t = experiments::experiment<std::string, uint32_t, uint32_t, std::string>;

experiment_t exp_res( "contest_xag", "benchmark", "#gates", "depth", "method" );

struct contest_parameters
{
  bool verbose{false};
};


/* 100 */
std::unordered_set<uint32_t> id_skipped = {0,1,2,5,8,9,11,15,16,18,20,22,24,25,26,27,28,29,31,33,34,36,37,38,40,42,43,45,47,48,50,51,54,56,58,59,60,61,64,66,67,69,70,72,73,74,75,76,77,78,80,82,85,86,87,88,91,92,94,98};

/* 1000 */

#pragma region mutex
std::atomic<uint32_t> exp_id{0};
std::mutex exp_mutex;
#pragma endregion

void thread_run( contest_parameters const& ps, std::string const& run_only_one )
{
  const std::string benchmark_path = "../experiments/contest_benchmarks/";
  const std::string output_path = "../experiments/contest_results/xags/";

  exp_mutex.lock();
  uint32_t id = exp_id++;
  exp_mutex.unlock();

  while ( id < 100 )
  {
    /* Step 1: Read benchmarks */
    std::string benchmark = fmt::format( "ex{:02}", id );
    if ( run_only_one != "" && benchmark != run_only_one )
    {
      exp_mutex.lock();
      id = exp_id++;
      exp_mutex.unlock();
      continue;
    }

    if ( id_skipped.find( id ) != id_skipped.end() )
    {
      exp_mutex.lock();
      id = exp_id++;
      exp_mutex.unlock();
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
    
    xag_network xag;
    auto start = std::chrono::high_resolution_clock::now();

    /* Step 2: define method here*/
    auto method = contest::contest_method_xag();
    xag = method.run( klut );

    if ( !abc_cec_truth( xag, klut, benchmark ) )
    {
      std::cout << "[w] cec = false!\n";
      id = exp_id++;
      continue;
    }

    /* Step 3: Evaluation */
    auto dxag = depth_view{xag};
    if ( xag.num_gates() < current_best )
    {
      fmt::print( "[i] obtained better result on {}: {} < {}\n", benchmark, xag.num_gates(), current_best );
      exp_mutex.lock();
      exp_res( benchmark, xag.num_gates(), dxag.depth(), method.name() );
      exp_res.update( "best" );
      exp_mutex.unlock();
      aig_network aig = cleanup_dangling<xag_network, aig_network>( xag );
      write_aiger( aig, output_path + benchmark + ".aig" );
    }
    else
    {
      fmt::print( "[i] obtained worse result on {}: {} >= {}\n", benchmark, xag.num_gates(), current_best );
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    if ( ps.verbose )
    {
      std::cout << ".b " << benchmark << "\n";
      std::cout << ".g " << xag.num_gates() << "\n";
      std::cout << ".d " << dxag.depth() << "\n";
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
