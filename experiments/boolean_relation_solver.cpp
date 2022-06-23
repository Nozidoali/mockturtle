#include "experiments.hpp"
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/experimental/br_solver.hpp>
#include <mockturtle/networks/aig.hpp>

#include <iostream>
#include <string>

int main()
{
  using namespace mockturtle;
  using namespace mockturtle::experimental;

  std::vector<std::string> br_inputs = {
      // benchmark 1
      // "00001111",
      // "00110011",
      // "01010101"

      // benchmark 2
      "10101010",
      "11001100",
      "11110000"

  };
  std::vector<std::string> br_divisors = {};
  std::vector<std::string> br_outputs = {
      // benchmark 1
      // "11110000", // real input
      // "00011110",
      // "00011110",
      // "00010001"

      // benchmark 2
      "00010111",
      "11101000" };
  uint32_t tt_size = 8;

  using TT = kitty::partial_truth_table;

  std::vector<TT> targets;
  std::vector<TT> divs;
  for ( std::string br_input : br_inputs )
  {
    TT tt( tt_size );
    kitty::create_from_binary_string( tt, br_input );
    divs.emplace_back( tt );
  }
  for ( std::string br_divisor : br_divisors )
  {
    TT tt( tt_size );
    kitty::create_from_binary_string( tt, br_divisor );
    divs.emplace_back( tt );
  }
  for ( std::string br_target : br_outputs )
  {
    TT tt( tt_size );
    kitty::create_from_binary_string( tt, br_target );
    targets.emplace_back( tt );
  }

  br_solver engine;
  auto res = engine( divs, targets );
  if ( res )
  {
    std::cout << "solution : " << to_index_list_string( *res ) << std::endl;
  }
  else
  {
    std::cout << "no solution found" << std::endl;
  }
  return 0;
}