#include <lorina/blif.hpp>
#include <mockturtle/mockturtle.hpp>

using namespace mockturtle;
int main()
{
  names_view<sequential<klut_network>> klut;
  const auto result = lorina::read_blif( "./error.blif", blif_reader( klut ) );
  if ( result != lorina::return_code::success )
  {
    std::cout << "Read benchmark failed\n";
    return -1;
  }
  mig_npn_resynthesis resyn;
  names_view<sequential<mig_network>> named_dest;
  node_resynthesis( named_dest, klut, resyn );
  named_dest = cleanup_dangling( named_dest );
  write_blif_params ps;
  ps.skip_feedthrough = false;
  write_blif( named_dest, "output1.blif", ps );
  return 0;
}