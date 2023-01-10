#include <lorina/blif.hpp>
#include <mockturtle/mockturtle.hpp>

int main()
{
  mockturtle::names_view<mockturtle::sequential<mockturtle::klut_network>> klut;
  const auto result = lorina::read_blif( "./error.blif", mockturtle::blif_reader( klut ) );
  if ( result != lorina::return_code::success )
  {
    std::cout << "Read benchmark failed\n";
    return -1;
  }
  mockturtle::mig_npn_resynthesis resyn;
  mockturtle::names_view<mockturtle::sequential<mockturtle::mig_network>> named_dest;
  mockturtle::node_resynthesis( named_dest, klut, resyn );
  mockturtle::write_blif( named_dest, "output.blif" );
  return 0;
}