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
#include "../experimental/cost_generic_resub.hpp"

#include "cec.hpp"

namespace mockturtle
{
namespace contest
{

class contest_method_xag_params
{

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
    auto costfn = xag_size_cost_function<xag_network>();
    experimental::cost_generic_resub( xag, costfn );
    xag = cleanup_dangling( xag );
    return xag;
  }

  std::string name() const
  {
    return "convert_klut_to_graph";
  }
private:

  contest_method_xag_params ps;
};

} // namespace contest
} // namespace mockturtle
