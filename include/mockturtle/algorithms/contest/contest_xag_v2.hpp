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

#include "evolutionary_algorithm.hpp"
#include "xag_moves.hpp"
#include "xag_initializations.hpp"

namespace mockturtle
{
namespace contest
{

class contest_method_xag_params
{
public:
  uint32_t search_length = 20;
  uint32_t num_search = 20;
};
class contest_method_xag
{
public:
  contest_method_xag(contest_method_xag_params const& ps = {})
    : ps(ps)
  {

  }

  xag_network run(klut_network const& _klut)
  {

    evolutionary_algorithm_params params;
    params.mutation_rate = 0.5;
    params.num_parents = 4;
    params.num_offsprings = 5;
    params.num_generations = 100;
    params.size_limit = 10000;

    auto xags = xag_initializations(_klut);

    auto xag = evolutionary_algorithm<xag_network>(xags.begin(), xags.end(), xag_compression, xag_decompression, params);

    return xag;
  }

  std::string name() const
  {
    return "v2";
  }
private:

  contest_method_xag_params ps;
};

} // namespace contest
} // namespace mockturtle
