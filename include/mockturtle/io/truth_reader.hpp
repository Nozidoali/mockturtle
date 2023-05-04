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
  \file truth_reader.hpp
  \brief Lorina reader for TRUTH files

  \author Andrea Costamagna
*/

#pragma once

#include <kitty/kitty.hpp>
#include <lorina/truth.hpp>

#include <string>
#include <vector>

#include "../networks/klut.hpp"
#include "../traits.hpp"

namespace mockturtle
{

/*! \brief Lorina reader callback for TRUTH files.
 *
 * **Required network functions:**
 * - `create_pi`
 * - `create_po`
 * - `create_node` 
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      klut_network klut;
      lorina::read_truth( "file.truth", truth_reader( klut ) );

   \endverbatim
 */
template<typename Ntk>
class truth_reader : public lorina::truth_reader
{
public:
  explicit truth_reader( Ntk& ntk, bool reduce_support = false )
      : ntk_( ntk ),
      reduce( reduce_support )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi function" );
    static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po function" );
    static_assert( has_create_node_v<Ntk>, "Ntk does not implement the create_node function" );
    static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant function" );
  }

  ~truth_reader()
  {
    for ( auto const& o : outputs )
    {
      ntk_.create_po( o );
    }
  }

  virtual void on_input() const override
  {
    inputs.push_back( ntk_.create_pi() );
  }

  virtual void on_output( const std::string& tt_binary_string ) const override
  {
    auto num_inputs = log2( tt_binary_string.size() );
    kitty::dynamic_truth_table tt( num_inputs );
    kitty::create_from_binary_string( tt, tt_binary_string );
    auto copy = tt;
    const auto support = kitty::min_base_inplace( copy );
    auto support_size = static_cast<unsigned int>( support.size() );

    if( ( support_size < inputs.size() ) && reduce )
    { 
      std::vector<signal<Ntk>> small_inputs;
      for( auto i : support )
        small_inputs.push_back( inputs[i] );

      const auto small_func = kitty::shrink_to( copy, support_size );    
      outputs.emplace_back( ntk_.create_node( small_inputs, small_func ) );
    }
    else
      outputs.emplace_back( ntk_.create_node( inputs, tt ) );
  }

private:
  Ntk& ntk_;

  mutable std::vector<signal<Ntk>> inputs;
  mutable std::vector<signal<Ntk>> outputs;
  mutable bool reduce;
}; /* truth_reader */

} /* namespace mockturtle */