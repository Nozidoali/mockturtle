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


/**
 * Author: Hanyu Wang
 */
#pragma once

#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/hash.hpp>
#include <fmt/format.h>
#include <unordered_map>
#include <unordered_set>
#include <functional>

namespace mockturtle
{

enum class npn_canon_method : uint32_t
{
  exact_npn = 0,
  flip_swap = 1,
  sifting = 2
};

struct heuristic_npn_params
{
  npn_canon_method method{ npn_canon_method::flip_swap };
};

template<class Ntk, class DatabaseNtk = xag_network>
class heuristic_npn_resynthesis
{
public:
  using truth_table_t = kitty::dynamic_truth_table;

  explicit heuristic_npn_resynthesis( xag_npn_resynthesis_params const& ps = {},
                                      xag_npn_resynthesis_stats* pst = nullptr,
                                      uint32_t max_inputs = 6,
                                      heuristic_npn_params const& npn_ps = {} )
    : _ps( ps ), _pst( pst ), _max_inputs( max_inputs ), _npn_ps( npn_ps ),
      _xag4( ps, &_xag4_st )
  {
    build_database();
    if ( _ps.verbose )
    {
      const char* method_name = _npn_ps.method == npn_canon_method::exact_npn ? "exact" :
                                _npn_ps.method == npn_canon_method::sifting ? "sifting" : "flip_swap";
      fmt::print( "[i] Heuristic NPN resynthesis: {} inputs, method: {}\n", _max_inputs, method_name );
      fmt::print( "[i]   4-input complete DB: {} entries, {} circuits\n",
                  _xag4_st.covered_classes, _xag4_st.db_size );
      fmt::print( "[i]   5-6 input heuristic: {} entries, {} circuits\n", 
                  _database.size(), _st.db_size );
      _st.report();
    }
  }

  ~heuristic_npn_resynthesis()
  {
    if ( _pst ) *_pst = _st;
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function,
                   LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    if ( function.num_vars() > _max_inputs ) return;
    if ( function.num_vars() <= 4 )
    {
      _xag4( ntk, function, begin, end, fn );
      return;
    }
    (void)fn;
  }

  size_t get_database_size() const { return _database.size(); }
  size_t get_total_circuits() const { return _st.db_size; }
  size_t get_covered_classes() const { return _st.covered_classes; }

  std::tuple<truth_table_t, uint32_t, std::vector<uint8_t>> npn_canonize( const truth_table_t& tt ) const
  {
    switch ( _npn_ps.method )
    {
      case npn_canon_method::exact_npn:
        return kitty::exact_npn_canonization( tt );
      case npn_canon_method::sifting:
        return kitty::sifting_npn_canonization( tt );
      case npn_canon_method::flip_swap:
      default:
        return kitty::flip_swap_npn_canonization( tt );
    }
  }

  void print_stats() const
  {
    size_t total_nodes = 0, total_gates = 0;
    for ( const auto& [hash, circuits] : _database )
      for ( const auto& c : circuits )
        if ( !c.empty() )
        {
          total_gates += ( c[0] >> 16 );
          total_nodes += ( 1 + (c[0] & 0xff) + (c[0] >> 16) );
        }

    fmt::print( "\n[i] Database: {} entries, {} circuits, {} nodes, {} gates\n",
                _database.size(), _st.db_size, total_nodes, total_gates );
    fmt::print( "[i] Status: POPULATED\n\n" );
  }

private:
  void build_database()
  {
    stopwatch t( _st.time_db );
    fmt::print( "[i] Building {}-input database...\n", _max_inputs );
    _database.clear();
    _st.db_size = 0;
    _st.covered_classes = 0;
    import_4input_database();
    populate_common_functions();
    print_stats();
  }

  void import_4input_database()
  {
    fmt::print( "[i] Importing 4-input complete database...\n" );
    size_t imported = 0, circuits_imported = 0;
    
    kitty::static_truth_table<4u> tt;
    do {
      kitty::dynamic_truth_table tt_dyn( 4u );
      kitty::create_from_words( tt_dyn, tt.cbegin(), tt.cend() );
      auto tt_ext = kitty::extend_to( tt_dyn, _max_inputs );
      const auto [repr, phase, perm] = npn_canonize( tt_ext );
      uint64_t hash = kitty::hash<truth_table_t>{}( repr );
      
      bool is_new = ( _database.find( hash ) == _database.end() );
      auto circuits = extract_4input_circuits( tt );
      
      for ( const auto& circuit : circuits )
      {
        _database[hash].push_back( circuit );
        _st.db_size++;
        circuits_imported++;
      }
      
      if ( is_new && !circuits.empty() )
      {
        _st.covered_classes++;
        imported++;
      }
      
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );
    
    fmt::print( "[i] Imported {} NPN classes, {} circuits\n", imported, circuits_imported );
  }

  std::vector<std::vector<uint32_t>> extract_4input_circuits( kitty::static_truth_table<4u> const& tt ) const
  {
    std::vector<std::vector<uint32_t>> circuits;
    Ntk temp_ntk;
    std::vector<signal<Ntk>> pis;
    for ( uint32_t i = 0; i < 4; ++i )
      pis.push_back( temp_ntk.create_pi() );
    
    kitty::dynamic_truth_table tt_dyn( 4u );
    kitty::create_from_words( tt_dyn, tt.cbegin(), tt.cend() );
    
    _xag4( temp_ntk, tt_dyn, pis.begin(), pis.end(), [&]( signal<Ntk> const& f ) {
      auto encoded = encode_from_network( temp_ntk, f, 4 );
      if ( !encoded.empty() )
        circuits.push_back( encoded );
      return false;
    });
    
    return circuits;
  }

  template<typename TempNtk>
  std::vector<uint32_t> encode_from_network( TempNtk const& ntk, signal<TempNtk> const& output, uint32_t num_pis ) const
  {
    std::vector<uint32_t> circuit;
    std::unordered_map<node<TempNtk>, uint32_t> node_to_idx;
    node_to_idx[0] = 0;
    uint32_t next_idx = 1;
    
    ntk.foreach_pi( [&]( auto n, uint32_t ) { node_to_idx[n] = next_idx++; });
    
    std::vector<std::pair<uint32_t, uint32_t>> gates;
    std::unordered_set<node<TempNtk>> visited;
    
    std::function<void(node<TempNtk>)> visit = [&]( node<TempNtk> n ) {
      if ( visited.count( n ) || node_to_idx.count( n ) ) return;
      ntk.foreach_fanin( n, [&]( auto const& fi ) { visit( ntk.get_node( fi ) ); });
      visited.insert( n );
      
      uint32_t lit0 = 0, lit1 = 0, fanin_count = 0;
      ntk.foreach_fanin( n, [&]( auto const& fi ) {
        uint32_t idx = node_to_idx[ntk.get_node( fi )];
        uint32_t lit = idx * 2 + ( ntk.is_complemented( fi ) ? 1 : 0 );
        if ( fanin_count == 0 ) lit0 = lit;
        else if ( fanin_count == 1 ) lit1 = lit;
        fanin_count++;
      });
      
      if ( fanin_count == 2 )
      {
        if ( ntk.is_xor( n ) && lit0 < lit1 ) std::swap( lit0, lit1 );
        else if ( !ntk.is_xor( n ) && lit0 > lit1 ) std::swap( lit0, lit1 );
        gates.push_back( { lit0, lit1 } );
        node_to_idx[n] = next_idx++;
      }
    };
    
    visit( ntk.get_node( output ) );
    
    if ( gates.empty() ) return circuit;
    
    circuit.push_back( num_pis | (1 << 8) | (gates.size() << 16) );
    for ( const auto& [lit0, lit1] : gates )
    {
      circuit.push_back( lit0 );
      circuit.push_back( lit1 );
    }
    
    uint32_t out_idx = node_to_idx[ntk.get_node( output )];
    circuit.push_back( out_idx * 2 + ( ntk.is_complemented( output ) ? 1 : 0 ) );
    
    return circuit;
  }

  void populate_common_functions()
  {
    size_t count_2 = 0, count_3 = 0, count_4 = 0, start = _st.db_size;

    for ( uint32_t i = 0; i < _max_inputs; ++i )
      for ( uint32_t j = i + 1; j < _max_inputs; ++j )
      {
        add_2input_variants( i, j );
        count_2 += 6;
      }

    if ( _max_inputs >= 3 )
      for ( uint32_t i = 0; i < _max_inputs; ++i )
        for ( uint32_t j = i + 1; j < _max_inputs; ++j )
          for ( uint32_t k = j + 1; k < _max_inputs; ++k )
          {
            add_3input_functions( i, j, k );
            count_3 += 10;
          }

    if ( _max_inputs >= 4 )
      for ( uint32_t i = 0; i < std::min(_max_inputs, 5u); ++i )
        for ( uint32_t j = i + 1; j < std::min(_max_inputs, 5u); ++j )
          for ( uint32_t k = j + 1; k < std::min(_max_inputs, 5u); ++k )
            for ( uint32_t l = k + 1; l < std::min(_max_inputs, 5u); ++l )
            {
              add_4input_patterns( i, j, k, l );
              count_4 += 5;
            }

    add_complex_functions();

    fmt::print( "[i] Populated database:\n" );
    fmt::print( "[i]   2-input: {}, 3-input: {}, 4-input: {}, complex: {}\n", 
                count_2, count_3, count_4, _st.db_size - start - count_2 - count_3 - count_4 );
    fmt::print( "[i]   Total: {} circuits\n", _st.db_size - start );
  }

  void add_2input_variants( uint32_t i, uint32_t j )
  {
    truth_table_t tti( _max_inputs ), ttj( _max_inputs );
    kitty::create_nth_var( tti, i );
    kitty::create_nth_var( ttj, j );

    add_internal( tti & ttj, encode_2gate( i, j, false, false ) );
    add_internal( tti ^ ttj, encode_2gate( i, j, false, true ) );
    add_internal( tti | ttj, encode_or( i, j ) );
    add_internal( ~(tti & ttj), encode_nand( i, j ) );
    add_internal( ~(tti | ttj), encode_nor( i, j ) );
    add_internal( ~(tti ^ ttj), encode_xnor( i, j ) );
  }

  void add_3input_functions( uint32_t i, uint32_t j, uint32_t k )
  {
    truth_table_t tti( _max_inputs ), ttj( _max_inputs ), ttk( _max_inputs );
    kitty::create_nth_var( tti, i );
    kitty::create_nth_var( ttj, j );
    kitty::create_nth_var( ttk, k );

    add_internal( tti & ttj & ttk, encode_and_chain( i, j, k ) );
    add_internal( tti ^ ttj ^ ttk, encode_xor_chain( i, j, k ) );
    add_internal( tti | ttj | ttk, encode_or_chain( i, j, k ) );
    add_internal( (tti & ttj) | ttk, encode_and_or( i, j, k ) );
    add_internal( (tti | ttj) & ttk, encode_or_and( i, j, k ) );
    add_internal( (tti ^ ttj) & ttk, encode_xor_and( i, j, k ) );
    add_internal( (tti & ttj) ^ ttk, encode_and_xor( i, j, k ) );
    add_internal( (tti & ttj) | (tti & ttk) | (ttj & ttk), encode_maj( i, j, k ) );
    add_internal( tti ^ ttj ^ (ttj & ttk), encode_mux( i, j, k ) );
    add_internal( (tti & ttj) | (~tti & ttk), encode_mux2( i, j, k ) );
  }

  void add_4input_patterns( uint32_t i, uint32_t j, uint32_t k, uint32_t l )
  {
    truth_table_t tti( _max_inputs ), ttj( _max_inputs ), ttk( _max_inputs ), ttl( _max_inputs );
    kitty::create_nth_var( tti, i );
    kitty::create_nth_var( ttj, j );
    kitty::create_nth_var( ttk, k );
    kitty::create_nth_var( ttl, l );

    add_internal( tti & ttj & ttk & ttl, encode_4and( i, j, k, l ) );
    add_internal( tti ^ ttj ^ ttk ^ ttl, encode_4xor( i, j, k, l ) );
    add_internal( (tti & ttj) | (ttk & ttl), encode_and_or_tree( i, j, k, l ) );
    add_internal( (tti ^ ttj) ^ (ttk ^ ttl), encode_xor_tree( i, j, k, l ) );
    add_internal( (tti & ttj) | (ttk & ttl), encode_or_of_ands( i, j, k, l ) );
  }

  void add_complex_functions()
  {
    if ( _max_inputs >= 3 )
      for ( uint32_t i = 0; i + 2 < _max_inputs && i < 3; ++i )
      {
        truth_table_t a( _max_inputs ), b( _max_inputs ), c( _max_inputs );
        kitty::create_nth_var( a, i );
        kitty::create_nth_var( b, i + 1 );
        kitty::create_nth_var( c, i + 2 );
        add_internal( a ^ b ^ c, encode_xor_chain( i, i+1, i+2 ) );
        add_internal( (a & b) | (c & (a ^ b)), encode_carry( i, i+1, i+2 ) );
      }
  }

  std::vector<uint32_t> encode_2gate( uint32_t i, uint32_t j, bool complement, bool is_xor ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (1 << 16) );
    uint32_t lit0 = (i+1)*2 + (complement ? 1 : 0);
    uint32_t lit1 = (j+1)*2;
    if ( is_xor ) { circuit.push_back( lit1 ); circuit.push_back( lit0 ); }
    else { circuit.push_back( lit0 ); circuit.push_back( lit1 ); }
    circuit.push_back( (_max_inputs+1)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_or( uint32_t i, uint32_t j ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (i+1)*2 + 1 ); circuit.push_back( (j+1)*2 + 1 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+3)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_nand( uint32_t i, uint32_t j ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (1 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_nor( uint32_t i, uint32_t j ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (i+1)*2 + 1 ); circuit.push_back( (j+1)*2 + 1 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_xnor( uint32_t i, uint32_t j ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (1 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_and_chain( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_xor_chain( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (_max_inputs+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_or_chain( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 + 1 ); circuit.push_back( (j+1)*2 + 1 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+3)*2 + 1 ); circuit.push_back( (k+1)*2 + 1 );
    circuit.push_back( (_max_inputs+4)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_and_or( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 + 1 ); circuit.push_back( (k+1)*2 + 1 );
    circuit.push_back( (_max_inputs+2)*2 ); circuit.push_back( (_max_inputs+3)*2 );
    circuit.push_back( (_max_inputs+4)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_or_and( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 + 1 ); circuit.push_back( (j+1)*2 + 1 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+3)*2 + 1 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+4)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_xor_and( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_and_xor( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (2 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (_max_inputs+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_maj( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 + 1 );
    circuit.push_back( (_max_inputs+3)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_mux( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (_max_inputs+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_mux2( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (i+1)*2 + 1 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 + 1 );
    circuit.push_back( (_max_inputs+3)*2 + 1 );
    return circuit;
  }

  std::vector<uint32_t> encode_4and( uint32_t i, uint32_t j, uint32_t k, uint32_t l ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (l+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_4xor( uint32_t i, uint32_t j, uint32_t k, uint32_t l ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (l+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 ); circuit.push_back( (_max_inputs+1)*2 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_and_or_tree( uint32_t i, uint32_t j, uint32_t k, uint32_t l ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (l+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 + 1 ); circuit.push_back( (_max_inputs+2)*2 + 1 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_xor_tree( uint32_t i, uint32_t j, uint32_t k, uint32_t l ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (l+1)*2 ); circuit.push_back( (k+1)*2 );
    circuit.push_back( (_max_inputs+2)*2 ); circuit.push_back( (_max_inputs+1)*2 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_or_of_ands( uint32_t i, uint32_t j, uint32_t k, uint32_t l ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (l+1)*2 );
    circuit.push_back( (_max_inputs+1)*2 + 1 ); circuit.push_back( (_max_inputs+2)*2 + 1 );
    circuit.push_back( (_max_inputs+3)*2 );
    return circuit;
  }

  std::vector<uint32_t> encode_carry( uint32_t i, uint32_t j, uint32_t k ) const
  {
    std::vector<uint32_t> circuit;
    circuit.push_back( _max_inputs | (1 << 8) | (3 << 16) );
    circuit.push_back( (i+1)*2 ); circuit.push_back( (j+1)*2 );
    circuit.push_back( (j+1)*2 ); circuit.push_back( (i+1)*2 );
    circuit.push_back( (k+1)*2 ); circuit.push_back( (_max_inputs+2)*2 );
    circuit.push_back( (_max_inputs+1)*2 );
    return circuit;
  }

  void add_internal( const truth_table_t& tt, const std::vector<uint32_t>& circuit )
  {
    const auto [repr, phase, perm] = npn_canonize( tt );
    uint64_t hash = kitty::hash<truth_table_t>{}( repr );
    bool is_new = ( _database.find( hash ) == _database.end() );
    _database[hash].push_back( circuit );
    _st.db_size++;
    if ( is_new ) _st.covered_classes++;
  }

  xag_npn_resynthesis_params _ps;
  xag_npn_resynthesis_stats* _pst{ nullptr };
  uint32_t _max_inputs;
  heuristic_npn_params _npn_ps;
  mutable xag_npn_resynthesis_stats _st;
  std::unordered_map<uint64_t, std::vector<std::vector<uint32_t>>> _database;
  xag_npn_resynthesis_stats _xag4_st;
  xag_npn_resynthesis<Ntk, DatabaseNtk, xag_npn_db_kind::xag_complete> _xag4;
};

} /* namespace mockturtle */

