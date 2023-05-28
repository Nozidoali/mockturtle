template<class Ntk>
Ntk init_flow_hdp( klut_network& klut, int const& topology = 1 )
{

  size_t num_bits = pow( 2, klut.num_pis() );

  std::vector<kitty::partial_truth_table> examples;
  for( size_t i = 0; i < klut.num_pis(); ++i )
  {
    examples.push_back( kitty::partial_truth_table(num_bits) );
    create_nth_var( examples[i], i );
  }
  std::vector<kitty::partial_truth_table> targets;

  partial_simulator sim( examples );  
  unordered_node_map<kitty::partial_truth_table, klut_network> node_to_value( klut );
  simulate_nodes<kitty::partial_truth_table>( klut, node_to_value, sim );
    
  size_t i = 0;

  klut.foreach_po( [&]( auto const& node, auto index ) {
    targets.push_back( kitty::partial_truth_table(num_bits) );
    std::string tt_str = kitty::to_binary( node_to_value[ node ] );
    kitty::create_from_binary_string( targets[i], tt_str );
    i++;
  } );
  
  klut_network klut2 = hdc::project_in_hd( examples, targets, topology );

  auto ntk = convert_klut_to_graph<Ntk>( klut2 );

  ntk = cleanup_dangling( ntk );
  // ntk = abc_opto( ntk, benchmark, "&get; &deepsyn -I 3 -J 200; &put" );

  return ntk;
}
