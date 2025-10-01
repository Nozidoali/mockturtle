# Custom 6-Input Heuristic NPN Resynthesis

## Overview

This directory contains a 6-input heuristic NPN resynthesis implementation for cut rewriting.

## Key Features

- **6-input cuts** (vs standard 4)
- **Heuristic NPN** canonization (swap-adjacent)
- **Incomplete database** (hash-based, populated at build time)
- **Runtime extensible** (add circuits via `add_to_database()`)

## Files

- `heuristic_npn_resynthesis.hpp` - Main implementation (~280 lines)
- `CMakeLists.txt` - Build configuration

## Database

**Current state**: Populated with 42 circuits covering common 2-3 input functions

- AND/XOR pairs: 30 circuits
- 3-input chains: 12 circuits  
- Total nodes: ~352
- Covered NPN classes: 5

## Usage

```cpp
#include <algorithms/node_resynthesis/heuristic_npn_resynthesis.hpp>

heuristic_npn_resynthesis<aig_network> resyn( params, nullptr, npn_params, 6 );
cut_rewriting_with_compatibility_graph( aig, resyn, ps, &st );
```

## Building

```bash
cd build
make cut_rewriting_npn_experiment
./experiments/cut_rewriting_npn_experiment
```

## Next Steps

To get actual improvements, populate with optimal circuits from synthesis tools (e.g., percy).

