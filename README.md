# Algebra

A C++ library for working with finite combinatorial and algebraic structures **up to isomorphism**: graphs, permutation groups, and related objects. The core focus is fast canonical labeling, isomorphism testing, and enumeration of graph families modulo isomorphism.

The example programs apply this machinery to problems in **Ramsey theory** and **Turán-type extremal graph theory**.

## Features

### Library

| Component | Description |
|-----------|-------------|
| `Perm` | Permutations: composition, conjugation, cycle structure |
| `Group` | Permutation groups via Schreier–Sims; includes `S(n)`, `A(n)`, `Z(n)`, `D(n)`, and the Mathieu groups |
| `Structure` | Abstract base for combinatorial structures with canonical labeling |
| `Certificate` | Canonical form used as an isomorphism invariant |
| `Graph` | Simple graphs with standard constructors (`K`, `C`, `P`, `Q`) |
| `StructSet` / `GraphSet` | Thread-safe sets of pairwise non-isomorphic structures |

A structure is **certified** by computing a canonical `Certificate`. Two structures are isomorphic if and only if their certificates match. Automorphism groups are computed as a side effect of certification.

### Example programs

| Target | Source | Purpose |
|--------|--------|---------|
| `gen_trees` | `examples/Trees.cpp` | Enumerate trees up to isomorphism |
| `one_vertex` | `examples/RamseyExtension.cpp` | Build Ramsey graphs `R(G, k)` by one-vertex extension |
| `gluing` | `examples/RamseyGluing.cpp` | Build Ramsey graphs `R(G, k, ≥n)` by a gluing algorithm |
| `turan` | `examples/Turan.cpp` | Compute Turán extremal graphs `EX(n, G)` |

Computed catalogs are written under `data/` (gitignored). Summary tables for Ramsey counts live in `results/RAMSEY/`.

## Requirements

- C++17 compiler (GCC or Clang)
- CMake ≥ 3.5
- [Catch2](https://github.com/catchorg/Catch2) v2 (header-only; vendored under `catch/`)
- POSIX threads (`pthread`)

## Building

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Binaries are placed in the build directory. For a debug build with AddressSanitizer:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

## Tests

```bash
./build/test_permutation
./build/test_group
./build/test_graph
```

The tests cover permutations, group order and containment, graph constructors, isomorphism, and automorphism group orders.

## Example usage

### Enumerate trees

```bash
./build/gen_trees
```

Prints the number of non-isomorphic trees on `n` vertices for `n = 1, …, 20`.

### Ramsey graphs (one-vertex extension)

```bash
./build/one_vertex <G> <k> [threads]
```

- **G** — forbidden graph: `K3`, `K4`, `C3`, `C4`, `C5`, `3`, or `4`
- **k** — clique parameter
- **threads** — optional thread count (default: 8)

Output goes to `data/RAMSEY/R(G,k)/`.

### Ramsey graphs (gluing)

```bash
./build/gluing <G> <k> <n> [threads]
```

- **G** — forbidden graph: `K3`, `C3`, `C4`, or `3`
- **k**, **n** — compute `R(G, k, ≥n)`
- **threads** — optional thread count (default: 8)

### Turán extremal graphs

```bash
./build/turan <G> [threads]
```

- **G** — forbidden bipartite graph: `C4`, `C6`, `C8`, `K23`, `K24`, `K25`, `K26`, `K33`, `K34`, `K35`, `K36`, `K44`, or `Q3`
- **threads** — optional thread count (default: 8)

Output goes to `data/TURAN/<G>/`.

## Project layout

```
inc/           Public headers
src/           Library implementation
tests/         Catch2 unit tests
examples/      Enumeration tools (Ramsey, Turán, trees)
catch/         Vendored Catch2
cmake/         CMake modules (FindCatch.cmake)
results/       Precomputed summary tables
data/          Generated graph catalogs (not in git)
```

## Library sketch

```cpp
#include "Graph.h"

Graph G = C(5);
G.certify();                    // compute canonical certificate
Group aut = G.aut();            // automorphism group

Graph H = K(2, 3);
H.certify();

bool same = isomorphic(G, H);   // false for C(5) vs K(2,3)

GraphSet catalog;
catalog.insert(G.certify());    // store up to isomorphism
catalog.write("graphs.gr");     // binary certificate dump
```

Graphs are stored on disk as sequences of canonical certificates. Use `readGraph` to load them back.
