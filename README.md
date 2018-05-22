# Hash-based colored de Bruijn graph with wavelet trie and Bloom filter color compression

**Prerequisites**
- cmake 3.6.1
- C++14
- HTSlib
- GNU GMP (wavelet trie only)
- boost (wavelet trie only)

### Install
1. `git clone --recursive https://github.com/ratschlab/graph_annotation`
2. go to the **build** directory `mkdir -p annograph/build && cd annograph/build`
3. compile by `cmake .. && make && ./unit_tests`

Build types: `cmake .. <arguments>` where arguments are:

- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (Debug by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)

#### Install wavelet trie compressor
See [wavelet_trie](./wavelet_trie).

### Typical workflow
1. Generate graph and uncompressed annotations (`.precise.dbg` and `.rawrows.dbg` files)  
`./annograph build -o <OUTPREFIX> <FLAGS> <INPUTS>`
2. Compress annotation with Bloom filters  
`./annograph build -i <OUTPREFIXX> -o <BLOOMOUTPREFIX> <FLAGS> <INPUTS>`
3. Compress annotation with wavelet tries (see [wavelet_trie](./wavelet_trie))

## bioRxiv preprint
This code was used to produce the results in the bioRxiv preprint, _[Dynamic compression schemes for graph coloring](https://www.biorxiv.org/content/early/2018/03/17/239806)_ by Harun Mustafa, Ingo Schilken, Mikhail Karasikov, Carsten Eickhoff, Gunnar Ratsch, and Andre Kahles. 

The input data used for those experiments is located [here](https://public.bmi.inf.ethz.ch/projects/2018/graph-anno/).
