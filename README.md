# Hash-based colored de Bruijn graph with wavelet trie and Bloom filter color compression

**Prerequisites**
- cmake 3.6.1
- C++14
- HTSlib
- GNU GMP
- boost
- sdsl-lite

### Install
1. `git clone --recursive https://github.com/ratschlab/graph_annotation`
2. Build *sdsl-lite* by `pushd wavelet_trie/external-libraries/sdsl-lite; ./install.sh $(pwd); popd`
3. go to the **build** directory `mkdir -p annograph/build && cd annograph/build`
4. compile by `cmake .. && make && ./unit_tests`

Build types: `cmake .. <arguments>` where arguments are:

- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (Debug by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)

### Typical workflow
1. Generate graph and uncompressed annotations (`.precise.dbg` and optionally `.wtr.dbg` files)  
`./annograph build -o <OUTPREFIX> <FLAGS> <INPUTS>`
2. Compress annotation with Bloom filters  
`./annograph build -i <OUTPREFIXX> -o <BLOOMOUTPREFIX> <FLAGS> <INPUTS>`
3. Compress annotation with wavelet tries (if not done in step 1)  
`./annograph build -i <OUTPREFIXX> -o <WTROUTPREFIX> --wavelet-trie <FLAGS> <INPUTS>`

## bioRxiv preprint
This code was used to produce the results in our bioRxiv preprint
> [_Dynamic compression schemes for graph coloring_](https://www.biorxiv.org/content/early/2018/03/17/239806) by Harun Mustafa, Ingo Schilken, Mikhail Karasikov, Carsten Eickhoff, Gunnar Ratsch, and Andre Kahles. 

The input data used for those experiments is located [here](https://public.bmi.inf.ethz.ch/projects/2018/graph-anno/).
