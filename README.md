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
2. Build *sdsl-lite* by `pushd external-libraries/sdsl-lite; ./install.sh $(pwd); popd`
3. go to the **build** directory `mkdir -p build && cd build`
4. compile by `cmake .. && make && ./unit_tests`

Build types: `cmake .. <arguments>` where arguments are:

- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (Debug by default)
- `-DBUILD_STATIC=ON` -- link statically (OFF by default)

### Typical workflow
1. Generate graph and uncompressed annotations (`.precise.dbg` and optionally `.wtr.dbg` files)  
`./annograph build -o <OUTPREFIX> <FLAGS> <INPUTS>`
2. Compress annotation with Bloom filters  
`./annograph build -i <OUTPREFIX> -o <BLOOMOUTPREFIX> <FLAGS> <INPUTS>`
3. Compress annotation with wavelet tries (if not done in step 1)  
`./annograph build -i <OUTPREFIX> -o <WTROUTPREFIX> --wavelet-trie <FLAGS> <INPUTS>`

#### Example
```
./annograph build -k 9 -o tiny_example ../tests/data/test_vcfparse.fa

./annograph build -i tiny_example --bloom-false-pos-prob 0.01 -o tiny_example ../tests/data/test_vcfparse.fa
./annograph map -i tiny_example TCGCGCGCTA TCGCGCGCTA TCGCGCGCTC TCGCGCGCTN TCGCGCGCTANA TCGCGCGCTC

./annograph build -i tiny_example --wavelet-trie -o tiny_example ../tests/data/test_vcfparse.fa
./annograph map --wavelet-trie -i tiny_example TCGCGCGCTA TCGCGCGCTA TCGCGCGCTC TCGCGCGCTN TCGCGCGCTANA TCGCGCGCTC
```

#### Other use cases
Constructing wavelet trie in blocks (slower, uses less RAM)  
`./annograph compress -i <OUTPREFIX> -o <WTROUTPREFIX>`

Annotation compressor query time  
`./annograph query -i <OUTPREFIX>`

Wavelet trie statistics  
`./annograph stats -i <OUTPREFIX> --wavelet-trie`

Compress wavelet tries with random column permutations  
`./annograph permutation -i <OUTPREFIX> --num-permutations <NUM_PERMS>`

## Reference
This code was used to produce the results in our paper
> [_Dynamic compression schemes for graph coloring, Bioinformatics, 2018_](https://doi.org/10.1093/bioinformatics/bty632) by Harun Mustafa, Ingo Schilken, Mikhail Karasikov, Carsten Eickhoff, Gunnar Rätsch, and André Kahles. 

The input data used for those experiments is located [here](https://public.bmi.inf.ethz.ch/projects/2018/graph-anno/).
