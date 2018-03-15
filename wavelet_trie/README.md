# Wavelet trie binary matrix compression
Multithreaded wavelet trie construction library

## Prerequisites
1. sdsl-lite
2. libmaus2
3. GNU GMP

## Install
0. Compile external libraries (see README in parent directory)
1. `make`

### Usage

## General form
`<SET_ENVIRONMENT_VARIABLES> ./wtr_compress <INPUT1> <INPUT2> ... <OUTPUT>`

## With raw uncompressed binary matrix (using rawrows.dbg file from graph)
`./wtr_compress <INPUT1> <INPUT2> ... <OUTPUT>`

## With serialized annotation map (precise.dbg file from graph)
`MAP=1 ./wtr_compress <INPUT1> <INPUT2> ... <OUTPUT>`

## With CSV file (comma-separated indices for set bits in row)
`COMMA=1 ./wtr_compress <INPUT1> <INPUT2> ... <OUTPUT>`

## Other command line flags
1. `SHUF_SEED=<N>`: set seed to `N` for randomly shuffling columns
2. `NJOBS=<N>`: number of threads
3. `INDEXSET=1`: store uncompressed matrix as set of indices instead of in binary format (slower, but saves RAM)

