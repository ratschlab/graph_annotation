CXX=g++ -std=c++14 -fopenmp -march=native -DNPRINT -O2 -DNDEBUG -pthread -O3 #-g #-O2 -pg
HOME=/cluster/home/hmustafa
EXT_LIB=../external-libraries
CXXFLAGS=-Wall -I$(EXT_LIB)/libmaus2/include -I$(EXT_LIB)/sdsl-lite/include -L$(EXT_LIB)/libmaus2/lib -L$(EXT_LIB)/sdsl-lite/lib -I$(EXT_LIB)/gmp/include -L$(EXT_LIB)/gmp/lib
LDFLAGS=-lsdsl -lpthread -lgmp -lboost_serialization

all: wtr_compress wtr_merge wtr_test_decompression pack_sd count_unique unpack_commas

wavelet_trie.o: wavelet_trie.hpp wavelet_trie.cpp
	$(CXX) $(CXXFLAGS) -c wavelet_trie.cpp

wtr_compress.o: wtr_compress.cpp
	$(CXX) $(CXXFLAGS) -c wtr_compress.cpp

unix_tools.o: unix_tools.cpp unix_tools.hpp
	$(CXX) $(CXXFLAGS) -c unix_tools.cpp

utils.o: utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

pack_sd: pack_sd.cpp
	$(CXX) $(CXXFLAGS) -o pack_sd pack_sd.cpp $(LDFLAGS)

count_unique: count_unique.cpp
	$(CXX) $(CXXFLAGS) -o count_unique count_unique.cpp $(LDFLAGS)

unpack_commas: unpack_commas.cpp
	$(CXX) $(CXXFLAGS) -o unpack_commas unpack_commas.cpp $(LDFLAGS)

wtr_compress: wavelet_trie.o unix_tools.o wtr_compress.o utils.o
	$(CXX) $(CXXFLAGS) -o wtr_compress wavelet_trie.o wtr_compress.o unix_tools.o utils.o $(LDFLAGS)

wtr_merge: wavelet_trie.o wtr_merge.cpp utils.o
	$(CXX) $(CXXFLAGS) -o wtr_merge wavelet_trie.o utils.o wtr_merge.cpp $(LDFLAGS)

wtr_test_decompression: wavelet_trie.o wtr_test_decompression.cpp utils.o
	$(CXX) $(CXXFLAGS) -o wtr_test_decompression wavelet_trie.o utils.o wtr_test_decompression.cpp $(LDFLAGS)

clean:
	rm pack_sd count_unique unpack_commas wtr_compress wtr_merge wtr_test_decompression *.o
