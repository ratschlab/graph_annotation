#include <fstream>
#include <sstream>
#include <ctime>
#include <map>
#include <memory>
#include <algorithm>

#include <zlib.h>

#include "dbg_hash.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "dbg_bloom_annotator.hpp"
#include "wavelet_trie_annotator.hpp"
#include "unix_tools.hpp"

KSEQ_INIT(gzFile, gzread);


template <class Graph, class Map>
void annotate_kmers(const std::vector<std::string> &kmers, const Graph &graph,
                                                           const Map &get_coloring) {
    for (const auto &kmer : kmers) {
        std::cout << kmer << "\t";
        if (kmer.size() != graph.get_k() + 1) {
            std::cout << "\n";
            std::cerr << "Error: Wrong k-mer size (" << kmer.size()
                                                     << " instead of "
                                                     << graph.get_k() + 1
                                                     << ")\n";
            continue;
        }
        std::vector<uint64_t> coloring;
        auto kmer_index = graph.map_kmer(kmer);
        if (kmer_index >= graph.first_edge()
                && kmer_index <= graph.last_edge()) {
            coloring = get_coloring(kmer_index);
        }

        if (coloring.size())
            std::cout << coloring[0];

        for (size_t i = 1; i < coloring.size(); ++i) {
            std::cout << "," << coloring[i];
        }
        std::cout << "\n";
    }
}

std::vector<size_t> wavelet_trie_test_permutations(
        const DBGHash &graph,
        const hash_annotate::PreciseHashAnnotator &precise,
        size_t num_perm,
        size_t p,
        bool verbose = false) {
    if (verbose) {
        std::cout << "Testing permutations: " << precise.num_prefix_columns() << " prefix columns" << std::endl;
    }
    std::vector<size_t> sizes(num_perm + 1);
    utils::ThreadPool thread_pool(p);
    //original permutation
    thread_pool.enqueue([&]() {
        std::ostringstream sout;
        annotate::WaveletTrieAnnotator wtr(
                precise,
                graph,
                1);
        sizes[0] = wtr.serialize(sout);
    });

    std::srand(42);

    //shuffles
    for (size_t _i = 1; _i < sizes.size(); ++_i) {
        std::vector<size_t> indices(precise.num_columns());
        std::iota(indices.begin(), indices.end(), 0);
        std::random_shuffle(indices.begin(), indices.end());
        std::unordered_map<size_t, size_t> permut_map;
        for (size_t i = 0; i < indices.size(); ++i) {
            permut_map.emplace(i, indices[i]);
        }
        thread_pool.enqueue(
            [&precise, &graph, &sizes](size_t num_perm,
                                       std::unordered_map<size_t, size_t> permut_map) {
                std::ostringstream sout;
                annotate::WaveletTrieAnnotator wtr(
                            precise,
                            graph,
                            1,
                            std::move(permut_map));
                sizes[num_perm] = wtr.serialize(sout);
            },
        _i, std::move(permut_map));
    }
    thread_pool.join();
    for (auto &size : sizes) {
        std::cout << size << std::endl;
    }
    return sizes;
}

template <class Callback>
void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::vector<std::string> *annotation,
                            Callback callback, Timer *timer = NULL) {

    //TODO: make this a configurable option
    //default list of tokens to extract as annotations
    //TODO: extract these guys directly from vcf parsed
    const std::vector<std::string> annots = {
      "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
      "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
    };

    vcf_parser vcf;
    if (!vcf.init(ref_filename, filename, k)) {
        std::cerr << "ERROR reading VCF " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Extracting sequences from file " << filename << std::endl;
    }
    size_t seq_count = 0;
    while (vcf.get_seq(annots, annotation)) {
        callback(vcf.seq, annotation);
        seq_count++;
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
}

template <class Callback>
void read_fasta_file_critical(const std::string &filename,
                              Callback callback,
                              bool with_reverse = false,
                              Timer *timer = NULL) {
    size_t seq_count = 0;

    gzFile input_p = gzopen(filename.c_str(), "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR no such file " << filename << std::endl;
        exit(1);
    }
    //TODO: handle read_stream->qual
    kseq_t *read_stream = kseq_init(input_p);
    if (read_stream == NULL) {
        std::cerr << "ERROR while opening input file " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Start extracting sequences from file " << filename << std::endl;
    }
    while (kseq_read(read_stream) >= 0) {
        callback(read_stream);
        if (with_reverse) {
            reverse_complement(read_stream->seq);
            callback(read_stream);
        }
        seq_count++;
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
    kseq_destroy(read_stream);
    gzclose(input_p);
}

int main(int argc, const char *argv[]) {

    // parse command line arguments and options
    std::unique_ptr<Config> config(new Config(argc, argv));

    if (config->verbose) {
        std::cout << "#############################\n"
                  << "### Welcome to AnnoGraph! ###\n"
                  << "#############################\n" << std::endl;
    }

    const auto &files = config->fname;

    std::unique_ptr<hash_annotate::BloomAnnotator> annotator;
    std::unique_ptr<hash_annotate::PreciseHashAnnotator> precise_annotator;
    std::unique_ptr<annotate::WaveletTrieAnnotator> wt_annotator;

    if (config->identity == Config::BUILD) {

        double graph_const_time = 0;
        double precise_const_time = 0;
        Timer result_timer;
        DBGHash hashing_graph(config->k);
        wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
        precise_annotator.reset(new hash_annotate::PreciseHashAnnotator(hashing_graph));
        std::unordered_map<std::string, size_t> annot_map;
        if (!config->infbase.empty()) {
            result_timer.reset();
            std::cout << "Loading graph file" << std::endl;
            hashing_graph.load(config->infbase + ".graph.dbg");
            graph_const_time += result_timer.elapsed();
        }

        if (config->bloom_fpp > -0.5
                || config->bloom_bits_per_edge > -0.5
                || config->bloom_num_hash_functions > 0) {
            if (config->bloom_fpp > -0.5) {
                // Expected FPP is set, optimize other parameters automatically
                annotator.reset(new hash_annotate::BloomAnnotator(hashing_graph,
                                                                  config->bloom_fpp,
                                                                  config->verbose));
            } else {
                assert(config->bloom_bits_per_edge >= 0);
                // Experiment mode, estimate FPP given other parameters,
                // optimize the number of hash functions if it's set to zero
                annotator.reset(new hash_annotate::BloomAnnotator(
                    hashing_graph,
                    config->bloom_bits_per_edge,
                    config->bloom_num_hash_functions,
                    config->verbose
                ));
            }
        }
        if (annotator.get()) {
            std::cout << "Bloom filter settings" << std::endl;
            std::cout << "\tBits per edge:\t" << annotator->size_factor() << std::endl;
            std::cout << "\tNum hash functions:\t" << annotator->num_hash_functions() << std::endl;
            std::cout << "\tApprox false pos prob:\t" << annotator->approx_false_positive_rate() << std::endl;
        }
        if (!config->infbase.empty()) {
            if (annotator.get()) {
                result_timer.reset();
                wt_annotator->load(config->infbase + ".wtr.dbg");
                //precise_annotator->load(config->infbase + ".precise.dbg");
                precise_const_time += result_timer.elapsed();
            } else {
                precise_annotator.reset();
                wt_annotator.reset();
            }
        }

        bool has_vcf = false;

        if (config->verbose) {
            std::cerr << "k is " << hashing_graph.get_k() << std::endl;
        }

        //one pass per suffix
        double file_read_time = 0;
        double bloom_const_time = 0;
        std::cout << "Start reading data and extracting k-mers..." << std::endl;
        Timer timer;
        timer.reset();

        // iterate over input files
        for (unsigned int f = 0; f < files.size(); ++f) {
            if (!annotator.get() && !precise_annotator.get()) {
            //if (!annotator.get() && !wt_annotator.get())
                //break;
            //if (!annotator.get() && (wt_annotator.get() || precise_annotator.get())) {
                break;
            }
            if (config->verbose) {
                std::cout << std::endl << "Parsing " << files[f] << std::endl;
            }

            if (utils::get_filetype(files[f]) == "VCF") {
                has_vcf = true;
                //READ FROM VCF
                Timer data_reading_timer;
                std::vector<std::string> annotation;
                std::map<size_t, std::string> variants;

                read_vcf_file_critical(files[f],
                                       config->refpath,
                                       hashing_graph.get_k(),
                                       &annotation,
                    [&](std::string &seq, std::vector<std::string> *annotation) {
                        //doesn't cover the no annot case
                        size_t cursize = annot_map.size();
                        for (size_t _i = 0; _i < annotation->size(); ++_i) {
                            auto insert_annot_map = annot_map.emplace(
                                    annotation->operator[](_i),
                                    cursize);
                            if (insert_annot_map.second) {
                                cursize++;
                            }
                            if (!_i
                                && precise_annotator.get()
                                //&& wt_annotator.get()
                                && config->infbase.empty()) {
                                precise_annotator->make_column_prefix(
                                        insert_annot_map.first->second);
                            }
                            auto insert_annot = variants.emplace(
                                    insert_annot_map.first->second,
                                    seq);
                            if (!insert_annot.second) {
                                insert_annot.first->second += std::string("$") + seq;
                            }
                        }
                        annotation->clear();
                    }, config->verbose ? &data_reading_timer : NULL);

                file_read_time += data_reading_timer.elapsed();
                for (auto &variant : variants) {
                    for (size_t j = 0; j < 2; ++j) {
                        if (config->infbase.empty()) {
                            data_reading_timer.reset();
                            hashing_graph.add_sequence(variant.second, true);
                            graph_const_time += data_reading_timer.elapsed();
                            if (precise_annotator.get()) {
                            //if (wt_annotator.get()) {
                                data_reading_timer.reset();
                                //wt_annotator->add_sequence(variant.second, variant.first, true);
                                precise_annotator->add_sequence(variant.second, variant.first, true);
                                precise_const_time += data_reading_timer.elapsed();
                            }
                        }
                        if (annotator.get()) {
                            data_reading_timer.reset();
                            annotator->add_sequence(variant.second, variant.first,
                                        (variant.second.length() - hashing_graph.get_k())
                                            * (static_cast<size_t>(config->reverse) + 1));
                            bloom_const_time += data_reading_timer.elapsed();
                        }
                        if (config->reverse) {
                            reverse_complement(variant.second.begin(), variant.second.end());
                        } else {
                            break;
                        }
                    }
                }
            } else if (utils::get_filetype(files[f]) == "FASTA"
                        || utils::get_filetype(files[f]) == "FASTQ") {
                result_timer.reset();
                bool fastq = utils::get_filetype(files[f]) == "FASTQ";
                std::pair<std::unordered_map<std::string, size_t>::iterator, bool> map_ins;
                std::vector<std::string> annotation;
                if (fastq) {
                    //TODO annotations for reference genome
                    size_t cursize = annot_map.size();
                    map_ins = annot_map.emplace(files[f], cursize);
                    annotation.push_back(files[f]);
                }

                read_fasta_file_critical(files[f],
                    [&](kseq_t *read_stream) {
                        file_read_time += result_timer.elapsed();
                        if (!fastq) {
                            if (config->verbose) {
                                std::cout << "Parsing " << read_stream->name.s << "\n";
                            }
                            annotation.clear();
                            char *sep = NULL;
                            if (!config->fasta_header_delimiter.empty())
                                sep = strchr(read_stream->name.s, config->fasta_header_delimiter[0]);
                            if (sep) {
                                annotation.emplace_back(sep);
                                annotation.emplace_back(std::string(read_stream->name.s, sep));
                            } else {
                                annotation.emplace_back(read_stream->name.s);
                            }
                        }

                        if (config->infbase.empty()) {
                            result_timer.reset();
                            hashing_graph.add_sequence(read_stream->seq.s);
                            graph_const_time += result_timer.elapsed();
                        }

                        if (annotation.empty() && precise_annotator.get() && config->infbase.empty()) {
                        //if (annotation.empty() && wt_annotator.get() && config->infbase.empty()) {
                            result_timer.reset();
                            //wt_annotator->add_sequence(read_stream->seq.s);
                            precise_annotator->add_sequence(read_stream->seq.s);
                            precise_const_time += result_timer.elapsed();
                        }

                        assert(annotation.size() <= 2);
                        for (size_t _i = 0; _i < annotation.size(); ++_i) {
                            if (!fastq) {
                                size_t cursize = annot_map.size();
                                map_ins = annot_map.emplace(annotation[_i], cursize);
                            }
                            if (annotator.get()) {
                                result_timer.reset();
                                annotator->add_sequence(read_stream->seq.s, map_ins.first->second,
                                        (read_stream->seq.l - hashing_graph.get_k())
                                        * (static_cast<size_t>(config->reverse) + 1));
                                bloom_const_time += result_timer.elapsed();
                            }
                            if (precise_annotator.get() && config->infbase.empty()) {
                            //if (wt_annotator.get() && config->infbase.empty()) {
                                result_timer.reset();
                                if (!_i && annotation.size() > 1) {
                                    precise_annotator->make_column_prefix(map_ins.first->second);
                                }
                                //wt_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                                precise_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                                precise_const_time += result_timer.elapsed();
                            }
                        }
                        result_timer.reset();
                    },
                    config->reverse,
                    &result_timer);
            } else {
                std::cerr << "ERROR: Filetype unknown for file "
                          << files[f] << std::endl;
                exit(1);
            }
            //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
        }
        get_RAM();
        std::cout << "Reading data finished\t" << timer.elapsed() << "sec" << std::endl;

        //Runtime stats
        std::cout << "Runtime statistics" << std::endl;
        std::cout << "File reading\t" << file_read_time << std::endl;
        std::cout << "Graph construction\t" << graph_const_time << std::endl;
        if (annotator.get())
            std::cout << "Bloom filter\t" << bloom_const_time << std::endl;
        if (precise_annotator.get()) {
            std::cout << "Index set\t" << precise_const_time << std::endl;
            std::cout << "# colors\t" << precise_annotator->num_columns() << std::endl;
            std::cout << "# class indicators\t" << precise_annotator->num_prefix_columns() << std::endl;
        }
        //if (annotator.get() && precise_annotator.get() && config->bloom_test_num_kmers) {
        if (annotator.get() && (wt_annotator.get() || precise_annotator.get()) && config->bloom_test_num_kmers) {
            //Check FPP
            std::cout << "Approximating FPP...\t" << std::flush;
            timer.reset();
            /*
            if (wt_annotator.get()) {
                annotator->test_fp_all(*wt_annotator, config->bloom_test_num_kmers, has_vcf);
            } else if (precise_annotator.get()) {
                annotator->test_fp_all(*precise_annotator, config->bloom_test_num_kmers, has_vcf);
            }
            */
            annotator->test_fp_all(*precise_annotator, config->bloom_test_num_kmers, has_vcf);
            std::cout << timer.elapsed() << "sec" << std::endl;
        }
        // output and cleanup

        // graph output
        if (!config->outfbase.empty() && config->infbase.empty()) {
            std::cout << "Serializing hash graph\t" << std::flush;
            std::cout << hashing_graph.serialize(config->outfbase + ".graph.dbg")
                      << " bytes" << std::endl;
        }
        if (!config->outfbase.empty() && annotator.get()) {
            std::cout << "Serializing bloom filters\t" << std::flush;
            std::cout << annotator->serialize(config->outfbase + ".anno.dbg")
                      << " bytes" << std::endl;
        }

        if (config->wavelet_trie) {
            std::cout << "Computing wavelet trie\t" << std::flush;
            timer.reset();
            if (precise_annotator.get()) {
                wt_annotator.reset(new annotate::WaveletTrieAnnotator(
                    *precise_annotator, hashing_graph, config->p
                ));
            } else {
                std::ifstream in(config->infbase + ".precise.dbg");
                if (!in.good()) {
                    std::cerr << "ERROR: corrupt precise annotator. Please reconstruct it."
                              << std::endl;
                    exit(1);
                }
                wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
                wt_annotator->load_from_precise_file(in, config->p);
            }
            std::cout << wt_annotator->serialize(config->outfbase + ".wtr.dbg")
                      << " bytes\t"
                      << timer.elapsed() << " s\t"
                      << config->p << " threads" << std::endl;
        }

        if (!config->outfbase.empty() && precise_annotator.get() && config->infbase.empty()) {
            std::cout << "Serializing index set\t" << std::flush;
            timer.reset();
            //precise_annotator->export_rows(config->outfbase + ".anno.rawrows.dbg");
            std::cout << precise_annotator->serialize(config->outfbase + ".precise.dbg")
                      << " bytes\t"
                      << timer.elapsed() << " s" << std::endl;
        }

        std::cout << "Annotation matrix size\t"
                  << hashing_graph.get_num_edges()
                  << " x "
                  << (annotator.get()
                            ? annotator->num_columns()
                            : (wt_annotator.get()
                                    ? wt_annotator->num_columns()
                                    : 0))
                  << std::endl;

    } else if (config->identity == Config::UPDATE) {

        Timer timer;
        DBGHash hashing_graph(0);
        if (!hashing_graph.load(config->infbase + ".graph.dbg")) {
            std::cerr << "Error: Graph loading failed for "
                      << config->infbase + ".graph.dbg" << std::endl;
            exit(1);
        }
        if (config->verbose)
            std::cout << "k is " << hashing_graph.get_k() << std::endl;

        if (config->verbose) {
            std::cout << "Graph loading: " << timer.elapsed() << "sec" << std::endl;
        }
        timer.reset();

        double graph_const_time = 0;
        double precise_const_time = 0;
        Timer result_timer;
        // precise_annotator.reset(new hash_annotate::PreciseHashAnnotator(hashing_graph));
        std::unordered_map<std::string, size_t> annot_map;

        if (config->wavelet_trie) {
            wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
            if (!wt_annotator->load(config->infbase + ".wtr.dbg")) {
                std::cerr << "Error: Can't load Wavelet Trie annotation from "
                          << config->infbase + ".wtr.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Wavelet Trie loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();
        } else {
            annotator.reset(new hash_annotate::BloomAnnotator(hashing_graph, 0.5));
            if (!annotator->load(config->infbase + ".anno.dbg")) {
                std::cerr << "Error: Can't load Bloom filter annotation from "
                          << config->infbase + ".anno.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Bloom filter loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            if (annotator.get() && config->verbose) {
                std::cout << "Bloom filter settings" << std::endl;
                std::cout << "\tBits per edge:\t" << annotator->size_factor() << std::endl;
                std::cout << "\tNum hash functions:\t" << annotator->num_hash_functions() << std::endl;
                std::cout << "\tApprox false pos prob:\t" << annotator->approx_false_positive_rate() << std::endl;
            }
        }

        bool has_vcf = false;

        //one pass per suffix
        double file_read_time = 0;
        double bloom_const_time = 0;
        std::cout << "Start reading data and extracting k-mers..." << std::endl;
        timer.reset();

        // iterate over input files
        for (size_t f = 0; f < files.size(); ++f) {
            if (!annotator.get() && !wt_annotator.get())
                break;

            if (config->verbose)
                std::cout << std::endl << "Parsing " << files[f] << std::endl;

            if (utils::get_filetype(files[f]) == "VCF") {
                has_vcf = true;
                //READ FROM VCF
                Timer data_reading_timer;

                std::vector<std::string> annotation;
                std::map<size_t, std::string> variants;

                read_vcf_file_critical(files[f],
                                       config->refpath,
                                       hashing_graph.get_k(),
                                       &annotation,
                    [&](std::string &seq, std::vector<std::string> *annotation) {
                        //doesn't cover the no annot case
                        size_t cursize = annot_map.size();
                        for (size_t _i = 0; _i < annotation->size(); ++_i) {
                            auto insert_annot_map = annot_map.emplace(
                                    annotation->operator[](_i),
                                    cursize);
                            if (insert_annot_map.second) {
                                cursize++;
                            }
                            if (!_i
                                && wt_annotator.get()
                                && config->infbase.empty()) {
                                //TODO
                                //precise_annotator->make_column_prefix(
                                //        insert_annot_map.first->second);
                            }
                            auto insert_annot = variants.emplace(
                                    insert_annot_map.first->second,
                                    seq);
                            if (!insert_annot.second) {
                                insert_annot.first->second += std::string("$") + seq;
                            }
                        }
                        annotation->clear();
                    }, config->verbose ? &data_reading_timer : NULL);

                file_read_time += data_reading_timer.elapsed();
                for (auto &variant : variants) {
                    for (size_t j = 0; j < 2; ++j) {

                        data_reading_timer.reset();
                        hashing_graph.add_sequence(variant.second, true);
                        graph_const_time += data_reading_timer.elapsed();
                        if (wt_annotator.get()) {
                            data_reading_timer.reset();
                            wt_annotator->add_sequence(variant.second, variant.first, true);
                            precise_const_time += data_reading_timer.elapsed();
                        }

                        if (annotator.get()) {
                            data_reading_timer.reset();
                            annotator->add_sequence(variant.second, variant.first,
                                        (variant.second.length() - hashing_graph.get_k())
                                            * (static_cast<size_t>(config->reverse) + 1));
                            bloom_const_time += data_reading_timer.elapsed();
                        }
                        if (config->reverse) {
                            reverse_complement(variant.second.begin(), variant.second.end());
                        } else {
                            break;
                        }
                    }
                }
            } else if (utils::get_filetype(files[f]) == "FASTA"
                        || utils::get_filetype(files[f]) == "FASTQ") {
                // open stream
                gzFile input_p = gzopen(files[f].c_str(), "r");
                if (input_p == Z_NULL) {
                    std::cerr << "ERROR no such file " << files[f] << std::endl;
                    exit(1);
                }

                //TODO: handle read_stream->qual
                kseq_t *read_stream = kseq_init(input_p);
                if (read_stream == NULL) {
                    std::cerr << "ERROR while opening input file "
                              << files[f] << std::endl;
                    exit(1);
                }
                result_timer.reset();
                bool fastq = utils::get_filetype(files[f]) == "FASTQ";
                std::pair<std::unordered_map<std::string, size_t>::iterator, bool> map_ins;
                std::vector<std::string> annotation;
                if (fastq) {
                    //TODO annotations for reference genome
                    size_t cursize = annot_map.size();
                    map_ins = annot_map.emplace(files[f], cursize);
                    annotation.push_back(files[f]);
                }

                while (kseq_read(read_stream) >= 0) {
                    file_read_time += result_timer.elapsed();
                    if (!fastq) {
                        if (config->verbose) {
                            std::cout << "Parsing " << read_stream->name.s << "\n";
                        }
                        annotation.clear();
                        char *sep = NULL;
                        if (!config->fasta_header_delimiter.empty())
                            sep = strchr(read_stream->name.s, config->fasta_header_delimiter[0]);
                        if (sep) {
                            annotation.emplace_back(sep);
                            annotation.emplace_back(std::string(read_stream->name.s, sep));
                        } else {
                            annotation.emplace_back(read_stream->name.s);
                        }
                    }
                    for (size_t _j = 0; _j < 2; ++_j) {
                        result_timer.reset();
                        hashing_graph.add_sequence(read_stream->seq.s);
                        graph_const_time += result_timer.elapsed();

                        if (annotation.empty() && wt_annotator.get()) {
                            result_timer.reset();
                            wt_annotator->add_sequence(read_stream->seq.s);
                            precise_const_time += result_timer.elapsed();
                        }

                        assert(annotation.size() <= 2);
                        for (size_t _i = 0; _i < annotation.size(); ++_i) {
                        //for (auto &annot : annotation) {
                            if (!fastq) {
                                size_t cursize = annot_map.size();
                                map_ins = annot_map.emplace(annotation[_i], cursize);
                            }
                            if (annotator.get()) {
                                result_timer.reset();
                                annotator->add_sequence(read_stream->seq.s, map_ins.first->second,
                                        (read_stream->seq.l - hashing_graph.get_k())
                                        * (static_cast<size_t>(config->reverse) + 1));
                                bloom_const_time += result_timer.elapsed();
                            }
                            if (wt_annotator.get()) {
                                result_timer.reset();
                                if (!_i && annotation.size() > 1) {
                                    //TODO
                                    //precise_annotator->make_column_prefix(map_ins.first->second);
                                }
                                wt_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                                precise_const_time += result_timer.elapsed();
                            }
                        }

                        if (config->reverse) {
                            reverse_complement(read_stream->seq);
                        } else {
                            break;
                        }
                    }
                    result_timer.reset();
                }
                kseq_destroy(read_stream);
                gzclose(input_p);
            } else {
                std::cerr << "ERROR: Filetype unknown for file "
                          << files[f] << std::endl;
                exit(1);
            }
            //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
        }
        get_RAM();
        std::cout << "Reading data finished\t" << timer.elapsed() << "sec" << std::endl;

        //Runtime stats
        std::cout << "Runtime statistics" << std::endl;
        std::cout << "File reading\t" << file_read_time << std::endl;
        std::cout << "Graph construction\t" << graph_const_time << std::endl;
        if (annotator.get())
            std::cout << "Bloom filter\t" << bloom_const_time << std::endl;
        /*
        if (precise_annotator.get()) {
            std::cout << "Index set\t" << precise_const_time << std::endl;
            std::cout << "# colors\t" << precise_annotator->num_columns() << std::endl;
            std::cout << "# class indicators\t" << precise_annotator->num_prefix_columns() << std::endl;
        }
        */
        if (annotator.get() && wt_annotator.get() && config->bloom_test_num_kmers) {
            //Check FPP
            std::cout << "Approximating FPP...\t" << std::flush;
            timer.reset();
            annotator->test_fp_all(*wt_annotator, config->bloom_test_num_kmers, has_vcf);
            std::cout << timer.elapsed() << "sec" << std::endl;
        }
        // output and cleanup

        // graph output
        if (!config->outfbase.empty()) {
            std::cout << "Serializing hash graph\t" << std::flush;
            std::cout << hashing_graph.serialize(config->outfbase + ".graph.dbg")
                      << " bytes" << std::endl;
        }
        if (!config->outfbase.empty() && annotator.get()) {
            std::cout << "Serializing bloom filters\t" << std::flush;
            std::cout << annotator->serialize(config->outfbase + ".anno.dbg")
                      << " bytes" << std::endl;
        }

        if (config->wavelet_trie) {
            std::cout << "Wavelet trie update time\t" << std::flush;
            std::cout << precise_const_time << "sec" << std::endl;
        }

        if (!config->outfbase.empty() && config->wavelet_trie) {
            /*
            std::cout << "Computing wavelet trie\t" << std::flush;
            timer.reset();
            if (precise_annotator.get()) {
                wt_annotator.reset(new annotate::WaveletTrieAnnotator(
                    *precise_annotator, hashing_graph, config->p
                ));
            } else {
                std::ifstream in(config->infbase + ".precise.dbg");
                if (!in.good()) {
                    std::cerr << "ERROR: corrupt precise annotator. Please reconstruct it."
                              << std::endl;
                    exit(1);
                }
                wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
                wt_annotator->load_from_precise_file(in, config->p);
            }
            */
            std::cout << wt_annotator->serialize(config->outfbase + ".wtr.dbg")
                      << " bytes\t"
                      << timer.elapsed() << " s" << std::endl;
        }

        /*
        if (!config->outfbase.empty() && precise_annotator.get()) {
            std::cout << "Serializing index set\t" << std::flush;
            timer.reset();
            //precise_annotator->export_rows(config->outfbase + ".anno.rawrows.dbg");
            std::cout << precise_annotator->serialize(config->outfbase + ".precise.dbg")
                      << " bytes\t"
                      << timer.elapsed() << " s" << std::endl;
        }
        */

        std::cout << "Annotation matrix size:\t"
                  << hashing_graph.get_num_edges()
                  << " x "
                  << (annotator.get()
                            ? annotator->num_columns()
                            : (wt_annotator.get()
                                    ? wt_annotator->num_columns()
                                    : 0))
                  << std::endl;

    } else if (config->identity == Config::MAP) {

        Timer timer;
        DBGHash hashing_graph(0);
        if (!hashing_graph.load(config->infbase + ".graph.dbg")) {
            std::cerr << "Error: Graph loading failed for "
                      << config->infbase + ".graph.dbg" << std::endl;
            exit(1);
        }
        if (config->verbose) {
            std::cout << "Graph loading: " << timer.elapsed() << "sec" << std::endl;
        }
        timer.reset();

        if (config->wavelet_trie) {
            wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
            if (!wt_annotator->load(config->infbase + ".wtr.dbg")) {
                std::cerr << "Error: Can't load Wavelet Trie annotation from "
                          << config->infbase + ".wtr.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Wavelet Trie loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            annotate_kmers(files, hashing_graph, [&](uint64_t kmer_index) {
                return wt_annotator->annotate_edge(kmer_index);
            });
        } else {
            annotator.reset(new hash_annotate::BloomAnnotator(hashing_graph, 0.5));
            if (!annotator->load(config->infbase + ".anno.dbg")) {
                std::cerr << "Error: Can't load Bloom filter annotation from "
                          << config->infbase + ".anno.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Bloom filter loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            annotate_kmers(files, hashing_graph, [&](uint64_t kmer_index) {
                return annotator->get_annotation_corrected(kmer_index, true, 50);
            });
        }

        if (config->verbose) {
            std::cout << "Mapping finished in " << timer.elapsed() << "sec" << std::endl;
        }
    } else if (config->identity == Config::PERMUTATION) {
        Timer timer;
        DBGHash hashing_graph(0);
        if (!hashing_graph.load(config->infbase + ".graph.dbg")) {
            std::cerr << "Error: Graph loading failed for "
                      << config->infbase + ".graph.dbg" << std::endl;
            exit(1);
        }
        if (config->verbose) {
            std::cout << "Graph loading: " << timer.elapsed() << "sec" << std::endl;
        }
        timer.reset();

        precise_annotator.reset(new hash_annotate::PreciseHashAnnotator(hashing_graph));
        precise_annotator->load(config->infbase + ".precise.dbg");
        if (config->verbose) {
            std::cout << "Annotation loading: " << timer.elapsed() << "sec" << std::endl;
        }

        auto sizes = wavelet_trie_test_permutations(
                hashing_graph,
                *precise_annotator,
                config->num_permutations,
                config->p,
                config->verbose);
        if (config->verbose) {
            std::cout << "Permutations: " << timer.elapsed() << "sec" << std::endl;
        }

    } else if (config->identity == Config::QUERY) {
        Timer timer;
        DBGHash hashing_graph(0);
        std::cout << config->infbase << "\n";
        timer.reset();
        if (config->wavelet_trie) {
            wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
            if (!wt_annotator->load(config->infbase + ".wtr.dbg")) {
                std::cerr << "Error: Can't load Wavelet Trie annotation from "
                          << config->infbase + ".wtr.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Wavelet Trie loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            for (size_t i = 0; i < std::min(static_cast<size_t>(1000000), wt_annotator->size()); ++i) {
                wt_annotator->annotate_edge(i);
            }

        } else {
            if (!hashing_graph.load(config->infbase + ".graph.dbg")) {
                std::cerr << "Error: Graph loading failed for "
                          << config->infbase + ".graph.dbg" << std::endl;
                exit(1);
            }
            if (config->verbose) {
                std::cout << "Graph loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            /*
            precise_annotator.reset(new hash_annotate::PreciseHashAnnotator(hashing_graph));
            precise_annotator->load(config->infbase + ".precise.dbg");
            if (config->verbose) {
                std::cout << "Annotation loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();
            */

            annotator.reset(new hash_annotate::BloomAnnotator(hashing_graph, 0.5));
            if (!annotator->load(config->infbase + ".anno.dbg")) {
                std::cerr << "Error: Can't load Bloom filter annotation from "
                          << config->infbase + ".anno.dbg" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Bloom filter loading: " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            for (auto i = hashing_graph.first_edge(); i <= hashing_graph.last_edge(); ++i) {
                annotator->get_annotation_corrected(i, true, 0);
            }
        }

        std::cout << "Query: " << timer.elapsed() << "sec" << std::endl;
    } else if (config->identity == Config::STATS) {
        DBGHash hashing_graph(0);
        std::cout << config->infbase << "\n";
        if (config->wavelet_trie) {
            wt_annotator.reset(new annotate::WaveletTrieAnnotator(hashing_graph, config->p));
            if (!wt_annotator->load(config->infbase + ".wtr.dbg")) {
                std::cerr << "Error: Can't load Wavelet Trie annotation from "
                          << config->infbase + ".wtr.dbg" << std::endl;
                exit(1);
            }
            auto stats = wt_annotator->stats();
            std::cout << "Annotation matrix size\t"
                      << std::get<0>(stats) << " x " << std::get<1>(stats) << std::endl;
            std::cout << "# set bits\t"
                      << std::get<2>(stats) << "\t"
                      << static_cast<double>(std::get<2>(stats) * 100) / (std::get<0>(stats) * std::get<1>(stats)) << " %" << std::endl;
            std::cout << "# unique edge colorings\t"
                      << std::get<3>(stats) << std::endl;
        }
    } else {
        std::cerr << "Error: Only BUILD, MAP, and PERMUTATION modes are currently supported" << std::endl;
        exit(1);
    }
    return 0;
}
