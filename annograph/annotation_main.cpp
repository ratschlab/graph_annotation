#include <fstream>
#include <ctime>
#include <map>
#include <memory>

#include <zlib.h>

#include "dbg_hash.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "dbg_bloom_annotator.hpp"
#include "unix_tools.hpp"

KSEQ_INIT(gzFile, gzread);


const std::vector<std::string> annots = {
  "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
  "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
};

int main(int argc, const char *argv[]) {

    Timer timer;

    // parse command line arguments and options
    std::unique_ptr<Config> config(new Config(argc, argv));

    if (config->verbose) {
        std::cout << "#############################\n"
                  << "### Welcome to AnnoGraph! ###\n"
                  << "#############################\n" << std::endl;
    }

    const auto &files = config->fname;

    if (config->identity != Config::BUILD) {
        std::cerr << "Error: Only BUILD mode is supported" << std::endl;
        exit(1);
    }

    if (config->fasta_header_delimiter.length() > 1) {
        std::cerr << "FASTA header delimiter must be at most one character" << std::endl;
        exit(1);
    }

    hash_annotate::BloomAnnotator *annotator = NULL;
    hash_annotate::PreciseHashAnnotator *precise_annotator = NULL;

    double graph_const_time = 0;
    double precise_const_time = 0;
    Timer result_timer;
    DBGHash hashing_graph(config->k);
    precise_annotator = new hash_annotate::PreciseHashAnnotator(hashing_graph);
    std::unordered_map<std::string, size_t> annot_map;
    if (!config->infbase.empty()) {
        result_timer.reset();
        std::cout << "Loading graph file" << std::endl;
        hashing_graph.load(config->infbase + ".graph.dbg");
        graph_const_time += result_timer.elapsed();
        result_timer.reset();
        precise_annotator->load(config->infbase + ".precise.dbg");
        precise_const_time += result_timer.elapsed();
    }

    if (config->bloom_fpp > -0.5
            || config->bloom_bits_per_edge > -0.5
            || config->bloom_num_hash_functions > 0) {
        if (config->bloom_fpp > -0.5) {
            // Expected FPP is set, optimize other parameters automatically
            annotator = new hash_annotate::BloomAnnotator(hashing_graph,
                                                          config->bloom_fpp,
                                                          config->verbose);
        } else {
            assert(config->bloom_bits_per_edge >= 0);
            // Experiment mode, estimate FPP given other parameters,
            // optimize the number of hash functions if it's set to zero
            annotator = new hash_annotate::BloomAnnotator(
                hashing_graph,
                config->bloom_bits_per_edge,
                config->bloom_num_hash_functions,
                config->verbose
            );
        }
    }
    if (annotator) {
        std::cout << "Bloom filter settings" << std::endl;
        std::cout << "\tBits per edge:\t" << annotator->size_factor() << std::endl;
        std::cout << "\tNum hash functions:\t" << annotator->num_hash_functions() << std::endl;
        std::cout << "\tApprox false pos prob:\t" << annotator->approx_false_positive_rate() << std::endl;
    }

    bool has_vcf = false;

    if (config->verbose)
        std::cerr << "k is " << config->k << std::endl;

    //one pass per suffix
    double file_read_time = 0;
    double bloom_const_time = 0;
    std::cout << "Start reading data and extracting k-mers..." << std::endl;

    // iterate over input files
    for (unsigned int f = 0; f < files.size(); ++f) {
        if (config->verbose) {
            std::cout << std::endl << "Parsing " << files[f] << std::endl;
        }

        if (utils::get_filetype(files[f]) == "VCF") {
            has_vcf = true;
            //READ FROM VCF
            Timer data_reading_timer;

            vcf_parser vcf;
            if (!vcf.init(config->refpath, files[f], config->k)) {
                std::cerr << "ERROR reading VCF " << files[f] << std::endl;
                exit(1);
            }
            std::cout << "Reading VCF" << std::endl;
            std::string sequence;
            std::vector<std::string> annotation;
            std::map<size_t, std::string> variants;
            data_reading_timer.reset();
            for (size_t i = 1; vcf.get_seq(annots, &sequence, annotation); ++i) {
                //doesn't cover the no annot case
                size_t cursize = annot_map.size();
                for (size_t _i = 0; _i < annotation.size(); ++_i) {
                //for (auto &annot : annotation) {
                    auto insert_annot_map = annot_map.emplace(annotation[_i], cursize);
                    if (insert_annot_map.second) {
                        cursize++;
                    }
                    if (!_i && precise_annotator && config->infbase.empty()) {
                        // assume first annotation is the reference
                        precise_annotator->make_column_prefix(insert_annot_map.first->second);
                    }
                    auto insert_annot = variants.emplace(insert_annot_map.first->second, sequence);
                    if (!insert_annot.second) {
                        insert_annot.first->second += std::string("$") + sequence;
                    }
                }
                annotation.clear();
            }
            file_read_time += data_reading_timer.elapsed();
            for (auto &variant : variants) {
                for (size_t j = 0; j < 2; ++j) {
                    if (config->infbase.empty()) {
                        data_reading_timer.reset();
                        hashing_graph.add_sequence(variant.second, true);
                        graph_const_time += data_reading_timer.elapsed();
                        if (precise_annotator) {
                            data_reading_timer.reset();
                            precise_annotator->add_sequence(variant.second, variant.first, true);
                            precise_const_time += data_reading_timer.elapsed();
                        }
                    }
                    if (annotator) {
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

                    if (config->infbase.empty()) {
                        result_timer.reset();
                        hashing_graph.add_sequence(read_stream->seq.s);
                        graph_const_time += result_timer.elapsed();
                    }

                    if (annotation.empty() && precise_annotator && config->infbase.empty()) {
                        result_timer.reset();
                        precise_annotator->add_sequence(read_stream->seq.s);
                        precise_const_time += result_timer.elapsed();
                    }

                    assert(annotation.size() <= 2);
                    for (size_t _i = 0; _i < annotation.size(); ++_i) {
                    //for (auto &annot : annotation) {
                        if (!fastq) {
                            size_t cursize = annot_map.size();
                            map_ins = annot_map.emplace(annotation[_i], cursize);
                        }
                        if (annotator) {
                            result_timer.reset();
                            annotator->add_sequence(read_stream->seq.s, map_ins.first->second,
                                    (read_stream->seq.l - hashing_graph.get_k())
                                    * (static_cast<size_t>(config->reverse) + 1));
                            bloom_const_time += result_timer.elapsed();
                        }
                        if (precise_annotator && config->infbase.empty()) {
                            result_timer.reset();
                            if (!_i && annotation.size() > 1)
                                precise_annotator->make_column_prefix(map_ins.first->second);
                            precise_annotator->add_sequence(read_stream->seq.s, map_ins.first->second);
                            precise_const_time += result_timer.elapsed();
                        }
                    }

                    if (config->reverse)
                        reverse_complement(read_stream->seq);
                    else
                        break;
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
    if (annotator)
        std::cout << "Bloom filter\t" << bloom_const_time << std::endl;
    if (precise_annotator) {
        std::cout << "Precise filter\t" << precise_const_time << std::endl;
        std::cout << "# Annotation types\t" << annot_map.size() << std::endl;
    }
    if (annotator && precise_annotator && config->bloom_test_num_kmers) {
        //Check FPP
        std::cout << "Approximating FPP...\t" << std::flush;
        timer.reset();
        annotator->test_fp_all(*precise_annotator, config->bloom_test_num_kmers, has_vcf);
        std::cout << timer.elapsed() << "sec" << std::endl;
    }

    std::cout << "Graph size(edges):\t" << hashing_graph.get_num_edges() << std::endl;

    // output and cleanup

    // graph output
    if (!config->outfbase.empty() && config->infbase.empty()) {
        hashing_graph.serialize(config->outfbase + ".graph.dbg");
    }
    if (!config->outfbase.empty() && annotator)
        annotator->serialize(config->outfbase + ".anno.dbg");

    if (!config->outfbase.empty() && precise_annotator && config->infbase.empty()) {
        std::cout << "Serializing precise annotation...\t" << std::flush;
        timer.reset();
        precise_annotator->serialize(config->outfbase + ".precise.dbg");
        precise_annotator->export_rows(config->outfbase + ".anno.rawrows.dbg");
        std::cout << timer.elapsed() << std::endl;
    }

    if (annotator)
        delete annotator;
    if (precise_annotator)
        delete precise_annotator;

    return 0;
}
