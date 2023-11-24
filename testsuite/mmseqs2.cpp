#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/constexpr_for.hpp"
#include "alignment/blosum62.hpp"
#include "filter/MMseqs2It.hpp"

#ifndef MAX_SUBQGRAM_LENGTH
#define MAX_SUBQGRAM_LENGTH 3
#endif

//#undef __SSSE3__
#define ALLSEED

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false,
       with_simd = false,
       list_seeds = false,
       show = false,
       short_header_option = false,
       mmseqs = false,
       correct = false;
  std::string seeds = "0", sensitivity = "1",correct_ratio = "1",num_threads = "1";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"run MMseqs2 on a query & a target database");
    options.set_width(80);
    options.custom_help(std::string("[options] query_db target_db"));
    options.set_tab_expansion();
    options.add_options()
        ("e,short_header", "show header up to and excluding the first blank",
        cxxopts::value<bool>(short_header_option)->default_value("false"))
        ("l,list_seeds", "output available seed list",
        cxxopts::value<bool>(list_seeds)->default_value("false"))
        ("d,seeds", "choose seeds to compute",
        cxxopts::value<std::string>(seeds)->default_value("0"))
        ("w,with_simd", "compute with SIMD",
        cxxopts::value<bool>(with_simd)->default_value("false"))
        ("c,correct", "compute with local score correction",
        cxxopts::value<bool>(correct)->default_value("false"))
        ("r,correct_ratio", "compute with local score correction",
        cxxopts::value<std::string>(correct_ratio)->default_value("1"))
        ("v,sensitivity", "set sensitivity",
        cxxopts::value<std::string>(sensitivity)->default_value("1")) 
        ("t,threads", "set number of threads",
        cxxopts::value<std::string>(num_threads)->default_value("1")) 
        ("s,show", "output matches",
        cxxopts::value<bool>(show)->default_value("false"))
        ("m,mmseqs", "mmseqs mode using unsorted q-grams only",
        cxxopts::value<bool>(mmseqs)->default_value("false"))
        ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      } else
      {
        const std::vector<std::string>& unmatched_args = result.unmatched();
        if(unmatched_args.size() != 2){
          throw cxxopts::OptionException("Exactly 2 dbs are needed.");
        }
        for (size_t idx = 0; idx < unmatched_args.size(); idx++)
        {
          inputfiles.push_back(unmatched_args[idx]);
        }
      }
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  bool short_header_option_is_set(void) const noexcept
  {
    return short_header_option;
  }
  bool with_simd_is_set(void) const noexcept
  {
    return with_simd;
  }
  bool list_seed_option_is_set(void) const noexcept
  {
    return list_seeds;
  }
  bool show_option_is_set(void) const noexcept
  {
    return show;
  }
  const std::string &seeds_get(void) const noexcept
  {
    return seeds;
  }
  const std::string &sensitivity_get(void) const noexcept
  {
    return sensitivity;
  }
  const std::string &correct_ratio_get(void) const noexcept
  {
    return correct_ratio;
  }
  const std::string &num_threads_get(void) const noexcept
  {
    return num_threads;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  bool mmseqs_mode_is_set(void) const noexcept
  {
    return mmseqs;
  }
  bool local_score_correction_is_set(void) const noexcept
  {
    return correct;
  }
};

#ifdef SEEDWEIGHT5
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  59UL /* 1, 5, 6 111011, MMseq2_proteins_1 */,
  107UL /* 2, 5, 7 1101011, MMseq2_proteins_2 */,
  3205UL /* 3, 5, 12 110010000101, MMseq2_proteins_3 */,
};
#endif


#ifdef SEEDWEIGHT4
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  29UL /* 0, 4, 5 11101, MMseq2_proteins_0 */,
};
#endif

#ifdef SEEDWEIGHT6
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  237UL /* 4, 6, 8 11101101, MMseq2_proteins_4 */,
};
#endif

#ifdef ALLSEED
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  29UL /* 0, 4, 5 11101, MMseq2_proteins_0 */,
  59UL /* 1, 5, 6 111011, MMseq2_proteins_1 */,
  107UL /* 2, 5, 7 1101011, MMseq2_proteins_2 */,
  3205UL /* 3, 5, 12 110010000101, MMseq2_proteins_3 */,
  237UL /* 4, 6, 8 11101101, MMseq2_proteins_4 */,
  851UL /* 5, 6, 10 1101010011, MMseq2_proteins_5 */,
  981UL /* 6, 7, 10 1111010101, MMseq2_proteins_6 */,
  1715UL /* 7, 7, 11 11010110011, MMseq2_proteins_7 */,
  3699UL, /* 8, 8, 12 111001110011*/
  7399UL /* 9, 9, 13 1110011100111*/
};
#endif

#ifdef TESTSEED
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  29UL,
};
#endif

constexpr const uint8_t seed_table_size = sizeof(gt_spaced_seed_spec_tab)/sizeof(size_t);

template<class ScoreClass,const uint8_t seed_idx>
void process(GttlMultiseq* query,GttlMultiseq* target,const double sensitivity,
            const bool with_simd,const bool show,const bool short_header, const bool mmseqs, const bool correct, const double correct_ratio, const size_t num_threads)
{
  constexpr const size_t seed = gt_spaced_seed_spec_tab[seed_idx];
  const MMseqs2<ScoreClass,InvIntHashFunc,seed> mmseqs2{query,target,sensitivity,with_simd,short_header,show,mmseqs,correct,correct_ratio,num_threads};
}

int main(int argc, char *argv[])
{
  /* Different variables used for the optionparser as well as multiseq variable
     multiseq */
  SpacedSeedOptions options;

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }

  if(options.list_seed_option_is_set())
  {
    for(size_t seed_idx = 0; seed_idx < seed_table_size; seed_idx++)
    {
      std::cout << gt_spaced_seed_spec_tab[seed_idx] << std::endl;
    }
  }

  GttlMultiseq *query = nullptr, *target = nullptr;

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  if(inputfiles.size() == 0)
  {
    return EXIT_SUCCESS;
  }

  try
  {
    query = new GttlMultiseq(inputfiles[0],true,UINT8_MAX);
    target = new GttlMultiseq(inputfiles[1],true,UINT8_MAX);
  }
  catch (std::string &msg)
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    delete query;
    delete target;
    return EXIT_FAILURE;
  }

  if (options.short_header_option_is_set())
  {
    query->short_header_cache_create<'|','|'>();
    target->short_header_cache_create<'|','|'>();
  }

  const std::string &seed = options.seeds_get();
  const size_t seed_idx = std::stoi(seed);

  const std::string &correct_ratio_str = options.correct_ratio_get();
  const double correct_ratio = std::stod(correct_ratio_str);

  const std::string &num_threads_str = options.num_threads_get();
  const size_t num_threads = std::stoi(num_threads_str);

  const bool show = options.show_option_is_set();
  const bool with_simd = options.with_simd_is_set();
  const bool short_header = options.short_header_option_is_set();
  const double sensitivity = std::stod(options.sensitivity_get());
  const bool mmseqs = options.mmseqs_mode_is_set();
  const bool correct = options.local_score_correction_is_set();
  

  //const auto threshold_string = options.threshold_get();
  //const int8_t threshold = stoi(threshold_string);

  RunTimeClass rt{};
  constexpr_for<0,seed_table_size,1>([&] (auto seed_idx_constexpr)
  {
    if(seed_idx_constexpr == seed_idx)
    {
      process<Blosum62,seed_idx_constexpr>(query,target,sensitivity,with_simd,show,short_header,mmseqs,correct,correct_ratio,num_threads);
    }
  });
  
  if(!options.show_option_is_set())
  {
    rt.show("Generated environment");
    for (auto &&inputfile : inputfiles)
    {
      std::cout << "# filename\t" << inputfile << std::endl;
    }
    std::cout << "Query db statistics" << std::endl;
    for (auto &msg : query->statistics())
    {
      std::cout << "# " << msg << std::endl;
    }
    std::cout << "Target db statistics" << std::endl;
    for (auto &msg : target->statistics())
    {
      std::cout << "# " << msg << std::endl;
    }
  }
  delete query;
  delete target;
}