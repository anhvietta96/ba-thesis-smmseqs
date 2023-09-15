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
#define TESTSEED

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
       show = false;
  std::string seeds = "0", threshold = "20";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"run MMseqs2 on a query & a target database");
    options.set_width(80);
    options.custom_help(std::string("[options] query_db target_db"));
    options.set_tab_expansion();
    options.add_options()
       ("l,list_seeds", "output available seed list",
        cxxopts::value<bool>(list_seeds)->default_value("false"))
       ("d,seeds", "choose seeds to compute",
        cxxopts::value<std::string>(seeds)->default_value("0"))
       ("w,with_simd", "compute with SIMD",
        cxxopts::value<bool>(with_simd)->default_value("false"))
       ("t,threshold", "set threshold",
        cxxopts::value<std::string>(threshold)->default_value("20")) 
       ("s,show", "output sum of intcode",
        cxxopts::value<bool>(show)->default_value("false"))
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
  const std::string &threshold_get(void) const noexcept
  {
    return threshold;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
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
  127UL,
  981UL /* 6, 7, 10 1111010101, MMseq2_proteins_6 */,
  1715UL /* 7, 7, 11 11010110011, MMseq2_proteins_7 */
};
#endif

#ifdef TESTSEED
static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  15UL,
};
#endif
constexpr const uint8_t seed_table_size = sizeof(gt_spaced_seed_spec_tab)/sizeof(size_t);

template<const char* char_spec, const size_t undefined_rank,const uint8_t seed_idx>
void process(GttlMultiseq* query,GttlMultiseq* target,const bool show, const bool with_simd)
{
  constexpr const size_t seed = gt_spaced_seed_spec_tab[seed_idx];
  MMseqs2<Blosum62,InvIntHashFunc,seed> mmseqs2{query,target,with_simd};
  
  if(show){
    std::cout << "#query_seq_num" << '\t' << "query_seq_pos" << '\t' << 
      "target_seq_num" << '\t' << "target_seq_pos" << '\t' << 
      "score" << '\t' << "code" << '\t' << std::endl;
    for(size_t i = 0; i < mmseqs2.size(); i++){
      const auto hit = mmseqs2.hit_get(i);
      std::cout << (int) hit.query_seq_num << '\t' << (int) hit.query_seq_pos << '\t' << 
      (int) hit.target_seq_num << '\t' << (int) hit.target_seq_pos << '\t' << 
      (int) hit.score << '\t' << (int) hit.code << '\t' << std::endl;
    }
  }
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

  const std::string &seeds = options.seeds_get();
  if(seeds.size() > 1)
  {
    std::cout << "Multiple seeds not supported" << std::endl;
    return EXIT_SUCCESS;
  }
  if(seeds[0] >= static_cast<uint8_t>('0' + seed_table_size) or seeds[0] < static_cast<uint8_t>('0'))
  {
    std::cout << "Invalid seed" << std::endl;
    return EXIT_SUCCESS;
  }
  const uint8_t seed_idx = seeds[0] - '0';

  const bool show = options.show_option_is_set();
  const bool with_simd = options.with_simd_is_set();

  //const auto threshold_string = options.threshold_get();
  //const int8_t threshold = stoi(threshold_string);

  RunTimeClass rt{};
  constexpr_for<0,seed_table_size,1>([&] (auto seed_idx_constexpr)
  {
    if(seed_idx_constexpr == seed_idx)
    {
      static constexpr const Blosum62 sc{};
      static constexpr const auto char_spec = sc.character_spec;
      static constexpr const auto undefined_rank = sc.num_of_chars;
      process<char_spec,undefined_rank,seed_idx_constexpr>(query,target,show,with_simd);
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