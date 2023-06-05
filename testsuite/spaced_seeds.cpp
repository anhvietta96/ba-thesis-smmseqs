#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "filter/spaced_seeds.hpp"
#include "utilities/constexpr_for.hpp"


/*TODO
finish spaced seed test (seed choice) ? template x
pull x
score matrix x environment
check err ? new compile x
1D map sorted qmer x
*/

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions
{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false,
       protein_option = false,
       list_seeds = false,
       show = false;
  std::string seeds;

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for sum of spaced seed intcode");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("p,protein", "handle protein sequences",
        cxxopts::value<bool>(protein_option)->default_value("false"))
       ("l,list_seeds", "output available seed list",
        cxxopts::value<bool>(list_seeds)->default_value("false"))
       ("d,seeds", "choose seeds to compute",
        cxxopts::value<std::string>(seeds)->default_value("0"))
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
  bool protein_option_is_set(void) const noexcept
  {
    return protein_option;
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
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};


static constexpr const size_t gt_spaced_seed_spec_tab[] = {
  29UL /* 0, 4, 5 11101, MMseq2_proteins_0 */,
  59UL /* 1, 5, 6 111011, MMseq2_proteins_1 */,
  107UL /* 2, 5, 7 1101011, MMseq2_proteins_2 */,
  3205UL /* 3, 5, 12 110010000101, MMseq2_proteins_3 */,
  237UL /* 4, 6, 8 11101101, MMseq2_proteins_4 */,
  851UL /* 5, 6, 10 1101010011, MMseq2_proteins_5 */,
  981UL /* 6, 7, 10 1111010101, MMseq2_proteins_6 */,
  1715UL /* 7, 7, 11 11010110011, MMseq2_proteins_7 */
};

constexpr const uint8_t seed_table_size = sizeof(gt_spaced_seed_spec_tab)/sizeof(size_t);

template<const char* char_spec, const size_t undefined_rank,const uint8_t seed_idx>
constexpr void process(GttlMultiseq* multiseq, const bool show)
{
  constexpr const size_t seed = gt_spaced_seed_spec_tab[seed_idx];
  constexpr const GttlSpacedSeed<char_spec,undefined_rank,seed> spaced_seeds{};
  
  const auto total_seq_num = multiseq->sequences_number_get();
  for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
  {
    size_t code = 0;
    const char* curr_seq = multiseq->sequence_ptr_get(seqnum);
    const size_t seq_len = multiseq->sequence_length_get(seqnum);
    const size_t seed_len = spaced_seeds.span_get();
    if(seq_len >= seed_len)
    {
      for(size_t i = 0; i < seq_len - seed_len + 1; i++)
      {
        code += spaced_seeds.encode(curr_seq+i);
      }
      if(show)
      {
      std::cout << code << std::endl;
      }
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

  GttlMultiseq *multiseq = nullptr;

  if(options.list_seed_option_is_set())
  {
    for(size_t seed_idx = 0; seed_idx < seed_table_size; seed_idx++)
    {
      std::cout << gt_spaced_seed_spec_tab[seed_idx] << std::endl;
    }
  }

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  if(inputfiles.size() == 0)
  {
    return EXIT_SUCCESS;
  }

  try
  {
    multiseq = new GttlMultiseq(inputfiles,true,UINT8_MAX);
  }
  catch (std::string &msg)
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    delete multiseq;
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

  RunTimeClass rt{};
  constexpr_for<0,seed_table_size,1>([&] (auto compile_time_idx)
  {
    if(seed_idx==compile_time_idx and options.protein_option_is_set())
    {
      static constexpr const char amino_acids[]
      = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
      process<amino_acids,20,compile_time_idx>(multiseq,show);
    }
    else if(seed_idx==compile_time_idx)
    {
      static constexpr const char nucleotides_upper_lower_ACTG[]
      = "Aa|Cc|TtUu|Gg";
      process<nucleotides_upper_lower_ACTG,4,compile_time_idx>(multiseq,show);
    }
  });
  
  if(!options.show_option_is_set())
  {
    rt.show("Encoding of spaced seeds");
    for (auto &&inputfile : inputfiles)
    {
      std::cout << "# filename\t" << inputfile << std::endl;
    }
    for (auto &msg : multiseq->statistics())
    {
      std::cout << "# " << msg << std::endl;
    }
  }
  delete multiseq;
}