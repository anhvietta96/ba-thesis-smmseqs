#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "filter/spaced_seeds.hpp"

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

constexpr const std::array<std::string_view,1> seed_list = {{ "10101010" }};

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
    for(size_t seed_idx = 0; seed_idx < seed_list.size(); seed_idx++)
    {
      std::cout << seed_list[seed_idx] << std::endl;
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

  auto total_seq_num = multiseq->sequences_number_get();
  static constexpr const char seed[] = "10101010";
  /*
  const std::string &seeds = options.seeds_get();
  if(seeds.size() > 1)
  {
    std::cout << "Multiple seeds not yet supported" << std::endl;
    return EXIT_SUCCESS;
  }
  if(seeds[0] > 57 or seeds[0] < 48)
  {
    std::cout << "Invalid seed" << std::endl;
    return EXIT_SUCCESS;
  }
  static constexpr const char* seed = seed_list[seeds[0]-48].data();
  */
  RunTimeClass rt{};
  if (options.protein_option_is_set())
  {
    static constexpr const char amino_acids[]
      = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
    constexpr const GttlSpacedSeed<amino_acids,20,seed> spaced_seed;
    constexpr const size_t seed_len = sizeof(seed)-1;
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
    {
      size_t code = 0;
      const char* curr_seq = multiseq->sequence_ptr_get(seqnum);
      const size_t seq_len = multiseq->sequence_length_get(seqnum);
      if(seq_len >= seed_len)
      {
        for(size_t i = 0; i < seq_len - seed_len + 1; i++)
        {
          code += spaced_seed.encode(curr_seq+i,seed_len);
        }
        if(options.show_option_is_set())
        {
          std::cout << code << std::endl;
        }
      }
    }
  }
  else
  {
    static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
    constexpr const GttlSpacedSeed<nucleotides_upper_lower_ACTG,4,seed> spaced_seed;
    constexpr const size_t seed_len = sizeof(seed)-1;
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
    {
      size_t code = 0;
      const char* curr_seq = multiseq->sequence_ptr_get(seqnum);
      const size_t seq_len = multiseq->sequence_length_get(seqnum);
      if(seq_len >= seed_len)
      {
        for(size_t i = 0; i < seq_len - seed_len + 1; i++)
        {
          code += spaced_seed.encode(curr_seq+i,seed_len);
        }
        if(options.show_option_is_set())
        {
          std::cout << code << std::endl;
        }
      }
    }
  }
  if(!options.show_option_is_set())
  {
    rt.show("Encoding runtime");
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
