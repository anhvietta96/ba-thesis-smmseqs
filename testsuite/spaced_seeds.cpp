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

class MultiseqOptions
{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false,
       protein_option = false,
       zipped_option = false,
       rankdist_option = false,
       short_header_option = false;
  int width_arg = -1;

 public:
  MultiseqOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for GttlMultiseq");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("p,protein", "handle protein sequences",
        cxxopts::value<bool>(protein_option)->default_value("false"))
       ("z,zipped", "expect two fastq  files with the same "
                    "number of sequences; show them "
                    "in zipped order, i.e. the "
                    "sequences at even indexes (when "
                    "counting from 0) are from the first "
                    "file and sequences at odd indexes are "
                    "from the second file",
        cxxopts::value<bool>(zipped_option)->default_value("false"))
       ("r,rankdist", "output distribution of ranks of "
                      "transformed sequences",
        cxxopts::value<bool>(rankdist_option)->default_value("false"))
       ("s,short_header", "show header up to and excluding the first blank",
        cxxopts::value<bool>(short_header_option)->default_value("false"))
       ("w,width", "output headers and sequences; "
                   "width specifies the linewidth of the"
                   "sequence output; 0 means to output\n"
                   "a sequence in a single line",
        cxxopts::value<int>(width_arg)->default_value("-1"))
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
        if (unmatched_args.size() == 0)
        {
          throw std::invalid_argument("at least one inputput file is required");
        }
        for (size_t idx = 0; idx < unmatched_args.size(); idx++)
        {
          inputfiles.push_back(unmatched_args[idx]);
        }
      }
      if (zipped_option && inputfiles.size() != 2)
      {
        throw std::invalid_argument("option -z/--zipped requires exactly "
                                    "two files");
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
  bool zipped_option_is_set(void) const noexcept
  {
    return zipped_option;
  }
  bool rankdist_option_is_set(void) const noexcept
  {
    return rankdist_option;
  }
  bool short_header_option_is_set(void) const noexcept
  {
    return short_header_option;
  }
  int width_option_get(void) const noexcept
  {
    return width_arg;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};


int main(int argc, char *argv[])
{
  /* Different variables used for the optionparser as well as multiseq variable
     multiseq */
  MultiseqOptions options;

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
  

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
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
  if (options.short_header_option_is_set())
  {
    multiseq->short_header_cache_create();
  }
  

  auto total_seq_num = multiseq->sequences_number_get();
  size_t code = 0;
  size_t seq_len;
  const char* curr_seq;
  static constexpr const char seed[] = "10101010";
  static constexpr const size_t seed_len = sizeof(seed)-1;

  RunTimeClass rt{};
  if (options.protein_option_is_set())
  {
    static constexpr const char amino_acids[]
      = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
    GttlSpacedSeed<amino_acids,20,seed> spaced_seed;
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
    {
      code = 0;
      curr_seq = multiseq->sequence_ptr_get(seqnum);
      seq_len = multiseq->sequence_length_get(seqnum);
      if(seq_len >= seed_len)
      {
        for(size_t i = 0; i < seq_len - seed_len + 1; i++)
        {
            code += spaced_seed.encode(curr_seq+i,seed_len);
            //std::cout << curr_seq << '\t' << code << std::endl;
        }
        std::cout << code << std::endl;
      }
    }
  }
  else
  {
    static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
    GttlSpacedSeed<nucleotides_upper_lower_ACTG,4,seed> spaced_seed;
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
    {
      code = 0;
      curr_seq = multiseq->sequence_ptr_get(seqnum);
      seq_len = multiseq->sequence_length_get(seqnum);
      if(seq_len >= seed_len)
      {
        for(size_t i = 0; i < seq_len - seed_len + 1; i++)
        {
            code += spaced_seed.encode(curr_seq+i,seed_len);
            //std::cout << curr_seq << '\t' << code << std::endl;
        }
        std::cout << code << std::endl;
      }
    }
  }
  rt.show("Runtime");
  
  for (auto &&inputfile : inputfiles)
  {
    std::cout << "# filename\t" << inputfile << std::endl;
  }
  for (auto &msg : multiseq->statistics())
  {
    std::cout << "# " << msg << std::endl;
  }
  delete multiseq;
}
