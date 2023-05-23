#include <iostream>
#include "utilities/cxxopts.hpp"
#include "filter/sorted_q_mer.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

constexpr const std::array<std::array<char,3>,20> nuc_3 = {{
    std::array<char,3>({{ 'A','A','A' }}),
    std::array<char,3>({{ 'A','A','C' }}),
    std::array<char,3>({{ 'A','A','T' }}),
    std::array<char,3>({{ 'A','A','G' }}),
    std::array<char,3>({{ 'A','C','C' }}),    
    std::array<char,3>({{ 'A','C','T' }}),    
    std::array<char,3>({{ 'A','C','G' }}),    
    std::array<char,3>({{ 'A','T','T' }}),    
    std::array<char,3>({{ 'A','T','G' }}),    
    std::array<char,3>({{ 'A','G','G' }}),    
    std::array<char,3>({{ 'C','C','C' }}),    
    std::array<char,3>({{ 'C','C','T' }}),    
    std::array<char,3>({{ 'C','C','G' }}),    
    std::array<char,3>({{ 'C','T','T' }}),    
    std::array<char,3>({{ 'C','T','G' }}),    
    std::array<char,3>({{ 'C','G','G' }}),    
    std::array<char,3>({{ 'T','T','T' }}),    
    std::array<char,3>({{ 'T','T','G' }}),    
    std::array<char,3>({{ 'T','G','G' }}),    
    std::array<char,3>({{ 'G','G','G' }}),    
}};

constexpr const std::array<std::array<char,4>,35> nuc_4 = {{
    std::array<char,4>({{ 'A','A','A','A' }}),
    std::array<char,4>({{ 'A','A','A','C' }}),
    std::array<char,4>({{ 'A','A','A','T' }}),
    std::array<char,4>({{ 'A','A','A','G' }}),
    std::array<char,4>({{ 'A','A','C','C' }}),    
    std::array<char,4>({{ 'A','A','C','T' }}),    
    std::array<char,4>({{ 'A','A','C','G' }}),    
    std::array<char,4>({{ 'A','A','T','T' }}),    
    std::array<char,4>({{ 'A','A','T','G' }}),    
    std::array<char,4>({{ 'A','A','G','G' }}),    
    std::array<char,4>({{ 'A','C','C','C' }}),    
    std::array<char,4>({{ 'A','C','C','T' }}),    
    std::array<char,4>({{ 'A','C','C','G' }}),    
    std::array<char,4>({{ 'A','C','T','T' }}),    
    std::array<char,4>({{ 'A','C','T','G' }}),    
    std::array<char,4>({{ 'A','C','G','G' }}),    
    std::array<char,4>({{ 'A','T','T','T' }}),    
    std::array<char,4>({{ 'A','T','T','G' }}),    
    std::array<char,4>({{ 'A','T','G','G' }}),    
    std::array<char,4>({{ 'A','G','G','G' }}),    
    std::array<char,4>({{ 'C','C','C','C' }}),    
    std::array<char,4>({{ 'C','C','C','T' }}),    
    std::array<char,4>({{ 'C','C','C','G' }}),    
    std::array<char,4>({{ 'C','C','T','T' }}),    
    std::array<char,4>({{ 'C','C','T','G' }}),    
    std::array<char,4>({{ 'C','C','G','G' }}),    
    std::array<char,4>({{ 'C','T','T','T' }}),    
    std::array<char,4>({{ 'C','T','T','G' }}),    
    std::array<char,4>({{ 'C','T','G','G' }}),    
    std::array<char,4>({{ 'C','G','G','G' }}),    
    std::array<char,4>({{ 'T','T','T','T' }}),    
    std::array<char,4>({{ 'T','T','T','G' }}),    
    std::array<char,4>({{ 'T','T','G','G' }}),    
    std::array<char,4>({{ 'T','G','G','G' }}),    
    std::array<char,4>({{ 'G','G','G','G' }}),
    }};

class SortedQmerOptions
{
 private:
  //std::vector<std::string> inputfiles{};
  bool help_option = false,
       protein_option = false,
       display_option = false;
  std::string qgram_length = "3";

 public:
  SortedQmerOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for indexing sorted qmer");
    options.set_width(80);
    options.custom_help(std::string("[options]"));
    options.set_tab_expansion();
    options.add_options()
       ("p,protein", "handle protein sequences",
        cxxopts::value<bool>(protein_option)->default_value("false"))
       ("q,qgram_length", "length of q-gram",
        cxxopts::value<std::string>(qgram_length)->default_value("3"))
       ("s,show", "output computed sorted qmer",
        cxxopts::value<bool>(display_option)->default_value("false"))
       ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      } /*else
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
      }*/
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
  bool display_option_is_set(void) const noexcept
  {
    return display_option;
  }
  const std::string get_qgram_length(void) const noexcept
  {
    return qgram_length;
  }
  /*const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }*/
};


int main(int argc, char *argv[])
{
  static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
  static constexpr const char amino_acids[]
      = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
  constexpr const SortedQmer<nucleotides_upper_lower_ACTG,4,2> map_nuc_2{};
  constexpr const SortedQmer<nucleotides_upper_lower_ACTG,4,3> map_nuc_3{};
  constexpr const SortedQmer<nucleotides_upper_lower_ACTG,4,4> map_nuc_4{};
  constexpr const SortedQmer<amino_acids,20,2> map_prot_2{};
  constexpr const SortedQmer<amino_acids,20,3> map_prot_3{};

  //static_assert(map.get_seq_num(0)[0] != 0);

  SortedQmerOptions options;

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

  if (options.protein_option_is_set())
  {
    if(options.get_qgram_length() == "2")
    {
      const size_t map_size = map_prot_2.get_size();
      for(size_t i = 0; i < map_size; i++)
      {
        const std::array<char,2> qgram = map_prot_2.get_seq_num(i);
        if(options.display_option_is_set())
        {
          const std::string output_qgram(std::begin(qgram),std::end(qgram));
          std::cout << output_qgram << std::endl;
        }
      }
    }
    else if(options.get_qgram_length() == "3")
    {
      const size_t map_size = map_prot_3.get_size();
      for(size_t i = 0; i < map_size; i++)
      {
        const std::array<char,3> qgram = map_prot_3.get_seq_num(i);
        if(options.display_option_is_set())
        {
          const std::string output_qgram(std::begin(qgram),std::end(qgram));
          std::cout << output_qgram << std::endl;
        }
      }
    }
    else
    {
      std::cout << "Unaccounted qgram length" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  else
  {
    if(options.get_qgram_length() == "3")
    {
      const size_t map_size = map_nuc_3.get_size();
      for(size_t i = 0; i < map_size; i++)
      {
        const std::array<char,3> qgram = map_nuc_3.get_seq_num(i);
        if(options.display_option_is_set())
        {
          const std::string output_qgram(std::begin(qgram),std::end(qgram));
          std::cout << output_qgram << std::endl;
        }
      }
    }
    else if(options.get_qgram_length() == "4")
    {
      const size_t map_size = map_nuc_4.get_size();
      for(size_t i = 0; i < map_size; i++)
      {
        const std::array<char,4> qgram = map_nuc_4.get_seq_num(i);
        if(options.display_option_is_set())
        {
          const std::string output_qgram(std::begin(qgram),std::end(qgram));
          std::cout << output_qgram << std::endl;
        }
      }
    }
    else
    {
      std::cout << "Unaccounted qgram length" << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}