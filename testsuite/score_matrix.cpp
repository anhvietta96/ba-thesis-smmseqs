#include <iostream>
#include "utilities/cxxopts.hpp"
#include "filter/score_matrix.hpp"
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_nuc.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class ScoreMatrixOptions
{
 private:
  //std::vector<std::string> inputfiles{};
  bool help_option = false;
  std::string qgram_length = "1",
       matrix_option = "0";

 public:
  ScoreMatrixOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for score matrix");
    options.set_width(80);
    options.custom_help(std::string("[options]"));
    options.set_tab_expansion();
    options.add_options()
       ("m,matrix", "score matrix, 0 for nucleotide unit score, 1 for Blosum62",
        cxxopts::value<std::string>(matrix_option)->default_value("0"))
       ("q,qgram_length", "length of q-gram",
        cxxopts::value<std::string>(qgram_length)->default_value("1"))
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
  std::string get_matrix_option(void) const noexcept
  {
    return matrix_option;
  }
  std::string get_qgram_length(void) const noexcept
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
  ScoreMatrixOptions options;

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

  if(options.get_qgram_length() == "1")
  {
    if(options.get_matrix_option() == "0")
    {
      static constexpr const QgramScoreMatrix<Unit_score_nuc,1> sm{};
      constexpr const size_t sorted_size = sm.sorted_size_get();
      constexpr const size_t unsorted_size = sm.unsorted_size_get();

      for(size_t sorted_code = 0; sorted_code < sorted_size; sorted_code++)
      {
        for(size_t unsorted_code = 0; unsorted_code < unsorted_size; unsorted_code++)
        {
          std::cout << sm.score_get(sorted_code,unsorted_code) << '\t';
        }
        std::cout << std::endl;
      }
    }
    else if(options.get_matrix_option() == "1")
    {
      static constexpr const QgramScoreMatrix<Blosum62,1> sm{};
      constexpr const size_t sorted_size = sm.sorted_size_get();
      constexpr const size_t unsorted_size = sm.unsorted_size_get();

      for(size_t sorted_code = 0; sorted_code < sorted_size; sorted_code++)
      {
        for(size_t unsorted_code = 0; unsorted_code < unsorted_size; unsorted_code++)
        {
          std::cout << sm.score_get(sorted_code,unsorted_code) << '\t';
        }
        std::cout << std::endl;
      }
    }
    else
    {
      std::cerr << "Unaccounted score matrix" << std::endl;
    }
  }
  else if(options.get_qgram_length() == "2")
  {
    if(options.get_matrix_option() == "0")
    {
      static constexpr const QgramScoreMatrix<Unit_score_nuc,2> sm{};
      constexpr const size_t sorted_size = sm.sorted_size_get();
      constexpr const size_t unsorted_size = sm.unsorted_size_get();

      for(size_t sorted_code = 0; sorted_code < sorted_size; sorted_code++)
      {
        for(size_t unsorted_code = 0; unsorted_code < unsorted_size; unsorted_code++)
        {
          std::cout << sm.score_get(sorted_code,unsorted_code) << '\t';
        }
        std::cout << std::endl;
      }
    }/*
    else if(options.get_matrix_option() == "2")
    {
      static constexpr const QgramScoreMatrix<Blosum62,2> sm{};
      constexpr const size_t sorted_size = sm.sorted_size_get();
      constexpr const size_t unsorted_size = sm.unsorted_size_get();

      for(size_t sorted_code = 0; sorted_code < sorted_size; sorted_code++)
      {
        for(size_t unsorted_code = 0; unsorted_code < unsorted_size; unsorted_code++)
        {
          std::cout << sm.score_get(sorted_code,unsorted_code) << '\t';
        }
        std::cout << std::endl;
      }
    }*/
    else
    {
      std::cerr << "Unaccounted score matrix" << std::endl;
    }
  }
  else
  {
    std::cerr << "Unaccounted qgram length" << std::endl;
  }
  return EXIT_SUCCESS;
}
