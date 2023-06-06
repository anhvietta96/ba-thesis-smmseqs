#include "filter/qgram_environment.hpp"
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "sequences/alphabet.hpp"
#include <iostream>

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class EnvOptions
{
 private:
  //std::vector<std::string> inputfiles{};
  bool help_option = false;
  std::string qgram_length = "1",
       matrix_option = "0";

 public:
  EnvOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for score matrix");
    options.set_width(80);
    options.custom_help(std::string("[options]"));
    options.set_tab_expansion();
    options.add_options()
       ("m,matrix", "score matrix, 0 for Blosum62, 1 for nucleotide unit score func",
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

int main(int argc, char argv[])
{
  constexpr const Blosum62 b62{};
  constexpr const auto char_spec = b62.characters;
  constexpr const auto undefined_rank = b62.num_of_chars;
  constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  
  constexpr const size_t qgram_length = 2;
  constexpr const size_t threshold = 10;

  constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  const QgramEnvironment<Blosum62,qgram_length,threshold> env{};
  const size_t sorted_size = env.sorted_size_get();
  for(size_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
  {
    const auto refqgram = sorted_q.extern_qgram_get(sorted_idx);
    const std::string output_refqgram(std::begin(refqgram),std::end(refqgram));
    std::cout << output_refqgram << ":\t";
    for(size_t i = 0; i < env.size_get(sorted_idx); i++)
    {
      int8_t score = env.env_get(sorted_idx,i).first;
      std::cout << (int) score << '\t';
      /*
      auto qgram = env.env_get(sorted_idx,i).second;
      for(size_t j = 0; j < qgram_length; j++)
      {
        qgram[j] += 65;
      }
      const std::string output_qgram(std::begin(qgram),std::end(qgram));
      std::cout << output_qgram << '\t';*/
    }
    std::cout << '\n';
  }
}