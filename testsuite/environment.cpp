#include "filter/qgram_environment.hpp"
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "sequences/alphabet.hpp"
#include <iostream>
/*
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
  std::string get_matrix_option(void) const noexcept
  {
    return matrix_option;
  }
  std::string get_qgram_length(void) const noexcept
  {
    return qgram_length;
  }
};
*/
/*
template<const uint16_t arr_len,const uint16_t sorted_idx>
std::array<uint16_t,arr_len> threshold_idx_arr(std::vector<std::pair<int8_t,uint16_t>> sorted_env) const
{
  const auto max_score = sorted_env[0].first;
  const uint16_t env_size = sorted_env.size()

  std::array<uint16_t,arr_len> threshold_idx_arr{};
  int8_t curr_score = max_score;
  uint16_t threshold_idx_arr_idx = static_cast<uint16_t>(arr_len-1);

  for(uint16_t env_idx = 0; env_idx < sorted_env.size(); env_idx++)
  {
    if(sorted_env[env_idx].first != curr_score)
    {
      const uint16_t diff = curr_score - sorted_env[env_idx].first;
      for(uint16_t i = 0; i < diff; i++)
      {
        threshold_idx_arr[threshold_idx_arr_idx] = env_idx;
        threshold_idx_arr_idx--;
      }
      curr_score = sorted_env[env_idx].first;
    }
  }
}*/

int main(/*int argc, char argv[]*/)
{
  constexpr const Blosum62 b62{};
  constexpr const auto char_spec = b62.character_spec;
  constexpr const auto undefined_rank = b62.num_of_chars;
  constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  
  constexpr const size_t qgram_length = 2;
  constexpr const size_t threshold = 1;

  //constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  //constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  /*for(size_t unsorted_idx = 0; unsorted_idx < unsorted_q.size_get(); unsorted_idx++)
  {
    auto qgram = unsorted_q.extern_qgram_get(unsorted_idx);
    const std::string output_qgram(std::begin(qgram),std::end(qgram));
    std::cout << output_qgram << std::endl;
  }*/
  const QgramEnvironment<Blosum62,qgram_length,threshold> env{};
  const size_t sorted_size = env.sorted_size_get();
  for(size_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
  {/*
    const auto refqgram = sorted_q.extern_qgram_get(sorted_idx);
    const std::string output_refqgram(std::begin(refqgram),std::end(refqgram));
    std::cout << output_refqgram << ":\t";*/
    std::cout << (int) sorted_idx << ":\t";
    for(size_t i = 0; i < env.size_get(sorted_idx); i++)
    {
      const auto score = env.env_get(sorted_idx,i).score;
      const auto unsorted_qgram_code = env.env_get(sorted_idx,i).code;
      std::cout << '(' << (int) score << ',' << (int) unsorted_qgram_code << ')' << '\t';
     /*
      const auto unsorted_qgram_code = env.env_get(sorted_idx,i).second;
      std::cout << unsorted_qgram_code << '\t';
      const auto qgram = unsorted_q.extern_qgram_get(unsorted_qgram_code);
      const std::string output_qgram(std::begin(qgram),std::end(qgram));
      std::cout << output_qgram << '\t';*/
    }
    std::cout << '\n';
    /*
    constexpr const uint16_t sq_code = 1;
    const auto qgram_env = env.qgram_env(sq_code);
    const int8_t max_score = qgram_env[0].first;
    const int8_t min_score = qgram_env[qgram_env.size()-1].first;
    const uint16_t score_diff = static_cast<uint16_t>(max_score-min_score);
    const auto threshold_arr_idx = threshold_arr_idx<score_diff,>()*/
  } 
}