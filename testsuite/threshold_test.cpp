#include <iostream>
#include "alignment/blosum62.hpp"
#include "filter/utils.hpp"
#include "sequences/literate_multiseq.hpp"
#include <array>
#include "utilities/cxxopts.hpp"
#include "utilities/constexpr_for.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false;
  std::string qgram_length = "3";
  std::string threshold = "18";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for sum of spaced seed intcode");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
        ("q,qgram_length", "set qgram length",
        cxxopts::value<std::string>(qgram_length)->default_value("3"))
        ("t,threshold", "set threshold",
        cxxopts::value<std::string>(threshold)->default_value("18"))
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
  const std::string& threshold_get(void) const noexcept
  {
    return threshold;
  }
  const std::string& qgram_length_get(void) const noexcept
  {
    return qgram_length;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};


template<const uint8_t template_qgram_length, const uint64_t undefined_rank>
inline void reset_qgram(std::array<uint8_t,template_qgram_length>& qgram){
  for(uint8_t idx = template_qgram_length - 1; idx != 0; idx--){
    if(qgram[idx] == undefined_rank){
      qgram[idx] = 0;
      qgram[idx-1]++;
    }
    else break;
  }
}

template<const uint8_t template_qgram_length, class ScoreClass>
inline int8_t get_score(const std::array<uint8_t,template_qgram_length>& qgram_1,
                const std::array<uint8_t,template_qgram_length>& qgram_2){
  static constexpr const ScoreClass sc{};
  int8_t score = 0;
  for(uint8_t i = 0; i < template_qgram_length; i++){
    for(uint8_t j = 0; j < template_qgram_length; j++){
      score += sc.score_matrix[qgram_1[i]][qgram_2[j]];
    }
  }
  return score;
}

int main(int argc, char *argv[]){
  SpacedSeedOptions options;
  options.parse(argc,argv);

  static constexpr const Blosum62 sc{};
  static constexpr const uint64_t undefined_rank = sc.num_of_chars;
  uint64_t count = 0;
  
  const uint8_t input_qgram_length = std::stoi(options.qgram_length_get());
  const int8_t threshold = std::stoi(options.threshold_get());
  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  GttlMultiseq* multiseq = new GttlMultiseq(inputfiles,true,UINT8_MAX);
  const LiterateMultiseq<sc.character_spec,sc.num_of_chars> literate_multiseq{multiseq};
  const auto& target_distribution = literate_multiseq.rank_dist_get();
  
  size_t sum_count = 0;
  for(uint8_t i = 0; i < undefined_rank; i++){
    sum_count += target_distribution[i];
  }
  std::array<double,undefined_rank> freq{};
  for(uint8_t i = 0; i < undefined_rank; i++){
    freq[i] = static_cast<double>(target_distribution[i]) / sum_count;
  }

  constexpr_for<1,5,1>([&] (auto qgram_length){
    if(qgram_length == input_qgram_length){
      static constexpr const uint64_t hist_width = (sc.highest_score - sc.smallest_score + 1)*qgram_length;
      static constexpr const int8_t lowest_score = sc.smallest_score * qgram_length;
      static constexpr const uint64_t size = constexpr_pow(undefined_rank,qgram_length);

      std::array<double,hist_width> hist{};
      std::array<uint8_t,qgram_length> qgram_1{};
      std::array<uint8_t,qgram_length> qgram_2{};
      
      for(size_t i = 0; i < size; i++){
        for(size_t j = 0; j < size; j++){
          int8_t score = 0;
          double sfreq = 1;
          for(uint8_t k = 0; k < qgram_length; k++){
            score += sc.score_matrix[qgram_1[k]][qgram_2[k]];
            sfreq *= freq[qgram_1[k]] * freq[qgram_2[k]];
          }
          count += score >= threshold;
          hist[score-lowest_score] += sfreq;

          qgram_2[qgram_length-1]++;
          if(qgram_2[qgram_length-1]==undefined_rank) reset_qgram<qgram_length,undefined_rank>(qgram_2);
        }
        qgram_1[qgram_length-1]++;
        if(qgram_1[qgram_length-1]==undefined_rank) reset_qgram<qgram_length,undefined_rank>(qgram_1);
        qgram_2[0] = 0;
      }
      std::cout << "Distribution:" << std::endl;
      for(uint64_t i = 0; i < hist_width; i++){
        std::cout << hist[i] << std::endl;
      }
      std::cout << "Expected sensitivity: " << std::endl;
      std::cout << static_cast<double>(count) / constexpr_pow(undefined_rank,qgram_length) << std::endl;
    }
  });

  delete multiseq;
}