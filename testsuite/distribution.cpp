#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utilities/cxxopts.hpp"
#include "sequences/literate_multiseq.hpp"
#include "filter/distribution.hpp"
#include "alignment/blosum62.hpp"
#include "utilities/constexpr_for.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false;
  std::string sensitivity = "1";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for sum of spaced seed intcode");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("s,sensitivity", "set sensitivity",
        cxxopts::value<std::string>(sensitivity)->default_value("1"))
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
  const std::string& sensitivity_get(void) const noexcept
  {
    return sensitivity;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};

int main(int argc, char *argv[]){
  /*constexpr const uint8_t max_num_convolution = 6;
  constexpr const Distribution<Blosum62,max_num_convolution> distribution{};
  for(uint8_t i = 2; i <= max_num_convolution; i++){
    std::cout << (int)(i+1) << '\t' << (int) distribution.threshold_get(i) << std::endl;
  }*/

  SpacedSeedOptions options;
  options.parse(argc,argv);

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  const std::string& sensitivity_str = options.sensitivity_get();
  const double sensitivity = std::stod(sensitivity_str);
  std::cout << "Read sensitivity: " << sensitivity << std::endl;
  GttlMultiseq* multiseq = new GttlMultiseq(inputfiles,true,UINT8_MAX);
  constexpr const Blosum62 sc{};
  constexpr const auto char_spec = sc.character_spec;
  constexpr const auto undefined_rank = sc.num_of_chars;
  const LiterateMultiseq<char_spec,undefined_rank> lit_multiseq{*multiseq};

  constexpr_for<1,8,1>([&] (auto weight){
    const BGDistribution<Blosum62,weight> distribution{};
    const auto custom_threshold_wo_dis = distribution.custom_threshold_get(sensitivity);
    const auto custom_threshold_with_dis = distribution.custom_threshold_get(lit_multiseq.rank_dist_get(),sensitivity);
    std::cout << (int) weight << '\t' << (int) custom_threshold_wo_dis << '\t' << custom_threshold_with_dis << std::endl;
  });

  delete multiseq;

  //std::cout << (int) distribution.threshold(0.95) << std::endl;
}