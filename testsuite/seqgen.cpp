#include <iostream>
#include <ctime>
#include <cstdlib>
#include "sequences/literate_multiseq.hpp"
#include "alignment/blosum62.hpp"
#include "utilities/cxxopts.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false;
  std::string seqlen = "3000000";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for sum of spaced seed intcode");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("l,seqlen", "set sequence length",
        cxxopts::value<std::string>(seqlen)->default_value("3000000"))
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
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  const std::string& seqlen_get() const noexcept
  {
    return seqlen;
  }
};

int main(int argc, char *argv[]){
  SpacedSeedOptions options;
  options.parse(argc,argv);

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  const std::string& seqlen_str = options.seqlen_get();
  const size_t seqlen = std::stoi(seqlen_str);
  GttlMultiseq* multiseq = new GttlMultiseq(inputfiles,true,UINT8_MAX);
  constexpr const Blosum62 sc{};
  constexpr const auto char_spec = sc.character_spec;
  constexpr const auto undefined_rank = sc.num_of_chars;
  const LiterateMultiseq<char_spec,undefined_rank> lit_multiseq{*multiseq};
  const auto& target_distribution = lit_multiseq.rank_dist_get();
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  size_t sum_count = 0;
  for(uint8_t i = 0; i < undefined_rank; i++){
    sum_count += target_distribution[i];
  }
  std::array<double,undefined_rank> freq{};
  for(uint8_t i = 0; i < undefined_rank; i++){
    freq[i] = static_cast<double>(target_distribution[i]) / sum_count;
  }
  for(uint8_t i = 1; i < undefined_rank; i++){
    freq[i] += freq[i-1];
  }
  freq[undefined_rank-1] = 1;


  uint8_t* seq = new uint8_t[seqlen]{};
  std::srand(time(0));
  for(size_t i = 0; i < seqlen; i++){
    const double random = static_cast<double>(rand())/RAND_MAX;
    if(random < freq[0]){
      seq[i] = alpha.rank_to_char(0);
      continue;
    }
    uint8_t j = 1;
    for(; j < undefined_rank-1 and (freq[j-1] >= random or freq[j] < random); j++){};
    seq[i] = alpha.rank_to_char(j);
  }

  std::cout << std::string(seq,seq+seqlen-1) << std::endl;
  delete seq;
  delete multiseq;
}