#include "filter/InvIntHash.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "alignment/blosum62.hpp"
#include <iostream>
#include <algorithm>

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SpacedSeedOptions{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false,
       correctness = false;
  std::string seeds = "0";

 public:
  SpacedSeedOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"run test on seqdb: recursive hashing vs normal hashing");
    options.set_width(80);
    options.custom_help(std::string("[options] query_db target_db"));
    options.set_tab_expansion();
    options.add_options()
        ("c,correctness", "test for correctness",
        cxxopts::value<bool>(correctness)->default_value("false"))
        ("d,seeds", "choose seeds to compute",
        cxxopts::value<std::string>(seeds)->default_value("0"));
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
        if(unmatched_args.size() != 1){
          throw cxxopts::OptionException("Exactly 1 db are needed.");
        }
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
  bool test_for_correctness(void) const noexcept
  {
    return correctness;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  const std::string &seeds_get(void) const noexcept
  {
    return seeds;
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
  1715UL /* 7, 7, 11 11010110011, MMseq2_proteins_7 */,
  31709UL,
  56173UL
};

constexpr const uint8_t seed_table_size = sizeof(gt_spaced_seed_spec_tab)/sizeof(size_t);

size_t sizeof_unit_get(const size_t total_bits){
  if(total_bits <= 64) return 8;
  if(total_bits % 8 == 0) return total_bits/8;
  return total_bits/8+1;
}

template<class ScoreClass,const uint8_t seed_idx>
void process(GttlMultiseq* multiseq,const bool correctness){
  constexpr const size_t seed = gt_spaced_seed_spec_tab[seed_idx];
  static constexpr const Blosum62 sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;

  LiterateMultiseq<char_spec,undefined_rank> literate_multiseq{multiseq};
  literate_multiseq.perform_sequence_encoding();

  const Multiseq_Hash<ScoreClass,InvIntHashFunc,seed> multiseq_hash{};

  const size_t seq_len_bits = multiseq->sequences_length_bits_get();
  const size_t seq_num_bits = multiseq->sequences_number_bits_get();
  const size_t hashbits = multiseq_hash.hashbits_get();

  const auto sizeof_unit = sizeof_unit_get(hashbits+seq_num_bits+seq_len_bits);
  if(sizeof_unit>9){
    std::cerr << "sizeof_unit too large for computation" << std::endl;
  }

  constexpr_for<8,10,1>([&] (auto sizeof_unit_constexpr){
    if(sizeof_unit_constexpr == sizeof_unit){
      std::vector<BytesUnit<sizeof_unit_constexpr,3>> rechash_data;
      std::vector<BytesUnit<sizeof_unit_constexpr,3>> nhash_data;
      const GttlBitPacker<sizeof_unit_constexpr,3> packer{{hashbits,seq_num_bits,seq_len_bits}};
      
      RunTimeClass rt{};
      multiseq_hash.template linear_hash<sizeof_unit_constexpr>(multiseq,nhash_data,packer);
      rt.show("Normal hashing");
      multiseq_hash.template hash<sizeof_unit_constexpr>(multiseq,rechash_data,packer);  
      rt.show("Recursive hashing");

      if(correctness){
        assert(rechash_data.size() == nhash_data.size());
        for(size_t i = 0; i < rechash_data.size(); i++){
          assert(rechash_data[i] == nhash_data[i]);
        }
      }
    }
  });
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

  /*if(options.list_seed_option_is_set())
  {
    for(size_t seed_idx = 0; seed_idx < seed_table_size; seed_idx++)
    {
      std::cout << gt_spaced_seed_spec_tab[seed_idx] << std::endl;
    }
  }*/

  GttlMultiseq *multiseq = nullptr;

  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  if(inputfiles.size() == 0)
  {
    return EXIT_SUCCESS;
  }

  try
  {
    multiseq = new GttlMultiseq(inputfiles[0],true,UINT8_MAX);
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
  /*if(seeds.size() > 1)
  {
    std::cout << "Multiple seeds not supported" << std::endl;
    return EXIT_SUCCESS;
  }
  if(seeds[0] >= static_cast<uint8_t>('0' + seed_table_size) or seeds[0] < static_cast<uint8_t>('0'))
  {
    std::cout << "Invalid seed" << std::endl;
    return EXIT_SUCCESS;
  }*/
  const uint8_t seed_idx = std::stoi(seeds.c_str());

  const bool correctness = options.test_for_correctness();
  constexpr_for<0,seed_table_size,1>([&] (auto seed_idx_constexpr)
  {
    if(seed_idx_constexpr == seed_idx)
    {
      process<Blosum62,seed_idx_constexpr>(multiseq,correctness);
    }
  });

  for (auto &&inputfile : inputfiles)
  {
    std::cout << "# filename\t" << inputfile << std::endl;
  }
  std::cout << "DB statistics" << std::endl;
  for (auto &msg : multiseq->statistics())
  {
    std::cout << "# " << msg << std::endl;
  }
  delete multiseq;
}