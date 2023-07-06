#include<iostream>
#include"filter/env_constructor.hpp"
#include"alignment/blosum62.hpp"
#include"filter/spaced_seeds.hpp"

int main()
{
  constexpr const uint16_t sq_code_1 = 1683;
  constexpr const uint16_t sq_code_2 = 1857;
  constexpr const int8_t threshold = 30;
/*
  const Env4Constructor_1<Blosum62,sq_code_1,sq_code_2,threshold> env4{};

  constexpr const Blosum62 sc{};
  constexpr const auto char_spec = sc.characters;
  constexpr const auto undefined_rank = sc.num_of_chars;
  constexpr const UnsortedQmer<char_spec,undefined_rank,2> uq{};

  std::cout << "Size = " << (int)env4.size() << std::endl;

  for(uint16_t idx = 0; idx < env4.size(); idx++)
  {
    const auto elem = env4.elem_get(idx);
    const auto uq_1 = uq.extern_qgram_get(std::get<1>(elem));
    const std::string output_qgram_1(std::begin(uq_1),std::end(uq_1));
    const auto uq_2 = uq.extern_qgram_get(std::get<2>(elem));
    const std::string output_qgram_2(std::begin(uq_2),std::end(uq_2));
    std::cout << (int) std::get<0>(elem) << "\t" << output_qgram_1 << output_qgram_2 << std::endl;
  }*/

  //static constexpr const uint8_t qgram[] = {'L','H','G','P'};
#define MER4
#ifdef MER4
  const Env4Constructor_2<Blosum62,sq_code_1,sq_code_2,threshold> env4{};
  
  constexpr const Blosum62 sc{};
  constexpr const auto char_spec = sc.character_spec;
  constexpr const auto undefined_rank = sc.num_of_chars;
  constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};

  std::cout << "Size = " << (int)env4.size() << std::endl;
  for(uint16_t idx = 0; idx < env4.size(); idx++)
  {
    const auto elem = env4.elem_get(idx);
    std::cout << (int) std::get<0>(elem) << "\t";
    const auto uqgram = std::get<1>(elem);
    for(uint8_t idx = 0; idx < 4; idx++)
    {
      std::cout << alpha.rank_to_char(uqgram[idx]);
    }
    std::cout << std::endl;
  }
#endif
#ifdef MER6
  const Env6Constructor_2<Blosum62,sq_code_1,sq_code_2,threshold> env6{};
  
  constexpr const Blosum62 sc{};
  constexpr const auto char_spec = sc.character_spec;
  constexpr const auto undefined_rank = sc.num_of_chars;
  constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  constexpr const SortedQmer<char_spec,undefined_rank,3> sorted_q{};

  const auto sqgram_1 = sorted_q.extern_qgram_get(sq_code_1);
  const std::string output_sqgram_1(std::begin(sqgram_1),std::end(sqgram_1));
  std::cout << "Qgram 1: " << output_sqgram_1 << std::endl;

  const auto sqgram_2 = sorted_q.extern_qgram_get(sq_code_2);
  const std::string output_sqgram_2(std::begin(sqgram_2),std::end(sqgram_2));
  std::cout << "Qgram 2: " << output_sqgram_2 << std::endl;

  std::cout << "Threshold: " << (int)threshold << std::endl;

  std::cout << "Size = " << (int)env6.size() << std::endl;
  for(uint16_t idx = 0; idx < env6.size(); idx++)
  {
    const auto elem = env6.elem_get(idx);
    std::cout << (int) std::get<0>(elem) << "\t";
    const auto uqgram = std::get<1>(elem);
    for(uint8_t idx = 0; idx < 6; idx++)
    {
      std::cout << alpha.rank_to_char(uqgram[idx]);
    }
    std::cout << std::endl;
  }
#endif
}