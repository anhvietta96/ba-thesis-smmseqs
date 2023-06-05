#include "filter/qgram_environment.hpp"
#include "alignment/blosum62.hpp"
#include "sequences/alphabet.hpp"
#include <iostream>

int main()
{
  constexpr const Blosum62 b62{};
  constexpr const auto char_spec = b62.characters;
  constexpr const auto undefined_rank = b62.num_of_chars;
  constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  
  constexpr const size_t qgram_length = 2;
  constexpr const size_t threshold = 10;
  constexpr const uint8_t const_char = 1;

  constexpr const std::array<uint8_t,2> ref_qgram{ const_char,const_char };
  const QgramEnvironment<Blosum62,qgram_length,threshold> env{ref_qgram};
  const size_t size = env.size();
  for(size_t i = 0; i < size; i++)
  {
    int8_t score = env.env_get(i).first;
    std::cout << (int) score << std::endl;
    auto qgram = env.env_get(i).second;
    for(size_t j = 0; j < qgram_length; j++)
    {
      qgram[j] += 65;
    }
    const std::string output_qgram(std::begin(qgram),std::end(qgram));
    std::cout << output_qgram << std::endl;
  }
}