#include "sequences/alphabet.hpp"
#include <cassert>
#include <iostream>

size_t calc_weight(const size_t span, const char *seed)
{
  size_t weight = 0;
  std::cout << span << std::endl;
  for(size_t i = 0; i < span; i++)
  {
    if(seed[i]-'0')
    {
      weight++;
    }
  }
  return weight;
}

template<const char* char_spec,uint8_t undefined_rank,const char *seed,const size_t seed_len>
class GttlSpacedSeed
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const size_t span = seed_len-1;
  static constexpr const size_t weight = calc_weight(span,seed);

  public:
  GttlSpacedSeed(){};

  size_t encode(const char* seq_ptr, const size_t seq_len)
  {
    assert((seq_len == this->span));
    size_t code = 0;
    size_t alphabet_size = this->alpha.size();
    for(size_t idx = 0; idx < span; idx++)
    {
      if(seed[idx]-'0')
      {
        code *= alphabet_size;
        code += this->alpha.char_to_rank(seq_ptr[idx]);
      }
    }

    return code;
  }
};