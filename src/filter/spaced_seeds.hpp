#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <cassert>
#include <iostream>
/*
constexpr size_t get_span(const char *seed, size_t index)
{
    return !(seed[index] - '\0') ? index : get_span(seed,index+1);
}

constexpr size_t get_weight(const char *seed, size_t index, size_t count)
{
    return !(seed[index] - '\0') ? count : !(seed[index] - '0') ? get_weight(seed,index+1,count) : get_weight(seed,index+1,count+1);
}
*/
/*
template<const size_t seed,size_t index>
constexpr size_t get_span(size_t copy)
{
  return (copy == 0) ? index : get_span<seed,index+1>(seed>>1);
}

template<const size_t seed, size_t index>
constexpr size_t get_weight(size_t copy)
{
  return (copy == 0) ? index : ((copy % 2 == 0) ? get_weight<seed,index>(copy>>1) : get_weight<seed,index+1>(copy>>1));
}*/
/*
template<size_t index>
constexpr size_t get_span(const std::bitset<16> seed_bitset)
{
  return (index >= 16) ? 0 : (!seed_bitset[index] ? (size_t(16)-index) : get_span<index+1>(seed_bitset));
}

template<size_t index, size_t count>
constexpr size_t get_weight(const std::bitset<16> seed_bitset)
{
  return (index >= 16) ? count : ((seed_bitset[index] == 1) ? 
  get_weight<index+1,count+1>(seed_bitset) : get_weight<index+1,count>(seed_bitset));
}
*/

constexpr size_t get_span(const std::bitset<16> seed_bitset)
{
  size_t span = 0;
  constexpr_for<0,16,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      span = idx;
    }
  });
  return span+1;
}

template<const size_t span>
constexpr size_t get_weight(const std::bitset<16> seed_bitset)
{
  size_t weight = 0;
  constexpr_for<0,span,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      weight++;
    }
  });
  return weight;
}

template<const char* char_spec,uint8_t undefined_rank,const size_t seed>
class GttlSpacedSeed
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const std::bitset<16> seed_bitset{seed};
  static constexpr const size_t span = get_span(seed_bitset);
  static constexpr const size_t weight = get_weight<span>(seed_bitset);

  public:
  constexpr GttlSpacedSeed(){};

  size_t encode(const char* seq_ptr, const size_t seq_len) const
  {
    //assert(seq_len == span);
    /*
    if(seq_len != span)
    {
      std::cout << seq_len << "      " << span << std::endl;
    }
    */
    size_t code = 0;
    constexpr const size_t alphabet_size = alpha.size();
    constexpr_for<0,span,1>([&] (auto idx)
    {
      if constexpr(seed_bitset[span-1-idx] == 1)
      {
        code *= alphabet_size;
        code += alpha.char_to_rank(seq_ptr[idx]);
      }
    });
    return code;
  }

  constexpr size_t span_get() const
  {
    return span;
  }

  constexpr size_t weight_get() const
  {
    return weight;
  }
};