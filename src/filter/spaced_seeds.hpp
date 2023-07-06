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

template<const uint8_t span> 
struct SpacedSeed2Qcode{
  uint16_t _code;
  bool _sorted;
  uint8_t _permutation[span];

  SpacedSeed2Qcode(const uint16_t code, const bool sorted, const uint8_t* permutation)
  {
    _code = code;
    _sorted = sorted;
    for(uint8_t i = 0; i < span; i++)
    {
      _permutation[i] = permutation[i];
    }
  }

  uint16_t code() const
  {
    return _code;
  }

  bool sorted() const
  {
    return _sorted;
  }

  uint8_t* permutation() const
  {
    return _permutation;
  }
};

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

template<const uint8_t qgram_length>
struct EncodeInfo {
  size_t code;
  bool sorted;
  std::array<uint8_t,qgram_length> permutation{};

  EncodeInfo(const size_t _code, const bool _sorted, const uint8_t* _permutation)
  {
    code = _code;
    sorted = _sorted;
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      permutation[i] = _permutation[i];
    }
  };
};

template<const char* char_spec,uint8_t undefined_rank,const size_t seed>
class SpacedSeed2SortedCode
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const std::bitset<16> seed_bitset{seed};
  static constexpr const size_t span = get_span(seed_bitset);
  static constexpr const size_t weight = get_weight<span>(seed_bitset);

  public:
  constexpr SpacedSeed2SortedCode(){};

  bool sort(const uint8_t* qgram_ptr,uint8_t* sorted_qgram, 
  uint8_t* permutation, size_t sort_length) const
  {
    for(size_t idx = 0; idx < sort_length; idx++)
    {
      sorted_qgram[idx] = qgram_ptr[idx];
      permutation[idx] = idx;
    }
    bool swapped = false;
    for(size_t pm = 0; pm < sort_length; pm++)
    {
      for(size_t pl = pm; pl > 0 and alpha.char_to_rank(sorted_qgram[pl-1]) 
      > alpha.char_to_rank(sorted_qgram[pl]); pl--)
      //for(size_t pl = pm; pl > 0 and sorted_qgram[pl-1] > sorted_qgram[pl]; pl--)
      {
        const uint8_t tmp_cc = sorted_qgram[pl-1];
        sorted_qgram[pl-1] = sorted_qgram[pl];
        sorted_qgram[pl] = tmp_cc;

        const size_t tmp_t = permutation[pl-1];
        permutation[pl-1] = permutation[pl];
        permutation[pl] = tmp_t;

        swapped = true;
      }
    }
    return !swapped;
  }

  EncodeInfo<weight> encode(const char* seq_ptr) const
  {
    size_t code = 0;
    constexpr const auto alphabet_size = alpha.size();
    
    uint8_t ref_qgram[weight];
    uint8_t sorted_qgram[weight];
    uint8_t permutation[weight];
    uint8_t i = 0;
    
    constexpr_for<0,span,1>([&] (auto idx)
    {
      if constexpr(seed_bitset[span-1-idx] == 1)
      {
        ref_qgram[i] = static_cast<uint8_t>(seq_ptr[idx]);
        i++;
      }
    });

    const bool sorted = sort(ref_qgram,sorted_qgram,permutation,weight);

    constexpr_for<0,weight,1>([&] (auto idx)
    {
      code *= alphabet_size;
      code += alpha.char_to_rank(sorted_qgram[idx]);
    });

    return EncodeInfo<weight>(code,sorted,permutation);
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
