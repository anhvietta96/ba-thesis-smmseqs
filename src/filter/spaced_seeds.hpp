#ifndef SPACED_SEED_HPP
#define SPACED_SEED_HPP
#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include "filter/multiset_code.hpp"
#include "filter/distribution.hpp"

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

template<const uint8_t qgram_length, const uint8_t num_of_primary_env>
struct EncodeInfo {
  std::array<uint8_t,qgram_length> permutation;
  bool sorted;
  std::array<uint16_t,num_of_primary_env> codes;

  EncodeInfo(std::array<uint8_t,qgram_length> _permutation,const bool& _sorted,
              std::array<uint16_t,num_of_primary_env> _codes): 
              permutation(_permutation),sorted(_sorted),codes(_codes){};
};

template<const uint8_t num_of_primary_env,const uint8_t weight>
constexpr std::array<uint8_t,num_of_primary_env> create_subqgram_length_arr(){
  std::array<uint8_t,num_of_primary_env> qgram_length_arr{};
  for(uint8_t i = 0; i < num_of_primary_env; i++){
    qgram_length_arr[i] = 3;
  }

  constexpr const uint8_t overcount = (num_of_primary_env * 3) - weight;
  for(uint8_t i = 0; i < overcount; i++){
    qgram_length_arr[i]--;
  }

  return qgram_length_arr;
}

template<const uint8_t num_of_primary_env,const uint8_t weight>
constexpr std::array<int8_t,num_of_primary_env> create_env_threshold_arr(const std::array<uint8_t,num_of_primary_env>& qgram_length_arr,const int8_t& threshold,const int8_t& highest_score){
  std::array<int8_t,num_of_primary_env> env_threshold_arr{};
  for(uint8_t i = 0; i < num_of_primary_env; i++){
    env_threshold_arr[i] = threshold - (weight - qgram_length_arr[i]) * highest_score;
  }
  return env_threshold_arr;
}

template<class ScoreClass,const size_t seed, const uint8_t max_subqgram_length>
class SpacedSeedEncoder {
  private:
  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const std::bitset<16> seed_bitset{seed};
  static constexpr const size_t span = get_span(seed_bitset);
  static constexpr const size_t weight = get_weight<span>(seed_bitset);
  
  static constexpr const uint8_t num_of_primary_env = weight % max_subqgram_length ? (weight / max_subqgram_length + 1) : 
                                                      weight / max_subqgram_length;
  static constexpr const MultisetEncoder<char_spec,undefined_rank,max_subqgram_length> multiset_encoder{};
  
  //Threshold
  static constexpr const Distribution<ScoreClass,weight-1> dis{};
  static constexpr const int8_t threshold = dis.threshold_get(weight-1);

  //Init in constexpr constructor
  std::array<uint8_t,num_of_primary_env> subqgram_length_arr{};
  std::array<int8_t,num_of_primary_env> env_threshold_arr{};

  public:
  constexpr SpacedSeedEncoder(){
    uint8_t overcount = (num_of_primary_env * max_subqgram_length) - weight;
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      subqgram_length_arr[i] = max_subqgram_length;
    }
    
    while(overcount != 0){
      const uint8_t _overcount = overcount;
      for(uint8_t i = 0; i < _overcount; i++){
        subqgram_length_arr[num_of_primary_env-1-i]--;
        overcount--;
      }
    }
    
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      env_threshold_arr[i] = threshold - (weight - subqgram_length_arr[i]) * sc.highest_score;
    }
  };

  bool sort(const uint8_t* qgram_ptr,uint8_t* sorted_qgram, 
            uint8_t* permutation) const
  {
    for(size_t idx = 0; idx < weight; idx++)
    {
      sorted_qgram[idx] = qgram_ptr[idx];
      permutation[idx] = idx;
    }
    bool swapped = false;
    for(size_t pm = 0; pm < weight; pm++)
    {
      for(size_t pl = pm; pl > 0 and sorted_qgram[pl-1] > sorted_qgram[pl]; pl--)
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

  std::array<size_t,num_of_primary_env> encode(const char* seq, uint8_t* qgram, uint8_t* permutation, bool& sorted) const {
    std::array<size_t,num_of_primary_env> sorted_qgram_codes{};
    uint8_t extracted_qgram[weight];
    uint8_t qgram_idx = 0;
    
    constexpr_for<0,span,1>([&] (auto idx)
    {
      if constexpr(seed_bitset[span-1-idx] == 1)
      {
        extracted_qgram[qgram_idx] = static_cast<uint8_t>(seq[idx]);
        qgram_idx++;
      }
    });

    sorted = sort(extracted_qgram,qgram,permutation);

    qgram_idx = 0;
    size_t code;
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      code = 0;
      for(uint8_t j = 0; j < subqgram_length_arr[i]; j++){
        code += multiset_encoder.relative_encode(subqgram_length_arr[i],j,qgram[qgram_idx]);
        qgram_idx++;
      }
      sorted_qgram_codes[i] = code;
    }

    return sorted_qgram_codes;
  }

  std::array<size_t,num_of_primary_env> encode_unsorted(const char* seq) const {
    std::array<size_t,num_of_primary_env> sorted_qgram_codes{};
    
    size_t code = 0;
    uint8_t extracted_qgram[weight];
    uint8_t qgram_idx = 0;
    
    constexpr_for<0,span,1>([&] (auto idx)
    {
      if constexpr(seed_bitset[span-1-idx] == 1)
      {
        extracted_qgram[qgram_idx] = static_cast<uint8_t>(seq[idx]);
        qgram_idx++;
      }
    });

    qgram_idx = 0;
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      code = 0;
      //std::cout << (int)subqgram_length_arr[i] << std::endl;
      for(uint8_t j = 0; j < subqgram_length_arr[i]; j++){
        code *= undefined_rank;
        code += extracted_qgram[qgram_idx];
        qgram_idx++;
        //std::cout << code << std::endl;
      }
      sorted_qgram_codes[i] = code;
    }

    return sorted_qgram_codes;
  }

  constexpr const std::array<uint8_t,num_of_primary_env>& subqgram_length_arr_get() const {
    return subqgram_length_arr;
  }

  constexpr const std::array<int8_t,num_of_primary_env>& env_threshold_arr_get() const {
    return env_threshold_arr;
  }

  constexpr uint8_t num_of_primary_env_get() const {
    return num_of_primary_env;
  }

  constexpr int8_t threshold_get() const {
    return threshold;
  }

  constexpr size_t span_get() const {
    return span;
  }

  constexpr size_t weight_get() const {
    return weight;
  }

  constexpr const std::bitset<16>& seed_bitset_get() const {
    return seed_bitset;
  }
};
#endif