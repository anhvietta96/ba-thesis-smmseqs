#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <algorithm>
#include <array>
#include <iterator>
#include <cassert>
#include <iostream>

constexpr size_t binom(size_t a,size_t b)
{
  if (b == 0 or a == b)
  {
    return size_t(1);
  }
  if (2*b > a and a != b)
  {
    return binom(a,a-b);
  }
  return binom(a-1,b-1) + binom(a-1,b);
}

constexpr size_t get_multiset_num(const size_t qgram_length, const size_t alphabet_size)
{
  return binom(qgram_length + alphabet_size -1,qgram_length);
}

template<const size_t qgram_length, const size_t alphabet_size,const size_t ms_size>
constexpr void multiset_rec(size_t m, size_t k, std::array<uint8_t,qgram_length*ms_size>* map_ptr, size_t* index_ptr, std::array<uint8_t,qgram_length> * mslist_ptr)
{
  if(m == 0)
  {
    //std::initializer_list
    constexpr_for<0,(size_t) qgram_length,1>([&] (auto char_idx)
    {
      (*map_ptr)[(*index_ptr)*qgram_length+char_idx] = (*mslist_ptr)[char_idx];
    });
    (*index_ptr)++;
  }
  else if(k > 0)
  {
    for(size_t take = 0; take < m+1; take++)
    {
      for(size_t d = qgram_length-m; d < qgram_length-m+take; d++)
      {
        (*mslist_ptr)[d] = alphabet_size - k;
      }
      multiset_rec<qgram_length,alphabet_size,ms_size>(m-take,k-1,map_ptr,index_ptr,mslist_ptr);
    }
  }
}

template<const size_t qgram_length, const size_t alphabet_size,const size_t ms_size>
constexpr void enum_multiset_rec(std::array<uint8_t,qgram_length*ms_size> * map_ptr)
{
  size_t index = 0;
  std::array<uint8_t,qgram_length> mslist = {{ 0 }};
  multiset_rec<qgram_length,alphabet_size,ms_size>(qgram_length,alphabet_size,map_ptr,&index,&mslist);
}

template<const char* char_spec,uint8_t undefined_rank,const size_t qgram_length,const size_t ms_size>
constexpr const std::array<uint8_t,qgram_length*ms_size> create_map()
{
  std::array<uint8_t,qgram_length*ms_size> map = {{ 0 }};
  enum_multiset_rec<qgram_length,undefined_rank,ms_size>(&map);
  
  /*constexpr_for<0,ms_size,1>([&] (auto qgram_idx)
  {
    constexpr_for<0,qgram_length,1>([&] (auto cc)
    {
      map[qgram_idx][cc] = alpha.rank_to_char(map[qgram_idx][cc]);
    });
  });*/
  //constexpr_reverse(map);
  //std::ranges::reverse(std::begin(map),std::end(map));

  constexpr_for<0,(size_t) ms_size/2,1>([&] (auto qgram_idx)
  {
    constexpr_for<0,qgram_length,1>([&] (auto char_idx)
    {
      const uint8_t tmp_cc = map[qgram_idx*qgram_length+char_idx];
      map[qgram_idx*qgram_length+char_idx] = map[(ms_size-qgram_idx-1)*qgram_length+char_idx];
      map[(ms_size-qgram_idx-1)*qgram_length+char_idx] = tmp_cc;
    });
  });

  return map;
}

constexpr size_t power_1(const size_t& base, const size_t& exponent)
{
  return (exponent != 1) ? (base*power_1(base,exponent-1)) : base;
}

template<const char* char_spec,const size_t undefined_rank,const uint8_t qgram_length>
class SortedQmer
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const uint16_t ms_size = get_multiset_num(qgram_length,undefined_rank);
  std::array<uint8_t,qgram_length*ms_size> map = {{ 0 }};
  std::array<uint16_t,power_1(undefined_rank,qgram_length)> unsorted_to_sorted_map{};
  
  public:
  constexpr SortedQmer()
  {
    enum_multiset_rec<qgram_length,undefined_rank,ms_size>(&map);
    //Reverse map
    /*constexpr_for<0,(size_t) ms_size/2,1>([&] (auto qgram_idx)
    {
      constexpr_for<0,qgram_length,1>([&] (auto char_idx)
      {
        const uint8_t tmp_cc = map[qgram_idx*qgram_length+char_idx];
        map[qgram_idx*qgram_length+char_idx] = map[(ms_size-qgram_idx-1)*qgram_length+char_idx];
        map[(ms_size-qgram_idx-1)*qgram_length+char_idx] = tmp_cc;
      });
    });*/
    for(uint16_t qgram_idx = 0; qgram_idx < ms_size/2; qgram_idx++)
    {
      for(uint8_t char_idx = 0; char_idx < qgram_length; char_idx++)
      {
        const uint8_t tmp_cc = map[qgram_idx*qgram_length+char_idx];
        map[qgram_idx*qgram_length+char_idx] = map[(ms_size-qgram_idx-1)*qgram_length+char_idx];
        map[(ms_size-qgram_idx-1)*qgram_length+char_idx] = tmp_cc;
      }
    }
    
    constexpr const auto alphasize = alpha.size();
    size_t code = 0;
    constexpr_for<0,(size_t) ms_size,1>([&] (auto qgram_idx){
      constexpr_for<0,qgram_length,1>([&] (auto char_idx){
        code *= alphasize;
        code += map[qgram_idx*qgram_length+char_idx];
      });
      unsorted_to_sorted_map[code] = qgram_idx;
      code = 0;
    });
  };

  constexpr uint16_t sorted_code_get(const uint16_t unsorted_code) const
  {
    if(unsorted_code != 0 and unsorted_to_sorted_map[unsorted_code] == 0)
    {
      //false calculation
      return UINT16_MAX;
    }
    return unsorted_to_sorted_map[unsorted_code];
  }
  
  constexpr size_t size_get() const
  {
    return ms_size;
  }

  constexpr std::array<uint8_t,qgram_length> qgram_get(const size_t qgram_idx) const
  {
    std::array<uint8_t,qgram_length> qgram{};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = map[qgram_idx*qgram_length+i];
    }
    return qgram;
  }

  std::array<uint8_t,qgram_length> extern_qgram_get(const size_t qgram_idx) const
  {
    std::array<uint8_t,qgram_length> qgram{};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = alpha.rank_to_char(map[qgram_idx*qgram_length+i]);
    }
    return qgram;
  }
};
