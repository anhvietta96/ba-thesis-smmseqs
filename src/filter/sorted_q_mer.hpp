#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <algorithm>
#include <array>
#include <iterator>
#include <cassert>
#include <iostream>

template<const char* char_spec,const uint8_t undefined_rank,const size_t qgram_length>
char compare_qgram(const std::array<char,qgram_length> qgram_ref,
const char* qgram, const GttlAlphabet<char_spec,undefined_rank> alpha)
{
  for(size_t idx = 0; idx < qgram_length; idx++)
  {
    if(alpha.char_to_rank(qgram_ref[idx]) > alpha.char_to_rank(qgram[idx]))
    {
      return (char) -1;
    }
    if(alpha.char_to_rank(qgram_ref[idx]) < alpha.char_to_rank(qgram[idx]))
    {
      return (char) 1;
    }
  }
  return 0;
}

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

template<const size_t qgram_length, const size_t ms_size>
static constexpr void constexpr_reverse(std::array<std::array<char,qgram_length>,ms_size> map)
{
  constexpr_for<0,(size_t) ms_size/2,1>([&] (auto qgram_idx)
  {
    constexpr const std::array<char,qgram_length> temp_qgram = map[qgram_idx];
    map[qgram_idx] = map[ms_size-qgram_idx];
    map[ms_size-qgram_idx] = temp_qgram;
  });
}

template<const size_t qgram_length, const size_t alphabet_size,const size_t ms_size>
constexpr void multiset_rec(size_t m, size_t k, std::array<char,qgram_length*ms_size>* map_ptr, size_t* index_ptr, std::array<char,qgram_length> * mslist_ptr)
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
constexpr void enum_multiset_rec(std::array<char,qgram_length*ms_size> * map_ptr)
{
  size_t index = 0;
  std::array<char,qgram_length> mslist = {{ 0 }};
  multiset_rec<qgram_length,alphabet_size,ms_size>(qgram_length,alphabet_size,map_ptr,&index,&mslist);
}

template<const char* char_spec,uint8_t undefined_rank,const size_t qgram_length,const size_t ms_size>
constexpr const std::array<char,qgram_length*ms_size> create_map()
{
  std::array<char,qgram_length*ms_size> map = {{ 0 }};
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
      const char tmp_cc = map[qgram_idx*qgram_length+char_idx];
      map[qgram_idx*qgram_length+char_idx] = map[(ms_size-qgram_idx-1)*qgram_length+char_idx];
      map[(ms_size-qgram_idx-1)*qgram_length+char_idx] = tmp_cc;
    });
  });

  return map;
}

template<const char* char_spec,uint8_t undefined_rank,const size_t qgram_length>
class SortedQmer
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const size_t ms_size = get_multiset_num(qgram_length,undefined_rank);
  static constexpr const std::array<char,qgram_length*ms_size> map
   = create_map<char_spec,undefined_rank,qgram_length,ms_size>();

  public:
  constexpr SortedQmer(){};
  
  std::array<std::array<char,qgram_length>,ms_size> get_map() const
  {
    return map;
  }

  constexpr size_t size_get() const
  {
    return ms_size;
  }

  std::array<char,qgram_length> get_seq_num(size_t i) const
  {
    return map[i];
  }

  constexpr std::array<char,qgram_length> qgram_get(size_t qgram_idx) const
  {
    std::array<char,qgram_length> qgram = {{ 0 }};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = map[qgram_idx*qgram_length+i];
    }
    return qgram;
  }

  std::array<char,qgram_length> extern_qgram_get(size_t qgram_idx) const
  {
    std::array<char,qgram_length> qgram = {{ 0 }};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = alpha.rank_to_char(map[qgram_idx*qgram_length+i]);
    }
    return qgram;
  }

  bool sort(const char* qgram_ptr,char* sorted_qgram, 
  size_t* permutation, size_t sort_length) const
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
    return swapped;
  }  
/*
  size_t encode(const char* qgram) const
  {
    size_t start=0,end=ms_size-1;
    while(start!=end)
    {
      const size_t mid = (start+end)/2;
      const char comp = compare_qgram(map[mid],qgram,alpha);
      if(comp < 0)
      {
        end = mid;
      }
      else if(comp > 0)
      {
        start = mid;
      }
      else
      {
        return mid;
      }
    }
    return start;
  }

  size_t sort_then_encode(const char* qgram)
  {
    return 0;
  }
  
  static char* decode(char* uninit_string, size_t code)
  {
    uninit_string = new char[qgram_length];

  }
  */
};
