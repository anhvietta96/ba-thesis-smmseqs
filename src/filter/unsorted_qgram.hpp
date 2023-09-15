#ifndef UNSORTED_QGRAM_HPP
#define UNSORTED_QGRAM_HPP
#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include "filter/utils.hpp"

template<const size_t undefined_rank, const size_t qgram_length, const size_t map_size>
constexpr std::array<char,qgram_length * map_size> map_create()
  {
    char tmp_qgram[qgram_length] = {{ 0 }};
    std::array<char,qgram_length * map_size> map = {{ 0 }};
    constexpr_for<0,map_size,1>([&] (auto qgram_idx)
    {
      constexpr_for<0,qgram_length,1>([&] (auto char_idx)
      {
        map[qgram_idx*qgram_length+char_idx] = tmp_qgram[char_idx];
      });
      tmp_qgram[qgram_length-1] += 1; 
      if(qgram_idx != map_size-1 and tmp_qgram[qgram_length-1] == undefined_rank)
      {
        for(int16_t i = qgram_length-1; i >= 0; i--)
        {
          if(tmp_qgram[i] < undefined_rank-1)
          {
            tmp_qgram[i]++;
            break;
          }
          else
          {
            tmp_qgram[i]=0;
          }
        }
      }
    });
    return map;
  }

template<const char* char_spec, const size_t undefined_rank, const size_t qgram_length>
class UnsortedQmer {
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const size_t map_size = constexpr_pow(undefined_rank,qgram_length);
  std::array<uint8_t,qgram_length * map_size> map = {{ 0 }};

  public:
  constexpr UnsortedQmer()
  {
    uint8_t tmp_qgram[qgram_length]{};
    for(size_t qgram_idx = 0; qgram_idx < map_size; qgram_idx++)
    {
      for(size_t char_idx = 0; char_idx < qgram_length; char_idx++)
      {
        map[qgram_idx*qgram_length+char_idx] = tmp_qgram[char_idx];
      }

      tmp_qgram[qgram_length-1] += 1;
      if(qgram_idx != map_size-1 and tmp_qgram[qgram_length-1] == undefined_rank)
      {
        for(int16_t i = qgram_length-1; i >= 0; i--)
        {
          if(static_cast<uint8_t>(tmp_qgram[i]) < undefined_rank-1)
          {
            tmp_qgram[i]++;
            break;
          }
          else
          {
            tmp_qgram[i]=0;
          }
        }
      }
    }
  };

  constexpr size_t size_get() const
  {
    return map_size;
  }
/*
  void qgram_get(size_t qgram_idx,char* qgram_ptr) const
  {
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram_ptr[i] = map[qgram_idx*qgram_length+i];
    }
  }*/

  constexpr std::array<uint8_t,qgram_length> qgram_get(size_t qgram_idx) const
  {
    std::array<uint8_t,qgram_length> qgram{};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = map[qgram_idx*qgram_length+i];
    }
    return qgram;
  }

  constexpr void get_qgram(size_t qgram_idx,std::array<uint8_t,qgram_length>& qgram) const
  {
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = map[qgram_idx*qgram_length+i];
    }
  }

  std::array<uint8_t,qgram_length> extern_qgram_get(size_t qgram_idx) const
  {
    std::array<uint8_t,qgram_length> qgram{};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = alpha.rank_to_char(map[qgram_idx*qgram_length+i]);
    }
    return qgram;
  }
};
#endif
