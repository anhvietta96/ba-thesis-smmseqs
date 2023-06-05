#include "filter/sorted_q_mer.hpp"
#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <iostream>

constexpr size_t power(const size_t base, const size_t exponent)
{
  return (exponent != 1) ? (base*power(base,exponent-1)) : base;
}

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
        for(size_t i = qgram_length-1; i >= 0; i--)
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
  static constexpr const size_t map_size = power(undefined_rank,qgram_length);
  std::array<uint8_t,qgram_length * map_size> map = {{ 0 }};

  public:
  constexpr UnsortedQmer()
  {
    uint8_t tmp_qgram[qgram_length] = {{ 0 }};
    for(size_t qgram_idx = 0; qgram_idx < map_size; qgram_idx++)
    {
      for(size_t char_idx = 0; char_idx < qgram_length; char_idx++)
      {
        map[qgram_idx*qgram_length+char_idx] = tmp_qgram[char_idx];
      }

      tmp_qgram[qgram_length-1] += 1;
      if(qgram_idx != map_size-1 and tmp_qgram[qgram_length-1] == undefined_rank)
      {
        for(size_t i = qgram_length-1; i >= 0; i--)
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
    std::array<uint8_t,qgram_length> qgram = {{ 0 }};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = map[qgram_idx*qgram_length+i];
    }
    return qgram;
  }

  std::array<uint8_t,qgram_length> extern_qgram_get(size_t qgram_idx) const
  {
    std::array<uint8_t,qgram_length> qgram = {{ 0 }};
    for(size_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = alpha.rank_to_char(map[qgram_idx*qgram_length+i]);
    }
    return qgram;
  }
};

template<class ScoreClass,const char* char_spec,const size_t undefined_rank,
const size_t qgram_length,const size_t sorted_size,const size_t unsorted_size>

constexpr const std::array<int16_t,sorted_size*unsorted_size> score_mat_create(
  const SortedQmer<char_spec,undefined_rank,qgram_length> SortedQgram,
  const UnsortedQmer<char_spec,undefined_rank,qgram_length> UnsortedQgram)
{
  constexpr const ScoreClass sc{};
  std::array<int16_t,sorted_size*unsorted_size> matrix = {{ 0 }};
  
  constexpr_for<0,sorted_size,1>([&] (auto sorted_idx)
  {
    constexpr const std::array<uint8_t,qgram_length> sorted_qgram = SortedQgram.qgram_get(sorted_idx);
    constexpr_for<0,unsorted_size,1>([&] (auto unsorted_idx)
    {
      constexpr const std::array<uint8_t,qgram_length> unsorted_qgram = UnsortedQgram.qgram_get(unsorted_idx);
      constexpr_for<0,qgram_length,1>([&] (auto char_idx)
      {
        if(sorted_qgram[char_idx] < undefined_rank and unsorted_qgram[char_idx] < undefined_rank)
        {
          matrix[sorted_idx*unsorted_size+unsorted_idx] += (int16_t) sc.score_matrix[sorted_qgram[char_idx]][unsorted_qgram[char_idx]];
        }
      });
    });
  });
  
  return matrix;
}

template<class ScoreClass,const size_t qgram_length>
class QgramScoreMatrix {
  private:
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = static_cast<size_t>(sc.num_of_chars);
  
  static constexpr const SortedQmer<sc.characters,undefined_rank,qgram_length> SortedQgram{};
  static constexpr const size_t sorted_size = SortedQgram.size_get();
  static constexpr const UnsortedQmer<sc.characters,undefined_rank,qgram_length> UnsortedQgram{};
  static constexpr const size_t unsorted_size = UnsortedQgram.size_get();

  std::array<int16_t,sorted_size*unsorted_size> matrix = {{ 0 }};

  public:
  constexpr QgramScoreMatrix()
  {
    constexpr_for<0,sorted_size,1>([&] (auto sorted_idx)
    {
      const std::array<uint8_t,qgram_length> sorted_qgram = SortedQgram.qgram_get(sorted_idx);
      constexpr_for<0,unsorted_size,1>([&] (auto unsorted_idx)
      {
        const std::array<uint8_t,qgram_length> unsorted_qgram = UnsortedQgram.qgram_get(unsorted_idx);
        constexpr_for<0,qgram_length,1>([&] (auto char_idx)
        {
          if(sorted_qgram[char_idx] < undefined_rank and unsorted_qgram[char_idx] < undefined_rank)
          {
            matrix[sorted_idx*unsorted_size+unsorted_idx] += (int16_t) sc.score_matrix[sorted_qgram[char_idx]][unsorted_qgram[char_idx]];
          }
        });
      });
    });
  };

  int16_t score_get(size_t sorted_qgram_code,size_t unsorted_qgram_code) const
  {
    return matrix[sorted_qgram_code * unsorted_size + unsorted_qgram_code];
  }

  constexpr size_t sorted_size_get() const
  {
    return sorted_size;
  }

  constexpr size_t unsorted_size_get() const
  {
    return unsorted_size;
  }
};