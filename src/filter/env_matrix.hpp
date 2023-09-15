#include <iostream>
#include "sorted_qgram.hpp"
#include "unsorted_qgram.hpp"


template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
constexpr std::array<int8_t,qgram_length> create_threshold_arr()
{
  constexpr const ScoreClass sc{};
  std::array<int8_t,qgram_length> threshold_arr{};
  
  constexpr const int8_t max_char_score = sc.highest_score;
/*
  constexpr_for<0,qgram_length,1>([&] (auto char_idx)
  {
    threshold_arr[char_idx] = static_cast<int8_t>(threshold - (qgram_length-char_idx-1) * max_char_score);
  });*/
  for(uint8_t char_idx = 0; char_idx < qgram_length; char_idx++)
  {
    threshold_arr[char_idx] = static_cast<int8_t>(threshold - (qgram_length-char_idx-1) * max_char_score);
  }

  return threshold_arr;
}

template<const uint8_t qgram_length, const size_t undefined_rank>
constexpr std::array<size_t,qgram_length> create_pow_arr()
{
  std::array<size_t, qgram_length> pow_arr{};
  pow_arr[0] = undefined_rank;
  for(uint8_t q_idx = 1; q_idx < qgram_length; q_idx++)
  {
    pow_arr[q_idx] = pow_arr[q_idx-1]*undefined_rank;
  }
  return pow_arr;
}

template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
struct EnvMatrix {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = power(undefined_rank,qgram_length);
  static constexpr const size_t matrix_size = sorted_size*unsorted_size;

  std::array<int8_t,matrix_size> matrix{};
  //int8_t matrix[matrix_size]{};
  std::array<uint16_t,sorted_size> size_arr{};

  constexpr EnvMatrix()
  {
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();
    std::array<uint8_t,qgram_length> sorted_qgram{};
    std::array<uint8_t,qgram_length> unsorted_qgram{};

    for(uint16_t sorted_idx = 0; sorted_idx < 100; sorted_idx++)
    //constexpr_for<0,2600,1>([&] (uint16_t sorted_idx)
    {
      sorted_q.get_qgram(sorted_idx,sorted_qgram);
      const size_t sorted_pos = sorted_idx*unsorted_size;
      
      for(uint16_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      //constexpr_for<0,unsorted_size,1>([&] (uint16_t unsorted_idx)
      {
        unsorted_q.get_qgram(unsorted_idx,unsorted_qgram);
        const size_t idx = sorted_pos + unsorted_idx;
        int8_t score = 0;

        for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
        //constexpr_for<0,qgram_length,1>([&] (uint8_t q_idx)
        {
          if(score != INT8_MIN)
          {
            score += sc.score_matrix[sorted_qgram[q_idx]][unsorted_qgram[q_idx]];
            if(score < threshold_arr[q_idx])
            {
              score = INT8_MIN;
            }
          }
        }
        
        matrix[idx] = score;
      }
    }
  };
};