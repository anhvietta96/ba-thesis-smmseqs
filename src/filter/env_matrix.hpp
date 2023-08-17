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

template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
struct EnvMatrix {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const size_t sorted_size = sorted_q.size_get();
  static constexpr const size_t unsorted_size = power(undefined_rank,qgram_length);
  static constexpr const size_t matrix_size = sorted_size*unsorted_size;

  //std::array<int8_t,total_arr_size> matrix{};
  int8_t matrix[matrix_size]{};
  std::array<uint16_t,sorted_size> size_arr{};

  EnvMatrix()
  {
    std::cout << "Init" << std::endl;
    /*
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();

    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      size_arr[sorted_idx] = unsorted_size;
      const std::array<uint8_t,qgram_length> sorted_qgram = sorted_q.qgram_get(sorted_idx);
      
      for(uint16_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      {
        const std::array<uint8_t,qgram_length> unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
        for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
        {
          if(matrix[sorted_idx*unsorted_size+unsorted_idx] != INT8_MIN)
          {
            matrix[sorted_idx*unsorted_size+unsorted_idx] += sc.score_matrix[sorted_qgram[q_idx]][unsorted_qgram[q_idx]];
            if((int8_t)matrix[sorted_idx*unsorted_size+unsorted_idx] < (int8_t)threshold_arr[q_idx])
            {
              matrix[sorted_idx*unsorted_size+unsorted_idx] = INT8_MIN;
              size_arr[sorted_idx]--;
            }
          }
        }
      }
    }*/
  };
};