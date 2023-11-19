#ifndef ENV_MATRIX
#define ENV_MATRIX
#include <iostream>
#include "filter/sorted_qgram.hpp"
#include "filter/unsorted_qgram.hpp"
#include "filter/utils.hpp"
/*
template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
constexpr std::array<int8_t,qgram_length> create_threshold_arr()
{
  constexpr const ScoreClass sc{};
  std::array<int8_t,qgram_length> threshold_arr{};
  
  constexpr const int8_t max_char_score = sc.highest_score;

  constexpr_for<0,qgram_length,1>([&] (auto char_idx)
  {
    threshold_arr[char_idx] = static_cast<int8_t>(threshold - (qgram_length-char_idx-1) * max_char_score);
  });
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
};*/
/*
template<class ScoreClass,const uint8_t qgram_length>
struct EnvMatrix {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = constexpr_pow(undefined_rank,qgram_length);
  //static constexpr const size_t matrix_size = sorted_size*unsorted_size;

  //std::array<int8_t,matrix_size> matrix{};
  std::array<std::array<ScoreQgramcodePair,unsorted_size>,sorted_size> matrix{};

  constexpr EnvMatrix()
  {
    //for(uint64_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    constexpr_for<0,sorted_size,1>([&] (uint16_t sorted_idx)
    {
      const auto sorted_qgram = sorted_q.qgram_get(sorted_idx);
      //for(uint64_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      constexpr_for<0,unsorted_size,1>([&] (uint16_t unsorted_idx)
      {
        const auto unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
        int8_t score = 0;
        //for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
        constexpr_for<0,qgram_length,1>([&] (uint8_t q_idx)
        {
          score += sc.score_matrix[sorted_qgram[q_idx]][unsorted_qgram[q_idx]];
        });
        matrix[sorted_idx*unsorted_size+unsorted_idx] = score;
      });
    });
  };
};*/

struct ScoreQgramcodePair2{
  int8_t score = 0;
  uint16_t code = 0;

  constexpr ScoreQgramcodePair2()
  : score { 0 }, code { 0 }{};

  constexpr ScoreQgramcodePair2(const int8_t& _score, const uint16_t& _code)
  : score { _score }, code { _code }{};

  constexpr ScoreQgramcodePair2(const ScoreQgramcodePair2& other_pair)
  : score { other_pair.score }, code { other_pair.code }{};

  constexpr void operator= (const ScoreQgramcodePair2& other_pair)
  {
    score = other_pair.score;
    code = other_pair.code;
  };

  bool operator<(const ScoreQgramcodePair2& other_pair){
    if(score > other_pair.score){
      return true;
    }
    return false;
  }
};

template<class ScoreClass,const uint8_t qgram_length>
struct EnvMatrix2 {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = constexpr_pow(undefined_rank,qgram_length);
  //static constexpr const size_t matrix_size = sorted_size*unsorted_size;

  //std::array<int8_t,matrix_size> matrix{};
  using Environment = std::vector<ScoreQgramcodePair2>;
  std::array<Environment,unsorted_size> matrix{};

  EnvMatrix2()
  {
    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    //constexpr_for<0,sorted_size,1>([&] (uint16_t sorted_idx)
    {
      const uint8_t* const sorted_qgram = sorted_q.qgram_get(sorted_idx);
      for(uint16_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      //constexpr_for<0,unsorted_size,1>([&] (uint16_t unsorted_idx)
      {
        const uint8_t* const unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
        int8_t score = 0;
        for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
        //constexpr_for<0,qgram_length,1>([&] (uint8_t q_idx)
        {
          score += sc.score_matrix[sorted_qgram[q_idx]][unsorted_qgram[q_idx]];
        }
        matrix[sorted_idx].push_back(ScoreQgramcodePair2(score,unsorted_idx));
      }
      std::sort(matrix[sorted_idx].begin(),matrix[sorted_idx].end());
    }
  };

  const ScoreQgramcodePair2* const sorted_env_get(const uint64_t sq_code) const {
    return matrix[sq_code].data();
  }
};

template<class ScoreClass,const uint8_t qgram_length>
struct FullMatrix {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const uint16_t unsorted_size = unsorted_q.size_get();
  
  using Environment = std::vector<ScoreQgramcodePair2>;
  std::array<Environment,unsorted_size> matrix{};

  FullMatrix(){
    for(size_t idx1 = 0; idx1 < unsorted_size; idx1++){
      const uint8_t* const qgram1 = unsorted_q.qgram_get(idx1);
      for(uint64_t idx2 = 0; idx2 < unsorted_size; idx2++){
        const uint8_t* const qgram2 = unsorted_q.qgram_get(idx2);
        int8_t score = 0;
        for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++){
          score += sc.score_matrix[qgram1[q_idx]][qgram2[q_idx]];
        }
        matrix[idx1].push_back(ScoreQgramcodePair2(score,idx2));
      }
      std::sort(matrix[idx1].begin(),matrix[idx1].end());
    }
  };

  const ScoreQgramcodePair2* const sorted_env_get(const uint64_t sq_code) const {
    return matrix[sq_code].data();
  }
};
#endif