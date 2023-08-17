#include "filter/unsorted_qgram.hpp"
#include "filter/sorted_qgram.hpp"

#include <vector>
#include<iostream>

/*sort matrix first*/

struct ScoreQgramcodePair{
  int8_t score = 0;
  uint16_t code = 0;

  constexpr ScoreQgramcodePair()
  : score { 0 }, code { 0 }{};

  constexpr ScoreQgramcodePair(const int8_t& _score, const uint16_t& _code)
  : score { _score }, code { _code }{};

  constexpr ScoreQgramcodePair(const ScoreQgramcodePair& other_pair)
  : score { other_pair.score }, code { other_pair.code }{};

  constexpr void operator= (const ScoreQgramcodePair& other_pair)
  {
    score = other_pair.score;
    code = other_pair.code;
  };
};

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

template<class ScoreClass,const size_t qgram_length, const int8_t threshold>
class QgramEnvironment {
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = static_cast<size_t>(sc.num_of_chars);
  static constexpr const UnsortedQmer<sc.character_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const SortedQmer<sc.character_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = unsorted_q.size_get();

  using ScoreQgramList = std::vector<ScoreQgramcodePair>;
  std::array<ScoreQgramList,sorted_size> environment{};

  public:
  constexpr QgramEnvironment()
  {
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();

    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      const std::array<uint8_t,qgram_length> sorted_qgram = sorted_q.qgram_get(sorted_idx);
      for(uint16_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      {
        const std::array<uint8_t,qgram_length> unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
        int8_t score = 0;
        for(size_t char_idx = 0; char_idx < qgram_length; char_idx++)
        {
          if(sorted_qgram[char_idx] < undefined_rank and unsorted_qgram[char_idx] < undefined_rank)
          {
            score += sc.score_matrix[sorted_qgram[char_idx]][unsorted_qgram[char_idx]];
          }
          if(score-threshold_arr[char_idx] < 0) break;
        }
        if(score-threshold >= 0)
        { 
          const ScoreQgramcodePair elem{score,unsorted_idx};
          size_t i=0;
          while(i < environment[sorted_idx].size())
          {
            if(environment[sorted_idx][i].score < score) break;
            i++;
          }/*
          if(i == environment[sorted_idx].size()) environment[sorted_idx].push_back(elem);
          else
          {
            const auto it = environment[sorted_idx].begin()+i;
            environment[sorted_idx].insert(it,elem);
          }*/
          environment[sorted_idx].push_back(elem);
          for(size_t tmp_idx = environment[sorted_idx].size()-1;tmp_idx > i; tmp_idx--)
          {
            environment[sorted_idx][tmp_idx] = environment[sorted_idx][tmp_idx-1];
          }
          environment[sorted_idx][i] = ScoreQgramcodePair(score,unsorted_idx);
        }
      }
    }
  };
/*
  template<const uint16_t arr_len,const uint16_t sorted_idx>
  std::array<uint16_t,arr_len> threshold_idx_arr(std::vector<std::pair<int8_t,uint16_t>> sorted_env) const
  {
    const auto max_score = sorted_env[0].first;
    const uint16_t env_size = sorted_env.size()

    std::array<uint16_t,arr_len> threshold_idx_arr{};
    int8_t curr_score = max_score;
    uint16_t threshold_idx_arr_idx = static_cast<uint16_t>(arr_len-1);

    for(uint16_t env_idx = 0; env_idx < sorted_env.size(); env_idx++)
    {
      if(sorted_env[env_idx].first != curr_score)
      {
        const uint16_t diff = curr_score - sorted_env[env_idx].first;
        for(uint16_t i = 0; i < diff; i++)
        {
          threshold_idx_arr[threshold_idx_arr_idx] = env_idx;
          threshold_idx_arr_idx--;
        }
        curr_score = sorted_env[env_idx].first;
      }
    }

    return threshold_idx_arr;
  }*/

  std::vector<ScoreQgramcodePair> qgram_env(const uint16_t sorted_idx) const
  {
    return environment[sorted_idx];
  }

  ScoreQgramcodePair env_get(const size_t sorted_idx,const size_t i) const
  {
    if(sorted_idx >= sorted_size or i >= environment[sorted_idx].size())
    {
      return ScoreQgramcodePair(0,0);
    }
    return environment[sorted_idx][i];
  }

  size_t sorted_size_get() const
  {
    return sorted_size;
  }

  size_t size_get(const size_t sorted_idx) const
  {
    if(sorted_idx >= sorted_size)
    {
      return 0;
    }
    return environment[sorted_idx].size();
  }
};

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

//Based on gsa materials, blast_slide.pdf
template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
struct EnvMatrix {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
  //static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = power_1(undefined_rank,qgram_length);

  std::array<int8_t,unsorted_size*sorted_size> matrix{};
  std::array<uint16_t,sorted_size> size_arr{};

  EnvMatrix()
  {
    std::cout << "Init" << std::endl;
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();
    constexpr const auto pow_arr = create_pow_arr<qgram_length,undefined_rank>();

    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      size_arr[sorted_idx] = unsorted_size;
      const std::array<uint8_t,qgram_length> sorted_qgram = sorted_q.qgram_get(sorted_idx);
      for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
      {
        const uint8_t cc = sorted_qgram[q_idx];
        for(uint16_t idx_1 = 0; idx_1 < pow_arr[q_idx]; idx_1++)
        {
          for(uint16_t idx_2 = 0; idx_2 < unsorted_size/pow_arr[q_idx]; idx_2++)
          {
            const uint16_t arr_idx = idx_1*unsorted_size/pow_arr[q_idx] + idx_2; 
            if(matrix[sorted_idx*unsorted_size+arr_idx] != INT8_MIN)
            {
              matrix[sorted_idx*unsorted_size+arr_idx] += sc.score_matrix[cc][idx_1 % undefined_rank];
              if(matrix[sorted_idx*unsorted_size+arr_idx] < (int8_t)threshold_arr[q_idx])
              {
                matrix[sorted_idx*unsorted_size+arr_idx] = INT8_MIN;
                size_arr[sorted_idx]--;
              }
            }
          }
        }
      }
      /*
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
      }*/
    }
    /*
    constexpr_for<0,sorted_size,1>([&] (auto sorted_idx)
    {
      constexpr const auto sorted_qgram = sorted_q.qgram_get(sorted_idx);
      constexpr_for<0,qgram_length,1>([&] (auto q_idx)
      {
        const auto cc = sorted_qgram[q_idx];
        constexpr_for<0,pow_arr[q_idx],1>([&] (auto idx_1)
        {
          constexpr_for<0,unsorted_size/pow_arr[q_idx],1>([&] (auto idx_2)
          {
            const uint16_t arr_idx = idx_1*unsorted_size/pow_arr[q_idx] + idx_2; 
            if(V[sorted_idx*unsorted_size+arr_idx])
            {
              Pos[sorted_idx*unsorted_size+arr_idx] += sc.score_matrix[cc][idx_1 % undefined_rank];
              if((int8_t) Pos[sorted_idx*unsorted_size+arr_idx] < (int8_t)threshold_arr[q_idx])
              {
                V[sorted_idx*unsorted_size+arr_idx] = false;
                size_arr[sorted_idx]--;
              }
            }
          });
        });
      });
    });*/
  }
};
/*
template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
class QgramEnvironment2 {
  static constexpr const EnvMatrix<ScoreClass,qgram_length,threshold> env{};
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto sorted_size = env.sorted_size;
  static constexpr const auto unsorted_size = env.unsorted_size;
  
  static constexpr const auto matrix = env.matrix;
  static constexpr const auto size_arr = env.size_arr;

  public:
  constexpr QgramEnvironment2(){};

  int8_t get_score(const uint32_t& sorted_idx, const uint32_t& unsorted_idx) const
  {
    return matrix[sorted_idx*unsorted_size+unsorted_idx];
  }

  constexpr std::array<int8_t,sorted_size*unsorted_size> get_matrix() const {
    return matrix;
  }

  constexpr uint16_t multiset_size() const
  {
    return sorted_size;
  }

  constexpr uint16_t unsorted() const
  {
    return unsorted_size;
  }

  //get size of prefiltered, unsorted environment
  template<const uint16_t sorted_idx>
  constexpr uint16_t env_size_get() const
  {
    if(sorted_idx >= sorted_size) return 0;
    return size_arr[sorted_idx];
  }

  //sort (with insertion sort) and return environment of sorted qgram with index sorted_idx
  template<const uint16_t sorted_idx>
  constexpr std::array<ScoreQgramcodePair,size_arr[sorted_idx]> sort() const
  {
    constexpr const uint16_t env_size = size_arr[sorted_idx];
    std::array<ScoreQgramcodePair,env_size> sorted_env{};

    uint16_t unsorted_env_idx = 0;
    
    for(uint16_t arr_idx = 0; arr_idx < env_size; arr_idx++)
    {
      while(matrix[sorted_idx*unsorted_size+unsorted_env_idx]==INT8_MIN) unsorted_env_idx++;
      
      const auto score = matrix[sorted_idx*unsorted_size+unsorted_env_idx];

      uint16_t position_to_insert = 0;
      while(sorted_env[position_to_insert].score - score > 0 and position_to_insert < arr_idx) position_to_insert++;
      
      for(uint16_t tmp_idx = arr_idx; tmp_idx > position_to_insert;tmp_idx--)
      {
        sorted_env[tmp_idx] = sorted_env[tmp_idx-1];
      }
      sorted_env[position_to_insert] = ScoreQgramcodePair(score,unsorted_env_idx);
      unsorted_env_idx++;
    }
    return sorted_env;
  }
};*/