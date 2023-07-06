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

  constexpr ScoreQgramcodePair(const int8_t _score, const uint16_t _code)
  : score { _score }, code { _code }{};

  constexpr ScoreQgramcodePair(const ScoreQgramcodePair& other_pair)
  : score { other_pair.score }, code { other_pair.code }{};

  constexpr ScoreQgramcodePair operator= (const ScoreQgramcodePair& other_pair)
  {
    return ScoreQgramcodePair(other_pair);
  };
};


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

  return threshold_arr;
}

bool CompareSCPair(ScoreQgramcodePair pair_1, ScoreQgramcodePair pair_2)
{
  return pair_1.score > pair_2.score;
}

template<class ScoreClass,const size_t qgram_length, const int8_t threshold>
class QgramEnvironment {
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = static_cast<size_t>(sc.num_of_chars);
  static constexpr const UnsortedQmer<sc.character_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const SortedQmer<sc.character_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const size_t sorted_size = sorted_q.size_get();
  static constexpr const size_t unsorted_size = unsorted_q.size_get();

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
            environment[sorted_idx][tmp_idx].score = environment[sorted_idx][tmp_idx-1].score;
            environment[sorted_idx][tmp_idx].code = environment[sorted_idx][tmp_idx-1].code;
          }
          environment[sorted_idx][i].score = score;
          environment[sorted_idx][i].code = unsorted_idx;
          /*if(score >= 11 and sorted_idx == 215)
          {
            std::cout << (int) i << '\t' << (int) score << '\t' << (int) unsorted_idx << '\t' << (int)environment[sorted_idx][i].score << '\t' << (int)environment[sorted_idx][i].code << std::endl;
          }*/
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

template<class ScoreClass,const uint8_t qgram_length, const int8_t threshold>
class QgramEnvironment2 {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;

  static constexpr const UnsortedQmer<sc.character_spec,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const SortedQmer<sc.character_spec,undefined_rank,qgram_length> sorted_q{};
  static constexpr const uint16_t sorted_size = sorted_q.size_get();
  static constexpr const uint16_t unsorted_size = unsorted_q.size_get();

  //Based on gsa materials, blast_slide.pdf
  std::array<bool,unsorted_size*sorted_size> V{};
  std::array<uint8_t,unsorted_size*sorted_size> Pos{};

  std::array<uint16_t,sorted_size> size_arr{};

  public:
  constexpr QgramEnvironment2()
  {
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();
    constexpr const auto pow_arr = create_pow_arr<qgram_length,undefined_rank>();
    /*
    constexpr_for<0,sorted_size,1>([&] (auto sorted_idx)
    {
      size_arr[sorted_idx] = unsorted_size;
      constexpr_for<0,unsorted_size,1>([&] (auto unsorted_idx)
      {
        V[sorted_idx*unsorted_size+unsorted_idx] = true;
        Pos[sorted_idx*unsorted_size+unsorted_idx] = 0;
      });
    });
    */
    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      size_arr[sorted_idx] = unsorted_size;
      for(uint16_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      {
        V[sorted_idx*unsorted_size+unsorted_idx] = true;
        Pos[sorted_idx*unsorted_size+unsorted_idx] = 0;
      }
    }
    
    
    for(uint16_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      const std::array<uint8_t,qgram_length> sorted_qgram = sorted_q.qgram_get(sorted_idx);
      for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++)
      {
        const uint8_t cc = sorted_qgram[q_idx];
        for(uint16_t idx_1 = 0; idx_1 < pow_arr[q_idx]; idx_1++)
        {
          for(uint16_t idx_2 = 0; idx_2 < unsorted_size/pow_arr[q_idx]; idx_2++)
          {
            const uint16_t arr_idx = idx_1*unsorted_size/pow_arr[q_idx] + idx_2; 
            if(V[sorted_idx*unsorted_size+arr_idx])
            {
              Pos[sorted_idx*unsorted_size+arr_idx] += sc.score_matrix[cc][idx_1 % undefined_rank];
              if((int8_t)Pos[sorted_idx*unsorted_size+arr_idx] < (int8_t)threshold_arr[q_idx])
              {
                V[sorted_idx*unsorted_size+arr_idx] = false;
                size_arr[sorted_idx]--;
              }
            }
          }
        }
      }
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
  };

  constexpr int8_t get_score(const uint16_t sorted_idx, const uint16_t unsorted_idx) const
  {
    return Pos[sorted_idx*unsorted_size+unsorted_idx];
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
  constexpr uint16_t env_size(const uint16_t sorted_idx) const
  {
    if(sorted_idx >= sorted_size) return 0;
    return size_arr[sorted_idx];
  }

  //sort (with insertion sort) and return environment of sorted qgram with index sorted_idx
  template<const uint16_t env_size>
  constexpr std::array<ScoreQgramcodePair,env_size> sort(const uint16_t sorted_idx) const
  //constexpr std::array<std::pair<int8_t,uint16_t>,env_size> sort(const uint16_t sorted_idx) const
  {
    std::array<ScoreQgramcodePair,env_size> sorted_env{};

    uint16_t unsorted_env_idx = 0;
    
    for(uint16_t arr_idx = 0; arr_idx < env_size; arr_idx++)
    {
      while(!V[sorted_idx*unsorted_size+unsorted_env_idx]) unsorted_env_idx++;
      
      uint16_t position_to_insert = 0;
      while(sorted_env[position_to_insert].score >= Pos[sorted_idx*unsorted_size+unsorted_env_idx] and position_to_insert < arr_idx) position_to_insert++;
      
      for(uint16_t tmp_idx = arr_idx; tmp_idx > position_to_insert;tmp_idx--)
      {
        sorted_env[tmp_idx] = sorted_env[tmp_idx-1];
      }
      sorted_env[position_to_insert].score = Pos[sorted_idx*unsorted_size+unsorted_env_idx];
      sorted_env[position_to_insert].code = unsorted_env_idx;
      unsorted_env_idx++;
    }/*
    constexpr_for<0,env_size,1>([&] (auto arr_idx)
    {
      while(!V[sorted_idx*unsorted_size+unsorted_env_idx]) unsorted_env_idx++;
      
      uint16_t position_to_insert = 0;
      while(sorted_env[position_to_insert].score >= Pos[sorted_idx*unsorted_size+unsorted_env_idx] and position_to_insert < arr_idx) position_to_insert++;
      
      for(uint16_t tmp_idx = arr_idx; tmp_idx > position_to_insert;tmp_idx--)
      {
        sorted_env[tmp_idx] = sorted_env[tmp_idx-1];
      }
      sorted_env[position_to_insert].score = Pos[sorted_idx*unsorted_size+unsorted_env_idx];
      sorted_env[position_to_insert].code = unsorted_env_idx;
      unsorted_env_idx++;
    });*/
    return sorted_env;
  }

  //create an array to index start position of score in environment
  //e.g array[pos]=i <=> env[j] <= pos - min_score for all j >= i, env[j] > pos - min_score for all j < i
  template<const uint16_t env_size,const uint16_t arr_len,const int8_t max_score,const uint16_t sorted_idx>
  constexpr std::array<uint16_t,arr_len> threshold_idx_arr(std::array<int8_t,env_size> sorted_env_score) const
  {
    std::array<uint16_t,arr_len> threshold_idx_arr{};
    int8_t curr_score = max_score;
    uint16_t threshold_idx_arr_idx = static_cast<uint16_t>(arr_len-2);

    for(uint16_t env_idx = 0; env_idx < env_size; env_idx++)
    {
      if(sorted_env_score[env_idx] != curr_score)
      {
        const uint16_t diff = curr_score - sorted_env_score[env_idx];
        for(uint16_t i = 0; i < diff-1; i++)
        {
          threshold_idx_arr[threshold_idx_arr_idx] = threshold_idx_arr[threshold_idx_arr_idx+1];
          threshold_idx_arr_idx--;
        }
        threshold_idx_arr[threshold_idx_arr_idx] = env_idx;
        threshold_idx_arr_idx--;
        curr_score = sorted_env_score[env_idx];
      }
    }

    return threshold_idx_arr;
  }
};