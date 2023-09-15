#include "filter/qgram_environment.hpp"
#include "utilities/constexpr_for.hpp"
#include <xmmintrin.h>
#include <x86intrin.h>

#include "alignment/simd.hpp"

#include <iostream>
#include <tuple>

template<class ScoreClass>
struct ScoreStats {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  std::array<int8_t,undefined_rank> min_of_each_char{};
  std::array<int8_t,undefined_rank> max_of_each_char{};

  constexpr ScoreStats()
  {
    constexpr_for<0,undefined_rank,1>([&] (auto char_idx_1)
    {
      int8_t min_score = 0;
      int8_t max_score = 0;
      constexpr_for<0,undefined_rank,1>([&] (auto char_idx_2)
      {
        if(sc.score_matrix[char_idx_1][char_idx_2] > max_score) max_score = sc.score_matrix[char_idx_1][char_idx_2];
        if(sc.score_matrix[char_idx_1][char_idx_2] < min_score) min_score = sc.score_matrix[char_idx_1][char_idx_2];
      });
      min_of_each_char[char_idx_1] = min_score;
      max_of_each_char[char_idx_1] = max_score;
    });
  };

  template<const uint8_t q>
  constexpr int8_t min_score_of_qgram(const std::array<uint8_t,q> qgram) const
  {
    int8_t min = 0;
    for(uint8_t i = 0; i < q; i++)
    {
      if(qgram[i] < undefined_rank)
      {
        min += min_of_each_char[qgram[i]];
      }
    }
    return min;
  }

  template<const uint8_t q>
  constexpr int8_t max_score_of_qgram(const std::array<uint8_t,q> qgram) const
  {
    int8_t max = 0;
    for(uint8_t i = 0; i < q; i++)
    {
      if(qgram[i] < undefined_rank)
      {
        max += max_of_each_char[qgram[i]];
      }
    }
    return max;
  }
};

//template<class ScoreClass, const char* 4gram, const int8_t threshold>
template<class ScoreClass,const uint16_t qgram_1_code,const uint16_t qgram_2_code,const int8_t threshold>
class EnvConstructor4 {
  static constexpr const ScoreClass sc{};
  std::vector<std::tuple<int8_t,uint16_t,uint16_t>> env{};
  
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

  public:
  EnvConstructor4()
  {
    constexpr const auto char_spec = sc.characters;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q{};
/*
    char sorted_qgram[4];
    char permutation[4];

    const bool 4qgram_sorted = sorted_q.sort(4gram,sorted_qgram,permutation,4);

    const uint16_t qgram_1_code = sorted_q.encode(qgram_1);
    const uint16_t qgram_2_code = sorted_q.encode(qgram_2);
*/
    constexpr const uint8_t q = 2;

    constexpr const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    constexpr const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    constexpr const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<q>(qgram_1);
    constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<q>(qgram_1);
    //constexpr const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<q>(qgram_2);
    constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<q>(qgram_2);

    constexpr const int8_t q_threshold = (threshold - max_score_qgram_1 > min_score_qgram_2) 
    ? (threshold - max_score_qgram_1 > min_score_qgram_1 ? min_score_qgram_1 : threshold - max_score_qgram_1)
    : (min_score_qgram_2 > min_score_qgram_1 ? min_score_qgram_1 : min_score_qgram_2) ;
    

    constexpr const QgramEnvironment2<ScoreClass,2,q_threshold> q_env{};

    constexpr const auto env_size_1 = q_env.env_size(qgram_1_code);
    constexpr const auto env_size_2 = q_env.env_size(qgram_2_code);
    
    constexpr const auto env_qgram_1 = q_env.template sort<env_size_1>(qgram_1_code);
    constexpr const auto env_qgram_1_score = env_qgram_1.first;
    constexpr const auto env_qgram_1_code = env_qgram_1.second;
    
    constexpr const auto env_qgram_2 = q_env.template sort<env_size_2>(qgram_2_code);
    constexpr const auto env_qgram_2_score = env_qgram_2.first;
    constexpr const auto env_qgram_2_code = env_qgram_2.second;

    const auto env_1_begin_it = env_qgram_1.begin();
    const auto env_2_end_it = env_qgram_2.end();

    constexpr const uint16_t score_diff_1 = static_cast<uint16_t>(env_qgram_1_score[0] - env_qgram_1_score[env_qgram_1.size()-1]);
    constexpr const uint16_t score_diff_2 = static_cast<uint16_t>(env_qgram_2_score[0] - env_qgram_2_score[env_qgram_2.size()-1]);

    const auto threshold_idx_arr_1 = threshold_idx_arr<env_size_1,score_diff_1,env_qgram_1_score[0],qgram_1_code>(env_qgram_1_score);
    const auto threshold_idx_arr_2 = threshold_idx_arr<env_size_2,score_diff_2,env_qgram_1_score[0],qgram_2_code>(env_qgram_2_score);
/*
    for(int8_t score_1 = min_score_qgram_1; score_1 <= max_score_qgram_1; score_1++)
    {
      const int8_t score_2 = threshold - score_1;
      if(score_2 < min_score_qgram_2 or score_2 > max_score_qgram_2) continue;

      const auto it_2_pos = (score_2 != max_score_qgram_2) ? (threshold_idx_arr_2[score_2 - min_score_qgram_2]) : 0;

      auto it_2 = env_qgram_2.begin() + it_2_pos;

      auto it_1 = (score_1 == min_score_qgram_1) ? env_qgram_1.end() : (env_qgram_1.begin() + threshold_idx_arr_1[score_1 - min_score_qgram_1 - 1] - 1);

      for(;it_1 != env_1_begin_it; it_1--)
      {
        const auto elem_1 = *(it_1);
        for(;it_2 != env_2_end_it; it_2++)
        {
          const auto elem_2 = *(it_2);
          env.push_back(std::tuple<int8_t,uint16_t,uint16_t>(elem_1.first+elem_2.first,elem_1.second,elem_2.second));
        }
      }
    }

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1; score_1--)
    {
      const int8_t min_score_for_2 = threshold - score_1;
      const uint16_t min_score_2_pos = 
      const uint16_t score_1_pos = 
      const uint16_t pos_1 = ;
      const uint16_t pos_2 = ;

    }*/
  };

  size_t size() const
  {
    return env.size();
  }
  
  std::tuple<int8_t,uint16_t,uint16_t> get_elem(const size_t env_idx) const
  {
    return env[env_idx];
  }
};


template<class ScoreClass,const uint16_t qgram_1_code,const uint16_t qgram_2_code,const int8_t threshold>
class Env4Constructor_1 {
  static constexpr const ScoreClass sc{};
  std::vector<std::tuple<int8_t,uint16_t,uint16_t>> env{};

  std::vector<std::pair<uint16_t,uint16_t>> create_threshold(std::vector<std::pair<int8_t,uint16_t>> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].first;
    uint16_t start = 0;
    std::vector<std::pair<uint16_t,uint16_t>> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.first < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.first -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.first;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using ScorePosition = std::pair<uint16_t,uint16_t>;
  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  public:
  Env4Constructor_1()
  {
    constexpr const auto char_spec = sc.characters;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q{};

    constexpr const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    constexpr const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    constexpr const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<2>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    constexpr const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<2>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);
    std::cout << (int) max_score_qgram_1 << '\t' << (int) max_score_qgram_2 << std::endl;
    /*constexpr const int8_t q_threshold = (threshold - max_score_qgram_1 > min_score_qgram_2) 
    ? (threshold - max_score_qgram_1 > min_score_qgram_1 ? min_score_qgram_1 : threshold - max_score_qgram_1)
    : (min_score_qgram_2 > min_score_qgram_1 ? min_score_qgram_1 : min_score_qgram_2) ;*/

    constexpr const int8_t q_threshold = (max_score_qgram_1 < max_score_qgram_2) ? (threshold - max_score_qgram_2) :
      (threshold - max_score_qgram_1);

    const QgramEnvironment<ScoreClass,2,q_threshold> q_env{};

    const auto qgram_env_1 = q_env.qgram_env(qgram_1_code);
    const auto qgram_env_2 = q_env.qgram_env(qgram_2_code);

    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].first;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].first;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);

    using ScorePosition = std::pair<uint16_t,uint16_t>;
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      std::cout << (int) threshold_arr_1[pos_1].first << std::endl;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        std::cout << (int) threshold_arr_2[pos_2].first << std::endl;
        
        const auto elem = std::tuple<int8_t,ScorePosition,ScorePosition>(score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]);
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          env.push_back(std::tuple<int8_t,uint16_t,uint16_t>(score,qgram_env_1[pos_1].second,qgram_env_2[pos_2].second));        
        }
      }
    }
  };

  using ReducedEnvElemType = std::tuple<int8_t,uint16_t,uint16_t>;

  uint32_t size() const
  {
    return env.size();
  }

  ReducedEnvElemType elem_get(const uint32_t idx) const
  {
    return env[idx];
  }
};


template<class ScoreClass,const uint16_t qgram_1_code,const uint16_t qgram_2_code,const int8_t threshold>
class Env4Constructor_2 {
  static constexpr const ScoreClass sc{};
  std::vector<std::pair<int8_t,std::array<uint8_t,4>>> env{};

  std::vector<std::pair<uint16_t,uint16_t>> create_threshold(std::vector<std::pair<int8_t,uint16_t>> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].first;
    uint16_t start = 0;
    std::vector<std::pair<uint16_t,uint16_t>> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.first < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.first -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.first;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using ScorePosition = std::pair<uint16_t,uint16_t>;
  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  std::array<uint8_t,4> process_qgram(std::array<uint8_t,2> u_qgram_1,std::array<uint8_t,2> u_qgram_2,
  const uint8_t* permutation, const bool orig_qgram_sorted)
  {
    if(orig_qgram_sorted)
    {
      return std::array<uint8_t,4>{ u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1]};
    }
    else
    {
      const uint8_t ref_qgram[4] = { u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1]};
      std::array<uint8_t,4> transformed_qgram{};
      for(uint8_t idx = 0; idx < 4; idx++)
      {
        transformed_qgram[permutation[idx]] = ref_qgram[idx];
      }
      return transformed_qgram;
    }
  }

  public:
  Env4Constructor_2()
  {
    constexpr const auto char_spec = sc.characters;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,2> unsorted_q{};

    constexpr const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    constexpr const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    uint8_t sorted_qgram[4];
    uint8_t permutation[4];

    constexpr const uint8_t qgram[] = {qgram_1[0],qgram_1[1],qgram_2[0],qgram_2[1]};
    const bool orig_qgram_sorted = sorted_q.sort(qgram,sorted_qgram,permutation,4);

    constexpr const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<2>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    constexpr const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<2>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);
    //std::cout << (int) max_score_qgram_1 << '\t' << (int) max_score_qgram_2 << std::endl;
    /*constexpr const int8_t q_threshold = (threshold - max_score_qgram_1 > min_score_qgram_2) 
    ? (threshold - max_score_qgram_1 > min_score_qgram_1 ? min_score_qgram_1 : threshold - max_score_qgram_1)
    : (min_score_qgram_2 > min_score_qgram_1 ? min_score_qgram_1 : min_score_qgram_2) ;*/

    constexpr const int8_t q_threshold = (max_score_qgram_1 < max_score_qgram_2) ? (threshold - max_score_qgram_2) :
      (threshold - max_score_qgram_1);

    const QgramEnvironment<ScoreClass,2,q_threshold> q_env{};

    const auto qgram_env_1 = q_env.qgram_env(qgram_1_code);
    const auto qgram_env_2 = q_env.qgram_env(qgram_2_code);

    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].first;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].first;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);

    using ScorePosition = std::pair<uint16_t,uint16_t>;
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        
        const auto elem = std::tuple<int8_t,ScorePosition,ScorePosition>(score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]);
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          const auto u_qgram_1 = unsorted_q.qgram_get(qgram_env_1[pos_1].second);
          const auto u_qgram_2 = unsorted_q.qgram_get(qgram_env_2[pos_2].second);
          env.push_back(std::make_pair(score,process_qgram(u_qgram_1,u_qgram_2,permutation,orig_qgram_sorted)));
        }
      }
    }
  };

  using ReducedEnvElemType = std::pair<int8_t,std::array<uint8_t,4>>;

  uint32_t size() const
  {
    return env.size();
  }

  ReducedEnvElemType elem_get(const uint32_t idx) const
  {
    return env[idx];
  }
};

template<class ScoreClass,const uint16_t qgram_1_code,const uint16_t qgram_2_code,const int8_t threshold>
class Env6Constructor_2 {
  static constexpr const ScoreClass sc{};
  std::vector<std::pair<int8_t,std::array<uint8_t,6>>> env{};

  std::vector<std::pair<uint16_t,uint16_t>> create_threshold(std::vector<std::pair<int8_t,uint16_t>> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].first;
    uint16_t start = 0;
    std::vector<std::pair<uint16_t,uint16_t>> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.first < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.first -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.first;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using ScorePosition = std::pair<uint16_t,uint16_t>;
  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  std::array<uint8_t,6> process_qgram(std::array<uint8_t,3> u_qgram_1,std::array<uint8_t,3> u_qgram_2,
  const uint8_t* permutation, const bool orig_qgram_sorted)
  {
    if(orig_qgram_sorted)
    {
      return std::array<uint8_t,6>{ u_qgram_1[0], u_qgram_1[1],u_qgram_1[2],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
    }
    else
    {
      const uint8_t ref_qgram[] = { u_qgram_1[0], u_qgram_1[1],u_qgram_1[2],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
      std::array<uint8_t,6> transformed_qgram{};
      for(uint8_t idx = 0; idx < 6; idx++)
      {
        transformed_qgram[permutation[idx]] = ref_qgram[idx];
      }
      return transformed_qgram;
    }
  }

  public:
  Env6Constructor_2()
  {
    constexpr const auto char_spec = sc.characters;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,3> sorted_q{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,3> unsorted_q{};

    constexpr const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    constexpr const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    uint8_t sorted_qgram[6];
    uint8_t permutation[6];

    constexpr const uint8_t qgram[] = { qgram_1[0], qgram_1[1],qgram_1[2],qgram_2[0],qgram_2[1],qgram_1[2]};
    const bool orig_qgram_sorted = sorted_q.sort(qgram,sorted_qgram,permutation,6);

    constexpr const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<3>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    constexpr const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<3>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);
    //std::cout << (int) max_score_qgram_1 << '\t' << (int) max_score_qgram_2 << std::endl;
    /*constexpr const int8_t q_threshold = (threshold - max_score_qgram_1 > min_score_qgram_2) 
    ? (threshold - max_score_qgram_1 > min_score_qgram_1 ? min_score_qgram_1 : threshold - max_score_qgram_1)
    : (min_score_qgram_2 > min_score_qgram_1 ? min_score_qgram_1 : min_score_qgram_2) ;*/

    constexpr const int8_t q_threshold = (max_score_qgram_1 < max_score_qgram_2) ? (threshold - max_score_qgram_2) :
      (threshold - max_score_qgram_1);

    const QgramEnvironment<ScoreClass,3,q_threshold> q_env{};

    const auto qgram_env_1 = q_env.qgram_env(qgram_1_code);
    const auto qgram_env_2 = q_env.qgram_env(qgram_2_code);

    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].first;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].first;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);

    using ScorePosition = std::pair<uint16_t,uint16_t>;
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        
        const auto elem = std::tuple<int8_t,ScorePosition,ScorePosition>(score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]);
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          const auto u_qgram_1 = unsorted_q.qgram_get(qgram_env_1[pos_1].second);
          const auto u_qgram_2 = unsorted_q.qgram_get(qgram_env_2[pos_2].second);
          env.push_back(std::make_pair(score,process_qgram(u_qgram_1,u_qgram_2,permutation,orig_qgram_sorted)));
        }
      }
    }
  };

  using ReducedEnvElemType = std::pair<int8_t,std::array<uint8_t,6>>;

  uint32_t size() const
  {
    return env.size();
  }

  ReducedEnvElemType elem_get(const uint32_t idx) const
  {
    return env[idx];
  }
};

template<const uint8_t qgram_length>
struct LocalEnvElem {
  size_t position;
  int8_t score;
  uint8_t qgram[qgram_length];

  LocalEnvElem(const size_t _position, const int8_t _score,
    const std::array<uint8_t,qgram_length> _qgram){
    position = _position;
    score = _score;
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = _qgram[i];
    }
  };

  void set_qgram(const uint8_t* ref_qgram)
  {
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = ref_qgram[i];
    }
  }
};

template<class ScoreClass,const int8_t threshold>
class Env4Constructor_3 {
  static constexpr const ScoreClass sc{};
  static constexpr const int8_t q_threshold = threshold - 2*sc.highest_score;
  const QgramEnvironment<ScoreClass,2,q_threshold> env_2{};

  std::vector<LocalEnvElem<4>> constructed_env{};

  using ScorePosition = std::pair<uint16_t,uint16_t>;
  std::vector<ScorePosition> create_threshold(std::vector<ScoreQgramcodePair> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].score;
    uint16_t start = 0;
    std::vector<ScorePosition> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.score < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.score -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.score;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  std::array<uint8_t,4> process_qgram(std::array<uint8_t,2> u_qgram_1,std::array<uint8_t,2> u_qgram_2,
  const std::array<uint8_t,4> permutation, const bool orig_qgram_sorted) const
  {
    if(orig_qgram_sorted)
    {
      return std::array<uint8_t,4>{ u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1]};
    }
    else
    {
      const uint8_t ref_qgram[4] = { u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1]};
      std::array<uint8_t,4> transformed_qgram{};
      for(uint8_t idx = 0; idx < 4; idx++)
      {
        transformed_qgram[permutation[idx]] = ref_qgram[idx];
      }
      return transformed_qgram;
    }
  }

  void process_qgram_with_simd(const std::array<uint8_t,4> permutation, const bool orig_qgram_sorted,
  const size_t position)
  {
    if(orig_qgram_sorted) return;
    size_t i = constructed_env.size()-1;
    
    while(i < constructed_env.size())
    {
      if(constructed_env[i].position != position) break;
      i--;
    }

    const size_t num_of_qgrams = constructed_env.size()-1-i;
    const size_t num_of_simd_operations = (num_of_qgrams) ? ((num_of_qgrams-1) / 4 + 1) : 0;
    size_t qgram_pos = i+1;
    size_t curr_qgram = 0;
    //std::cout << (int)position << '\t'<<(int) i << '\t' << (int) num_of_qgrams << '\t' << (int) num_of_simd_operations << '\t' << (int) qgram_pos << std::endl;
    
    //invert permutation
    std::array<uint8_t,4> inv_permutation{};
    for(uint8_t k = 0; k < 4; k++)
    {
      inv_permutation[permutation[k]] = k;
    }

    //create __m128i for transformation
    uint8_t transformation[16] = {
      static_cast<uint8_t>(inv_permutation[0]+0),static_cast<uint8_t>(inv_permutation[1]+0),static_cast<uint8_t>(inv_permutation[2]+0),static_cast<uint8_t>(inv_permutation[3]+0),
      static_cast<uint8_t>(inv_permutation[0]+4),static_cast<uint8_t>(inv_permutation[1]+4),static_cast<uint8_t>(inv_permutation[2]+4),static_cast<uint8_t>(inv_permutation[3]+4),
      static_cast<uint8_t>(inv_permutation[0]+8),static_cast<uint8_t>(inv_permutation[1]+8),static_cast<uint8_t>(inv_permutation[2]+8),static_cast<uint8_t>(inv_permutation[3]+8),
      static_cast<uint8_t>(inv_permutation[0]+12),static_cast<uint8_t>(inv_permutation[1]+12),static_cast<uint8_t>(inv_permutation[2]+12),static_cast<uint8_t>(inv_permutation[3]+12),
    };

    __m128i *transform_m128i = (__m128i *) transformation;
    __m128i loaded_transformation = _mm_load_si128(transform_m128i);
    
    for (size_t simd_op_no = 0; simd_op_no < num_of_simd_operations; simd_op_no++)
    {
      uint8_t buf[16]{};
      uint8_t filled_qgram = 0;
      while(curr_qgram < num_of_qgrams and filled_qgram < 4)
      {
        //fill buffer
        auto qgram = constructed_env[qgram_pos].qgram;
        for(uint8_t char_idx = 0; char_idx < 4; char_idx++)
        {
          buf[filled_qgram*4+char_idx] = qgram[char_idx];
        }
        //std::cout << (int)qgram_pos << std::endl;
        filled_qgram++;
        curr_qgram++;
        qgram_pos++;
      }

      __m128i *vector = (__m128i *) buf;
      __m128i loaded_vector = _mm_load_si128(vector);

      //apply on buffer
      __m128i loaded_result = _mm_shuffle_epi8(loaded_vector,loaded_transformation);
      _mm_store_si128(vector,loaded_result);

      uint8_t* result = (uint8_t*) vector;
      //std::cout << (int) filled_qgram << std::endl;
      //save back in env
      for(uint8_t saved_qgram = filled_qgram; saved_qgram > 0; saved_qgram--)
      {
        uint8_t curr_saved_qgram = qgram_pos - saved_qgram;
        //std::cout << '2' << '\t' << (int) curr_saved_qgram << std::endl;
        assert(curr_saved_qgram < constructed_env.size() and curr_saved_qgram >= 0);

        uint8_t extracted_from_result[4]{};
        for(uint8_t j = 0; j < 4; j++)
        {
          extracted_from_result[j] = result[(filled_qgram-saved_qgram)*4+j];
        }
        
        constructed_env[curr_saved_qgram].set_qgram(extracted_from_result);
        /*
        auto qgram = constructed_env[curr_saved_qgram].qgram;
        for(uint8_t k = 0; k < 4; k++)
        {
          std::cout << (int)qgram[k] << '\t'; 
        }
        std::cout << std::endl;*/
      }
    }
  }

  public:
  Env4Constructor_3(){};

  //using ReducedEnvElemType = std::pair<int8_t,std::array<uint8_t,4>>;
  void env_construct(const size_t unsorted_qgram_code,const std::array<uint8_t,4> &permutation, const bool sorted,
                      const size_t position, const bool with_simd)
  {
    constexpr const auto char_spec = sc.character_spec;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,2> unsorted_q{};

    const uint16_t unsorted_qgram_1_code = static_cast<uint16_t>(unsorted_qgram_code / (undefined_rank*undefined_rank));
    const uint16_t unsorted_qgram_2_code = static_cast<uint16_t>(unsorted_qgram_code % (undefined_rank*undefined_rank));

    const uint16_t qgram_1_code = sorted_q.sorted_code_get(unsorted_qgram_1_code);
    const uint16_t qgram_2_code = sorted_q.sorted_code_get(unsorted_qgram_2_code);

    //std::cout << "1" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
/*
    for(size_t i = 0; i < undefined_rank*undefined_rank; i++)
    {
      std::cout << (int) i << '\t' << (int) sorted_q.sorted_code_get(i) << std::endl;
    }
*/  
    const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<2>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<2>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);

    const auto qgram_env_1 = env_2.qgram_env(qgram_1_code);
    const auto qgram_env_2 = env_2.qgram_env(qgram_2_code);
    
    /*if(position == 0)
    {
      //std::cout << "1" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
      for(size_t i = 0; i < qgram_env_1.size(); i++)
      {
        const auto elem = qgram_env_1[i];
        std::cout << (int)elem.score << '\t' << (int)elem.code << std::endl;
      }
      for(size_t i = 0; i < qgram_env_2.size(); i++)
      {
        const auto elem = qgram_env_2[i];
        std::cout << (int)elem.score << '\t' << (int)elem.code << std::endl;
      }
    }*/
    
    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].score;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].score;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);
    //std::cout << "2" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        
        const std::tuple<int8_t,ScorePosition,ScorePosition> elem{score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]};
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    //std::cout << "3" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          const auto u_qgram_1 = unsorted_q.qgram_get(qgram_env_1[pos_1].code);
          const auto u_qgram_2 = unsorted_q.qgram_get(qgram_env_2[pos_2].code);
          if(!with_simd)
          {
            auto rearranged_qgram = process_qgram(u_qgram_1,u_qgram_2,permutation,sorted);
            LocalEnvElem<4> env_elem{position,score,rearranged_qgram};
            constructed_env.push_back(env_elem);
          }
          else
          {
            std::array<uint8_t,4> unrearranged_qgram = {u_qgram_1[0],u_qgram_1[1],u_qgram_2[0],u_qgram_2[1]};
            LocalEnvElem<4> env_elem{position,score,unrearranged_qgram};
            constructed_env.push_back(env_elem);
          }
        }
      }
    }

    if(with_simd)
    {
      process_qgram_with_simd(permutation,sorted,position);
    }
  }

  size_t size() const
  {
    return constructed_env.size();
  }

  LocalEnvElem<4> elem_get(const size_t idx) const
  {
    return constructed_env[idx];
  }
};

template<class ScoreClass,const int8_t threshold>
class Env6Constructor_3 {
  static constexpr const ScoreClass sc{};
  static constexpr const uint8_t qgram_length = 3;
  static constexpr const int8_t q_threshold = threshold - qgram_length*sc.highest_score;
  const QgramEnvironment<ScoreClass,qgram_length,q_threshold> env_3{};

  std::vector<LocalEnvElem<6>> constructed_env{};
  
  using ScorePosition = std::pair<uint16_t,uint16_t>;
  std::vector<ScorePosition> create_threshold(std::vector<ScoreQgramcodePair> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].score;
    uint16_t start = 0;
    std::vector<ScorePosition> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.score < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.score -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.score;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  std::array<uint8_t,6> process_qgram(std::array<uint8_t,qgram_length> u_qgram_1,std::array<uint8_t,qgram_length> u_qgram_2,
  const std::array<uint8_t,6> permutation, const bool orig_qgram_sorted) const
  {
    if(orig_qgram_sorted)
    {
      return std::array<uint8_t,6>{ u_qgram_1[0], u_qgram_1[1],u_qgram_1[2],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
    }
    else
    {
      const uint8_t ref_qgram[6] = { u_qgram_1[0], u_qgram_1[1],u_qgram_1[2],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
      std::array<uint8_t,6> transformed_qgram{};
      for(uint8_t idx = 0; idx < 6; idx++)
      {
        transformed_qgram[permutation[idx]] = ref_qgram[idx];
      }
      return transformed_qgram;
    }
  }

  void process_qgram_with_simd(const std::array<uint8_t,6> permutation, const bool orig_qgram_sorted,
  const size_t position)
  {
    if(orig_qgram_sorted) return;
    size_t i = constructed_env.size()-1;
    
    while(i < constructed_env.size())
    {
      if(constructed_env[i].position != position) break;
      i--;
    }

    const size_t num_of_qgrams = constructed_env.size()-1-i;
    const size_t num_of_simd_operations = (num_of_qgrams) ? ((num_of_qgrams-1) / 2 + 1) : 0;
    size_t qgram_pos = i+1;
    size_t curr_qgram = 0;
    //std::cout << (int)position << '\t'<<(int) i << '\t' << (int) num_of_qgrams << '\t' << (int) num_of_simd_operations << '\t' << (int) qgram_pos << std::endl;
    
    //invert permutation
    std::array<uint8_t,6> inv_permutation{};
    for(uint8_t k = 0; k < 6; k++)
    {
      inv_permutation[permutation[k]] = k;
    }

    //create __m128i for transformation
    uint8_t transformation[16] = {
      static_cast<uint8_t>(inv_permutation[0]+0),static_cast<uint8_t>(inv_permutation[1]+0),static_cast<uint8_t>(inv_permutation[2]+0),static_cast<uint8_t>(inv_permutation[3]+0),
      static_cast<uint8_t>(inv_permutation[4]+0),static_cast<uint8_t>(inv_permutation[5]+0),static_cast<uint8_t>(inv_permutation[0]+6),static_cast<uint8_t>(inv_permutation[1]+6),
      static_cast<uint8_t>(inv_permutation[2]+6),static_cast<uint8_t>(inv_permutation[3]+6),static_cast<uint8_t>(inv_permutation[4]+6),static_cast<uint8_t>(inv_permutation[5]+6),
      0x80,0x80,0x80,0x80
    };/*
    std::cout << "Permutation\t";
    for(uint8_t j = 0; j < 16; j++)
    {
      std::cout << (int)transformation[j] << '\t';
    }
    std::cout << std::endl;*/

    __m128i *transform_m128i = (__m128i *) transformation;
    __m128i loaded_transformation = _mm_load_si128(transform_m128i);
    
    for (size_t simd_op_no = 0; simd_op_no < num_of_simd_operations; simd_op_no++)
    {
      uint8_t buf[16]{};
      uint8_t filled_qgram = 0;
      while(curr_qgram < num_of_qgrams and filled_qgram < 2)
      {
        //fill buffer
        auto qgram = constructed_env[qgram_pos].qgram;
        for(uint8_t char_idx = 0; char_idx < 6; char_idx++)
        {
          buf[filled_qgram*6+char_idx] = qgram[char_idx];
        }
        //std::cout << (int)qgram_pos << std::endl;
        filled_qgram++;
        curr_qgram++;
        qgram_pos++;
      }/*
      std::cout << "Buffer\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)buf[j] << '\t';
      }
      std::cout << std::endl;*/

      __m128i *vector = (__m128i *) buf;
      __m128i loaded_vector = _mm_load_si128(vector);

      //apply on buffer
      __m128i loaded_result = _mm_shuffle_epi8(loaded_vector,loaded_transformation);
      _mm_store_si128(vector,loaded_result);

      uint8_t* result = (uint8_t*) vector;
      //std::cout << (int) filled_qgram << std::endl;
      //save back in env
      for(uint8_t saved_qgram = filled_qgram; saved_qgram > 0; saved_qgram--)
      {
        uint8_t curr_saved_qgram = qgram_pos - saved_qgram;
        //std::cout << '2' << '\t' << (int) curr_saved_qgram << std::endl;
        assert(curr_saved_qgram < constructed_env.size() and curr_saved_qgram >= 0);

        uint8_t extracted_from_result[6]{};
        for(uint8_t j = 0; j < 6; j++)
        {
          extracted_from_result[j] = result[(filled_qgram-saved_qgram)*6+j];
        }
        
        constructed_env[curr_saved_qgram].set_qgram(extracted_from_result);
        /*
        auto qgram = constructed_env[curr_saved_qgram].qgram;
        for(uint8_t k = 0; k < 6; k++)
        {
          std::cout << (int)qgram[k] << '\t'; 
        }
        std::cout << std::endl;*/
      }/*
      std::cout << "Result\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)result[j] << '\t';
      }
      std::cout << std::endl;*/
    }
  }

  public:
  Env6Constructor_3(){};

  //using ReducedEnvElemType = std::pair<int8_t,std::array<uint8_t,4>>;
  void env_construct(const size_t unsorted_qgram_code,const std::array<uint8_t,6> permutation, const bool sorted,
    const size_t position, const bool with_simd)
  {
    constexpr const auto char_spec = sc.character_spec;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};
    constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};

    const uint16_t unsorted_qgram_1_code = static_cast<uint16_t>(unsorted_qgram_code / (undefined_rank*undefined_rank*undefined_rank));
    const uint16_t unsorted_qgram_2_code = static_cast<uint16_t>(unsorted_qgram_code % (undefined_rank*undefined_rank*undefined_rank));

    const uint16_t qgram_1_code = sorted_q.sorted_code_get(unsorted_qgram_1_code);
    const uint16_t qgram_2_code = sorted_q.sorted_code_get(unsorted_qgram_2_code);

    //std::cout << "1" << '\t' << (int) unsorted_qgram_1_code << '\t' << (int) unsorted_qgram_2_code << std::endl;
/*
    for(size_t i = 0; i < undefined_rank*undefined_rank; i++)
    {
      std::cout << (int) i << '\t' << (int) sorted_q.sorted_code_get(i) << std::endl;
    }
*/
    const auto qgram_1 = sorted_q.qgram_get(qgram_1_code);
    const auto qgram_2 = sorted_q.qgram_get(qgram_2_code);

    const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<qgram_length>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<qgram_length>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);

    const auto qgram_env_1 = env_3.qgram_env(qgram_1_code);
    const auto qgram_env_2 = env_3.qgram_env(qgram_2_code);

    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].score;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].score;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);
    //std::cout << "2" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        
        const auto elem = std::tuple<int8_t,ScorePosition,ScorePosition>(score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]);
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    //std::cout << "3" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          const auto u_qgram_1 = unsorted_q.qgram_get(qgram_env_1[pos_1].code);
          const auto u_qgram_2 = unsorted_q.qgram_get(qgram_env_2[pos_2].code);
          //constructed_env.push_back(std::make_pair(score,process_qgram(u_qgram_1,u_qgram_2,permutation,sorted)));
          if(!with_simd)
          {
            const auto rearranged_qgram = process_qgram(u_qgram_1,u_qgram_2,permutation,sorted);
            const LocalEnvElem<6> env_elem{position,score,rearranged_qgram};
            constructed_env.push_back(env_elem);
          }
          else
          {
            const std::array<uint8_t,6> qgram = {u_qgram_1[0],u_qgram_1[1],u_qgram_1[2],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
            const LocalEnvElem<6> env_elem{position,score,qgram};
            constructed_env.push_back(env_elem);
          }
        }
      }
    }
    
    if(with_simd)
    {
      process_qgram_with_simd(permutation,sorted,position);
    }
  }

  size_t size() const
  {
    return constructed_env.size();
  }

  LocalEnvElem<6> elem_get(const size_t idx) const
  {
    return constructed_env[idx];
  }
};

template<class ScoreClass,const int8_t threshold>
class Env5Constructor {
  static constexpr const ScoreClass sc{};
  static constexpr const int8_t q_threshold_1 = threshold - 3*sc.highest_score;
  static constexpr const int8_t q_threshold_2 = threshold - 2*sc.highest_score;
  
  const QgramEnvironment<ScoreClass,2,q_threshold_1> env_2{};
  const QgramEnvironment<ScoreClass,3,q_threshold_2> env_3{};

  std::vector<LocalEnvElem<5>> constructed_env{};
  
  using ScorePosition = std::pair<uint16_t,uint16_t>;
  std::vector<ScorePosition> create_threshold(std::vector<ScoreQgramcodePair> qgram_env) const
  {
    int8_t curr_score = qgram_env[0].score;
    uint16_t start = 0;
    std::vector<ScorePosition> threshold_arr{};
    
    for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
    {
      const auto elem = qgram_env[idx];
      if(elem.score < curr_score)
      {
        threshold_arr.push_back(std::make_pair(start,idx-1));
        for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.score -1); j++)
        {
          threshold_arr.push_back(std::make_pair(__UINT16_MAX__,__UINT16_MAX__));  
        }
        start = idx;
        curr_score = elem.score;
      }
    }
    threshold_arr.push_back(std::make_pair(start,qgram_env.size()-1));

    return threshold_arr;
  }

  using EnvElemType = std::tuple<int8_t,ScorePosition,ScorePosition>;
  void insert_sort(std::vector<EnvElemType> temp_env, EnvElemType elem)
  {
    if(temp_env.size() == 0)
    {
      temp_env.push_back(elem);
    }
    else
    {
      auto it = temp_env.begin();
      while(std::get<0>(*it) >= std::get<0>(elem))
      {
        it++;
      }
      temp_env.insert(it,elem);
    }
  }

  std::array<uint8_t,5> process_qgram(std::array<uint8_t,2> u_qgram_1,std::array<uint8_t,3> u_qgram_2,
  const std::array<uint8_t,5> permutation, const bool orig_qgram_sorted) const
  {
    if(orig_qgram_sorted)
    {
      return std::array<uint8_t,5>{ u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
    }
    else
    {
      const uint8_t ref_qgram[5] = { u_qgram_1[0], u_qgram_1[1],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
      std::array<uint8_t,5> transformed_qgram{};
      for(uint8_t idx = 0; idx < 5; idx++)
      {
        transformed_qgram[permutation[idx]] = ref_qgram[idx];
      }
      return transformed_qgram;
    }
  }

  void process_qgram_with_simd(const std::array<uint8_t,5> permutation, const bool orig_qgram_sorted,
  const size_t position)
  {
    if(orig_qgram_sorted) return;
    size_t i = constructed_env.size()-1;
    
    while(i < constructed_env.size())
    {
      if(constructed_env[i].position != position) break;
      i--;
    }

    const size_t num_of_qgrams = constructed_env.size()-1-i;
    if(num_of_qgrams == 0) return;
    const size_t num_of_simd_operations = (num_of_qgrams) ? ((num_of_qgrams-1) / 3 + 1) : 0;
    size_t qgram_pos = i+1;
    size_t curr_qgram = 0;
    //std::cout << (int)position << '\t'<<(int) i << '\t' << (int) num_of_qgrams << '\t' << (int) num_of_simd_operations << '\t' << (int) qgram_pos << std::endl;
    
    //invert permutation
    std::array<uint8_t,5> inv_permutation{};
    for(uint8_t k = 0; k < 5; k++)
    {
      inv_permutation[permutation[k]] = k;
    }

    //create __m128i for transformation
    uint8_t transformation[16] = {
      static_cast<uint8_t>(inv_permutation[0]+0),static_cast<uint8_t>(inv_permutation[1]+0),static_cast<uint8_t>(inv_permutation[2]+0),static_cast<uint8_t>(inv_permutation[3]+0),
      static_cast<uint8_t>(inv_permutation[4]+0),static_cast<uint8_t>(inv_permutation[0]+5),static_cast<uint8_t>(inv_permutation[1]+5),static_cast<uint8_t>(inv_permutation[2]+5),
      static_cast<uint8_t>(inv_permutation[3]+5),static_cast<uint8_t>(inv_permutation[4]+5),static_cast<uint8_t>(inv_permutation[0]+10),static_cast<uint8_t>(inv_permutation[1]+10),
      static_cast<uint8_t>(inv_permutation[2]+10),static_cast<uint8_t>(inv_permutation[3]+10),static_cast<uint8_t>(inv_permutation[4]+10),0x80
    };/*
    std::cout << "Permutation\t";
    for(uint8_t j = 0; j < 16; j++)
    {
      std::cout << (int)transformation[j] << '\t';
    }
    std::cout << std::endl;*/

    __m128i *transform_m128i = (__m128i *) transformation;
    __m128i loaded_transformation = _mm_load_si128(transform_m128i);
    
    for (size_t simd_op_no = 0; simd_op_no < num_of_simd_operations; simd_op_no++)
    {
      uint8_t buf[16]{};
      uint8_t filled_qgram = 0;
      while(curr_qgram < num_of_qgrams and filled_qgram < 3)
      {
        //fill buffer
        auto qgram = constructed_env[qgram_pos].qgram;
        for(uint8_t char_idx = 0; char_idx < 5; char_idx++)
        {
          buf[filled_qgram*5+char_idx] = qgram[char_idx];
        }
        //std::cout << (int)qgram_pos << std::endl;
        filled_qgram++;
        curr_qgram++;
        qgram_pos++;
      }/*
      std::cout << "Buffer\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)buf[j] << '\t';
      }
      std::cout << std::endl;*/

      __m128i *vector = (__m128i *) buf;
      __m128i loaded_vector = _mm_load_si128(vector);

      //apply on buffer
      __m128i loaded_result = _mm_shuffle_epi8(loaded_vector,loaded_transformation);
      _mm_store_si128(vector,loaded_result);

      uint8_t* result = (uint8_t*) vector;
      //std::cout << (int) filled_qgram << std::endl;
      //save back in env
      for(uint8_t saved_qgram = filled_qgram; saved_qgram > 0; saved_qgram--)
      {
        uint8_t curr_saved_qgram = qgram_pos - saved_qgram;
        //std::cout << '2' << '\t' << (int) curr_saved_qgram << std::endl;
        assert(curr_saved_qgram < constructed_env.size() and curr_saved_qgram >= 0);

        uint8_t extracted_from_result[5]{};
        for(uint8_t j = 0; j < 5; j++)
        {
          extracted_from_result[j] = result[(filled_qgram-saved_qgram)*5+j];
        }
        
        constructed_env[curr_saved_qgram].set_qgram(extracted_from_result);
      }/*
      std::cout << "Result\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)result[j] << '\t';
      }
      std::cout << std::endl;*/
    }
  }

  public:
  Env5Constructor(){};

  void env_construct(const size_t unsorted_qgram_code,const std::array<uint8_t,5> permutation, const bool sorted,
    const size_t position, const bool with_simd)
  {
    constexpr const auto char_spec = sc.character_spec;
    constexpr const auto undefined_rank = sc.num_of_chars;
    constexpr ScoreStats<ScoreClass> stats{};

    constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q_2{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,2> unsorted_q_2{};
    constexpr const SortedQmer<char_spec,undefined_rank,3> sorted_q_3{};
    constexpr const UnsortedQmer<char_spec,undefined_rank,3> unsorted_q_3{};

    const uint16_t unsorted_qgram_1_code = static_cast<uint16_t>(unsorted_qgram_code / (undefined_rank*undefined_rank*undefined_rank));
    const uint16_t unsorted_qgram_2_code = static_cast<uint16_t>(unsorted_qgram_code % (undefined_rank*undefined_rank*undefined_rank));

    const uint16_t qgram_1_code = sorted_q_2.sorted_code_get(unsorted_qgram_1_code);
    const uint16_t qgram_2_code = sorted_q_3.sorted_code_get(unsorted_qgram_2_code);

    //std::cout << "1" << '\t' << (int) unsorted_qgram_1_code << '\t' << (int) unsorted_qgram_2_code << std::endl;
    //std::cout << "1" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
/*
    for(size_t i = 0; i < undefined_rank*undefined_rank; i++)
    {
      std::cout << (int) i << '\t' << (int) sorted_q.sorted_code_get(i) << std::endl;
    }
*/
    const auto qgram_1 = sorted_q_2.qgram_get(qgram_1_code);
    const auto qgram_2 = sorted_q_3.qgram_get(qgram_2_code);

    const int8_t max_score_qgram_1 = stats.template max_score_of_qgram<2>(qgram_1);
    //constexpr const int8_t min_score_qgram_1 = stats.template min_score_of_qgram<2>(qgram_1);
    const int8_t max_score_qgram_2 = stats.template max_score_of_qgram<3>(qgram_2);
    //constexpr const int8_t min_score_qgram_2 = stats.template min_score_of_qgram<2>(qgram_2);

    const auto qgram_env_1 = env_2.qgram_env(qgram_1_code);
    const auto qgram_env_2 = env_3.qgram_env(qgram_2_code);

    const int8_t min_score_qgram_1_in_env = qgram_env_1[qgram_env_1.size()-1].score;
    const int8_t min_score_qgram_2_in_env = qgram_env_2[qgram_env_2.size()-1].score;
    
    const auto threshold_arr_1 = create_threshold(qgram_env_1);
    const auto threshold_arr_2 = create_threshold(qgram_env_2);
    //std::cout << "2" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    
    std::vector<std::tuple<int8_t,ScorePosition,ScorePosition>> temp_env{};

    for(int8_t score_1 = max_score_qgram_1; score_1 >= min_score_qgram_1_in_env; score_1--)
    {
      const uint16_t pos_1 = static_cast<uint16_t>(max_score_qgram_1 - score_1);
      if(threshold_arr_1[pos_1].first == __UINT16_MAX__) continue;
      for(int8_t score_2 = max_score_qgram_2; score_2 >= threshold - score_1; score_2--)
      {
        if(score_2 < min_score_qgram_2_in_env) continue;
        const uint16_t pos_2 = static_cast<uint16_t>(max_score_qgram_2 - score_2);
        if(threshold_arr_2[pos_2].first == __UINT16_MAX__) continue;
        
        const auto elem = std::tuple<int8_t,ScorePosition,ScorePosition>(score_1+score_2,
          threshold_arr_1[pos_1],threshold_arr_2[pos_2]);
        if(temp_env.size() == 0)
        {
          temp_env.push_back(elem);
        }
        else
        {
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(std::get<0>(temp_env[idx]) < score_1+score_2)
            {
              break;
            }
            idx++;
          }
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }
      }
    }
    
    //std::cout << "3" << '\t' << (int) qgram_1_code << '\t' << (int) qgram_2_code << std::endl;
    for(auto const& elem : temp_env)
    {
      const auto score = std::get<0>(elem);
      const auto start_1 = std::get<1>(elem).first;
      const auto end_1 = std::get<1>(elem).second;
      const auto start_2 = std::get<2>(elem).first;
      const auto end_2 = std::get<2>(elem).second;
      for(uint16_t pos_1 = start_1; pos_1 <= end_1; pos_1++)
      {
        for(uint16_t pos_2 = start_2; pos_2 <= end_2; pos_2++)
        {
          const auto u_qgram_1 = unsorted_q_2.qgram_get(qgram_env_1[pos_1].code);
          const auto u_qgram_2 = unsorted_q_3.qgram_get(qgram_env_2[pos_2].code);
          //constructed_env.push_back(std::make_pair(score,process_qgram(u_qgram_1,u_qgram_2,permutation,sorted)));
          if(!with_simd)
          {
            const auto rearranged_qgram = process_qgram(u_qgram_1,u_qgram_2,permutation,sorted);
            const LocalEnvElem<5> env_elem{position,score,rearranged_qgram};
            constructed_env.push_back(env_elem);
          }
          else
          {
            const std::array<uint8_t,5> qgram = {u_qgram_1[0],u_qgram_1[1],u_qgram_2[0],u_qgram_2[1],u_qgram_2[2]};
            /*for(uint8_t i = 0; i < qgram.size(); i++)
            {
              std::cout << (int)qgram[i] << '\t';
            }
            std::cout << std::endl;*/
            const LocalEnvElem<5> env_elem{position,score,qgram};
            constructed_env.push_back(env_elem);
          }
        }
      }
    }
    
    if(with_simd)
    {
      process_qgram_with_simd(permutation,sorted,position);
    }
  }

  size_t size() const
  {
    return constructed_env.size();
  }

  LocalEnvElem<5> elem_get(const size_t idx) const
  {
    return constructed_env[idx];
  }
};