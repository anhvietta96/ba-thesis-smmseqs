#include "filter/unsorted_qgram.hpp"
#include "filter/sorted_qgram.hpp"

#include <vector>

template<class ScoreClass,const size_t qgram_length, const int8_t threshold>
constexpr std::array<int8_t,qgram_length> create_threshold_arr()
{
  constexpr const ScoreClass sc{};
  constexpr const auto undefined_rank = sc.num_of_chars;
  std::array<int8_t,qgram_length> threshold_arr{};
  
  int8_t max_char_score = 0;
  constexpr_for<0,undefined_rank,1>([&] (auto idx1)
  {
    constexpr_for<0,undefined_rank,1>([&] (auto idx2)
    {
      if(sc.score_matrix[idx1][idx2] > max_char_score)
      {
        max_char_score = sc.score_matrix[idx1][idx2];
      }
    });
  });

  constexpr_for<0,qgram_length,1>([&] (auto char_idx)
  {
    threshold_arr[char_idx] = threshold - (qgram_length-char_idx) * max_char_score;
  });
  return threshold_arr;
}



template<class ScoreClass,const size_t qgram_length, const int8_t threshold>
class QgramEnvironment {
  private:
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = static_cast<size_t>(sc.num_of_chars);
  static constexpr const UnsortedQmer<sc.characters,undefined_rank,qgram_length> unsorted_q{};
  static constexpr const SortedQmer<sc.characters,undefined_rank,qgram_length> sorted_q{};
  static constexpr const size_t sorted_size = sorted_q.size_get();
  static constexpr const size_t unsorted_size = unsorted_q.size_get();

  //std::array<std::pair<int8_t,std::array<uint8_t,qgram_length>>,SIZE_AT_COMPILE_TIME> environment{};
  
  std::array<std::vector<std::pair<int8_t,std::array<uint8_t,qgram_length>>>,sorted_size> environment{};

  public:
  constexpr QgramEnvironment()
  {
    constexpr const auto threshold_arr = create_threshold_arr<ScoreClass,qgram_length,threshold>();

    for(size_t sorted_idx = 0; sorted_idx < sorted_size; sorted_idx++)
    {
      const std::array<uint8_t,qgram_length> sorted_qgram = sorted_q.qgram_get(sorted_idx);
      for(size_t unsorted_idx = 0; unsorted_idx < unsorted_size; unsorted_idx++)
      {
        const std::array<uint8_t,qgram_length> unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
        int8_t score = 0;
        for(size_t char_idx = 0; char_idx < qgram_length; char_idx++)
        {
          if(score < threshold_arr[char_idx]) break;
          if(sorted_qgram[char_idx] < undefined_rank and unsorted_qgram[char_idx] < undefined_rank)
          {
            score += sc.score_matrix[sorted_qgram[char_idx]][unsorted_qgram[char_idx]];
          }
        }
        if(score >= threshold)
        {
          auto elem = std::make_pair(score,unsorted_qgram);
          if(environment[sorted_idx].empty())
          {
            environment[sorted_idx].push_back(elem);
          }
          else
          {
            size_t i;
            for(i = 0; i < environment[sorted_idx].size() and environment[sorted_idx][i].first >= score; i++){}
            const auto it = environment[sorted_idx].begin()+i;
            environment[sorted_idx].insert(it,elem);
          }
        }
      }
    }
  };

  std::pair<int8_t,std::array<uint8_t,qgram_length>> env_get(const size_t sorted_idx,const size_t i) const
  {
    if(sorted_idx >= sorted_size or i >= environment[sorted_idx].size())
    {
      constexpr const std::array<uint8_t,qgram_length> err_q = {{ undefined_rank }};
      return std::make_pair(0,err_q);
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