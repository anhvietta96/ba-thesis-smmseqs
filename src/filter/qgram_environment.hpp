#include "filter/score_matrix.hpp"

template<class ScoreClass,const size_t qgram_length, const int8_t threshold>
class QgramEnvironment {
  private:
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = static_cast<size_t>(sc.num_of_chars);
  static constexpr const UnsortedQmer<sc.characters,undefined_rank,qgram_length> unsorted_q{};
  std::vector<std::pair<int8_t,std::array<uint8_t,qgram_length>>> environment = {};

  public:
  constexpr QgramEnvironment(const std::array<uint8_t,qgram_length> ref_qgram)
  {
    constexpr const size_t unsorted_size = unsorted_q.size_get();
    constexpr_for<0,unsorted_size,1>([&] (auto unsorted_idx)
    {
      const std::array<uint8_t,qgram_length> unsorted_qgram = unsorted_q.qgram_get(unsorted_idx);
      int8_t score = 0;
      constexpr_for<0,qgram_length,1>([&] (auto char_idx)
      {
        if(ref_qgram[char_idx] < undefined_rank and unsorted_qgram[char_idx] < undefined_rank)
        {
          score += sc.score_matrix[ref_qgram[char_idx]][unsorted_qgram[char_idx]];
        }
      });
      if(score >= threshold)
      {
        auto elem = std::make_pair(score,unsorted_qgram);
        environment.push_back(elem);
      }
    });
  };

  std::pair<int8_t,std::array<uint8_t,qgram_length>> env_get(size_t i) const
  {
    if(i >= environment.size())
    {
      constexpr const std::array<uint8_t,qgram_length> err_q = {{ undefined_rank }};
      return std::make_pair(0,err_q);
    }
    return environment[i];
  }

  size_t size() const
  {
    return environment.size();
  }
};