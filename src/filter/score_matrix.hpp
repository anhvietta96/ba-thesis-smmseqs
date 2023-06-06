#include "filter/sorted_q_mer.hpp"
#include "filter/unsorted_qgram.hpp"
#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <iostream>

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