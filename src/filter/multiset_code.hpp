#include"filter/sorted_qgram.hpp"
#include "utilities/constexpr_for.hpp"
#include <iostream>

template<const char* char_spec,const size_t undefined_rank,const uint8_t qgram_length>
class MultisetEncoder {
  uint16_t multiset_weights[undefined_rank * qgram_length]{};

  public:
  MultisetEncoder(){
    for(size_t i = 0; i < undefined_rank; i++){
      multiset_weights[(qgram_length-1) * undefined_rank + i] = i;
    }

    constexpr_for<1,qgram_length,1>([&] (auto q_idx){
      size_t curr_char = undefined_rank;
      constexpr const SortedQmer<char_spec,undefined_rank,q_idx+1> sorted_q{};
      for(uint16_t sorted_idx = 0; sorted_idx < sorted_q.size_get(); sorted_idx++){
        const auto multiset = sorted_q.qgram_get(sorted_idx);
        if(multiset[0] == curr_char) continue;

        size_t suffixcode = 0;
        for(uint8_t pos = 0; pos < q_idx; pos++){
          suffixcode += multiset_weights[(qgram_length-1-pos)*undefined_rank+multiset[q_idx-pos]];
        }
        const uint16_t difference = sorted_idx - suffixcode;
        multiset_weights[(qgram_length-1-q_idx)*undefined_rank + multiset[0]] = difference;
        curr_char = multiset[0];
      }
    });
  };

  uint16_t encode(const std::array<uint8_t,qgram_length>& qgram) const {
    uint16_t code = 0;
    for(uint8_t q_idx = 0; q_idx < qgram_length; q_idx++){
      code += multiset_weights[(q_idx)*undefined_rank+qgram[q_idx]];
    }
    return code;
  }
};