#include <iostream>
#include "filter/multiset_code.hpp"
#include "alignment/blosum62.hpp"

int main(){
  constexpr const Blosum62 sc{};
  constexpr const auto undefined_rank = sc.num_of_chars;
  constexpr const auto char_spec = sc.character_spec;

  constexpr_for<2,4,1>([&] (auto qgram_length){
    const MultisetEncoder<char_spec,undefined_rank,qgram_length> encoder{};
    constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> sorted_q{};

    for(uint16_t sorted_idx = 0; sorted_idx < sorted_q.size_get(); sorted_idx++){
      const auto qgram = sorted_q.qgram_get(sorted_idx);
      const auto calculated_code = encoder.encode(qgram);
      assert(sorted_idx == calculated_code);
    }
  });
}