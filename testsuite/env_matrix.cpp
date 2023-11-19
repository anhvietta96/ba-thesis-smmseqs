#include "alignment/blosum62.hpp"
#include "filter/env_matrix.hpp"
#include "filter/utils.hpp"
#include <iostream>

int main(){
  static constexpr const uint16_t sq_code = 71;
  static constexpr const uint8_t qgram_length = 3;
  static constexpr const Blosum62 sc{};
  static constexpr const size_t undefined_rank = sc.num_of_chars;

  const EnvMatrix2<Blosum62,qgram_length> mat{};
  const ScoreQgramcodePair2 * const env = mat.sorted_env_get(sq_code);
  for(size_t i = 0; i < constexpr_pow(undefined_rank,qgram_length); i++){
    std::cout << env[i].code << '\t' << (int)env[i].score << std::endl;
  }
}