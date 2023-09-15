#include <iostream>
#include "filter/distribution.hpp"
#include "alignment/blosum62.hpp"

int main(){
  constexpr const uint8_t max_num_convolution = 6;
  constexpr const Distribution<Blosum62,max_num_convolution> distribution{};
  for(uint8_t i = 2; i <= max_num_convolution; i++){
    std::cout << (int)(i+1) << '\t' << (int) distribution.threshold_get(i) << std::endl;
  }
  //std::cout << (int) distribution.threshold(0.95) << std::endl;
}