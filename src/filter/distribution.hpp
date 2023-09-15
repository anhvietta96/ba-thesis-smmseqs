#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include<array>
#include<cassert>
#include "filter/utils.hpp"

#ifndef THRESHOLD_FACTOR
#define THRESHOLD_FACTOR 4
#endif


template<class ScoreClass,const uint8_t max_num_of_convolution>
class Distribution {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const size_t num_possible_values_single = sc.highest_score - sc.smallest_score + 1;
  static constexpr const size_t num_possible_values_max = num_possible_values_single*max_num_of_convolution;
  
  std::array<uint8_t,num_possible_values_single> original_distribution{};
  std::array<uint32_t,num_possible_values_max> convoluted_distribution{};
  std::array<int8_t,max_num_of_convolution> threshold{};

  constexpr int8_t calc_approx_threshold(const uint8_t curr_iteration){
    const int8_t min_val = sc.smallest_score * (curr_iteration+1);
    //size_t arr_idx = num_possible_values_single*(curr_iteration+1) - 1;
    //size_t total = 0;
    std::array<size_t,num_possible_values_max> psum_arr{};
    
    uint16_t i = 0;
    size_t psum = 0;
    for(i = num_possible_values_max-1; i < num_possible_values_max; i--){
      psum += convoluted_distribution[i];
      psum_arr[i] = psum;
    }
    //const size_t cutoff = constexpr_pow(undefined_rank,(curr_iteration+1))
    const size_t cutoff = static_cast<size_t>(psum_arr[0] * THRESHOLD_FACTOR / 1000000);
    
    for(i = num_possible_values_max-1; i < num_possible_values_max; i--){
      if(psum_arr[i] > cutoff) break;
    }

    //std::cout << (int) static_cast<int8_t>(i)+1+min_val << std::endl;
    return i + 1 + min_val;

    /*
    while(total < power_2(undefined_rank,2*(curr_iteration+1)) * (ACC_RATE)){
      total += convoluted_distribution[arr_idx];
      arr_idx--;
      if(arr_idx > num_possible_values_max-1) break;
    }
    return arr_idx + 1 + min_val;*/
  }

  constexpr void self_convolute(const uint8_t curr_iteration){
    //convolution 2 distributions
    std::array<uint32_t,num_possible_values_max> new_convoluted_distribution{};
    for(uint16_t i = 0; i < num_possible_values_single*curr_iteration; i++){
      for(uint8_t j = 0; j < num_possible_values_single; j++){
        new_convoluted_distribution[i+j] += convoluted_distribution[i] * original_distribution[j];
      }
    }
    convoluted_distribution = new_convoluted_distribution;

    //threshold
    threshold[curr_iteration] = calc_approx_threshold(curr_iteration);
  }

  public:
  constexpr Distribution(){
    constexpr const auto score_matrix = sc.score_matrix;
    constexpr const auto min_val = sc.smallest_score;

    for(uint8_t i = 0; i < undefined_rank; i++){
      for(uint8_t j = 0; j < undefined_rank; j++){
        original_distribution[score_matrix[i][j]-min_val]++;
        convoluted_distribution[score_matrix[i][j]-min_val]++;
      }
    }

    for(uint8_t iteration = 1; iteration < max_num_of_convolution; iteration++){
      self_convolute(iteration);
    }
  };

  constexpr int8_t threshold_get(const uint8_t num_convolution) const {
    return threshold[num_convolution-1];
  }
};

template<class ScoreClass, const uint8_t max_num_convolution>
class Distribution_2{
  static constexpr const uint8_t num_distributions = max_num_convolution+1;
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const uint8_t num_vals_single = sc.highest_score - sc.smallest_score + 1;
  static constexpr const uint8_t num_vals_max = num_vals_single * num_distributions;

  std::array<uint32_t,num_vals_single> original_distribution{};
  std::array<uint32_t,num_vals_max> convoluted_distribution{};
  
  public:
  constexpr Distribution_2(){
    constexpr const auto score_matrix = sc.score_matrix;
    constexpr const auto min_val = sc.smallest_score;

    for(uint8_t i = 0; i < undefined_rank; i++){
      for(uint8_t j = 0; j < undefined_rank; j++){
        original_distribution[score_matrix[i][j]-min_val]++;
      }
    }

    uint8_t wheel[num_distributions]{};
    uint8_t loopid = max_num_convolution;

    convoluted_distribution[0] = power_2(original_distribution[0],num_distributions);

    uint8_t score = 0;
    uint32_t count = 1;
    
    for (;;){
      wheel[loopid]++;
      if(wheel[loopid] == num_vals_single){
        wheel[loopid] = 0;
        if(loopid == 0) break;
        loopid--;
      }
      else {
        score = 0;
        count = 1;
        for(uint8_t i = 0; i < num_distributions; i++){
          score += wheel[i];
          count *= original_distribution[wheel[i]];
        }
        convoluted_distribution[score] += count;

        if(loopid < max_num_convolution) loopid = max_num_convolution;
      }
    }
  };

  int8_t threshold() const {
    constexpr const int8_t min_val = sc.smallest_score * num_distributions;
    size_t rel_score = num_vals_max - 1;
    size_t total = 0;

    for(uint8_t i = num_vals_max - 1; i < num_vals_max; i--){
      std::cout << (int)convoluted_distribution[i] << std::endl;
    }

    const size_t cutoff = static_cast<size_t>(constexpr_pow(undefined_rank,2*num_distributions) * THRESHOLD_FACTOR / 1000000);

    for(; rel_score < num_vals_max; rel_score--){
      total += convoluted_distribution[rel_score];
      if(total > cutoff) break;
    }

    return rel_score + 1 + min_val;
  }
};
#endif