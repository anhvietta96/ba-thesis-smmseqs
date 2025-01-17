#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include<array>
#include<cassert>
#include "filter/utils.hpp"
#include <math.h>
#include "filter/env_matrix.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/runtime_class.hpp"
#include "filter/spaced_seeds.hpp"

#ifndef THRESHOLD_FACTOR
#define THRESHOLD_FACTOR 1
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

template<class ScoreClass,const uint8_t weight>
class BGDistribution {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const size_t num_possible_values_single = sc.highest_score - sc.smallest_score + 1;
  static constexpr const size_t num_possible_values_max = num_possible_values_single*weight;
  static constexpr const int64_t smallest_score = sc.smallest_score * weight;

  public:
  /*BGDistribution(const std::array<uint64_t,undefined_rank+1>& target_distribution,
    const uint64_t num_bins,const uint64_t sensitivity){
    size_t sum_count = 0;
    for(uint8_t i = 0; i < alpha_size; i++){
      sum_count += target_distribution[i];
    }

    double lowest = 0, highest = 0;
    double corrected_score_matrix[undefined_rank][undefined_rank]{};

    for(size_t i = 0; i < undefined_rank; i++){
      for(size_t j = 0; j < undefined_rank; j++){
        corrected_score_matrix[i][j] = sc.score_matrix * target_distribution[i] 
        * target_distribution[j] / (sum_count * sum_count);

        if(lowest > corrected_score_matrix[i][j]) lowest = corrected_score_matrix[i][j];
        if(highest < corrected_score_matrix[i][j]) highest = corrected_score_matrix[i][j];
      }
    }

    std::array<uint64_t,num_bins> base_hist{};
    std::array<double,num_bins> intervals{};

    //calc_intervals(intervals,lowest,highest,__UINT64_MAX__);
    intervals[num_bins-1] = highest;
    double diff = (highest-lowest)/num_bins;
    for(size_t i = num_bins-1; i != 0; i--){
      intervals[i-1] = intervals[i] - diff;
    }
    const double base_diff = diff;
    const std::array<double,num_bins> base_intervals = intervals;

    for(size_t i = 0; i < undefined_rank; i++){
      for(size_t j = 0; j < undefined_rank; j++){
        for(size_t k = 0; k < num_bins; k++){
          if(intervals[k] > corrected_score_matrix[i][j]){
            base_hist[k]++;
            break;
          }
        }
      }
    }
    
    std::array<uint64_t,num_bins> hist = base_hist;

    for(uint8_t convolution_num = 0; convolution_num < weight-1; convolution_num++){
      std::array<double,num_bins> averages{};
      for(size_t i = num_bins-1; i != 0; i--){
        averages[i] = intervals[i] - diff/2;
      }
      averages[0] = intervals[0] - diff/2;
      
      std::array<double,num_bins> old_intervals = intervals;
      intervals[num_bins-1] = highest*(convolution_num+2);
      const double new_diff = (highest-lowest)*(convolution_num+2)/num_bins;
      for(size_t i = num_bins-2; i < num_bins; i--){
        intervals[i] = intervals[i+1] - diff;
      }
      
      std::array<double,num_bins> temp_hist{};
      for(size_t i = 0; i < num_bins; i++){
        for(size_t j = 0; j < num_bins; j++){
          const curr_score = (old_intervals[i] - diff/2) + (base_intervals[j] - base_diff/2);
          size_t idx = 0;
          for(;idx < num_bins;idx++) if(intervals[idx] > curr_score) break;
          temp_hist[idx] += hist[i] * base_hist[j];
        }
      }
      hist = temp_hist;
      diff = new_diff;  
    }

    extract_threshold(hist,sensitity);
    size_t idx = num_bins-1;
    
    for(;idx < num_bins;idx++){

    }
  }*/
  //BGDistribution(){};
  BGDistribution(){};

  int64_t custom_threshold_get(const double sensitivity) const {
    std::array<uint64_t,num_possible_values_max> base_hist{};
    for(size_t i = 0; i < undefined_rank; i++){
      for(size_t j = 0; j < undefined_rank; j++){
        base_hist[sc.score_matrix[i][j]-sc.smallest_score]++; 
      }
    }
    std::array<uint64_t,num_possible_values_max> hist = base_hist;

    for(uint8_t convolution_num = 1; convolution_num < weight; convolution_num++){
      std::array<uint64_t,num_possible_values_max> next_hist{};
      for(size_t i = 0; i < num_possible_values_single * convolution_num; i++){
        for(size_t j = 0; j < num_possible_values_single; j++){
          next_hist[i+j] += hist[i] * base_hist[j];
        }
      }
      hist = next_hist;
    }

    /*std::array<double,num_possible_values_max> test_hist;
    for(size_t i = 0; i < num_possible_values_max; i++){
      test_hist[i] = static_cast<double>(hist[i])/constexpr_pow(undefined_rank,weight);
      test_hist[i] /= constexpr_pow(undefined_rank,weight);
      std::cout << test_hist[i] << std::endl;
    }*/

    const uint64_t cutoff = constexpr_pow(undefined_rank,weight)*sensitivity;
    //std::cout << "Cutoff " << cutoff << std::endl;
    uint64_t count = 0;
    size_t idx = num_possible_values_max-1;
    for(;idx < num_possible_values_max;idx--){
      count += hist[idx];
      //std::cout << count << std::endl;
      if(count > cutoff) break;
    }

    return idx + sc.smallest_score * weight + 1;
  }

  int64_t custom_threshold_get(const std::array<uint64_t,undefined_rank+1>& target_distribution,
    const double sensitivity) const {
    
    size_t sum_count = 0;
    for(uint8_t i = 0; i < undefined_rank; i++){
      sum_count += target_distribution[i];
    }
    std::array<double,undefined_rank> freq{};
    for(uint8_t i = 0; i < undefined_rank; i++){
      freq[i] = static_cast<double>(target_distribution[i]) / sum_count;
      //std::cout << freq[i] << std::endl;
    }

    std::array<double,num_possible_values_max> base_hist{};
    for(size_t i = 0; i < undefined_rank; i++){
      for(size_t j = 0; j < undefined_rank; j++){
        base_hist[sc.score_matrix[i][j]-sc.smallest_score] += freq[i] * freq[j]; 
      }
    }
    std::array<double,num_possible_values_max> hist = base_hist;

    for(uint8_t convolution_num = 1; convolution_num < weight; convolution_num++){
      std::array<double,num_possible_values_max> next_hist{};
      for(size_t i = 0; i < num_possible_values_single * convolution_num; i++){
        for(size_t j = 0; j < num_possible_values_single; j++){
          next_hist[i+j] += hist[i] * base_hist[j];
        }
      }
      hist = next_hist;
    }
    /*std::cout << "Distribution" << std::endl;
    for(size_t i = 0; i < hist.size(); i++){
      std::cout << hist[i] << std::endl;
    }*/

    const double cutoff = sensitivity / constexpr_pow(undefined_rank,weight); 
    //std::cout << "Cutoff " << cutoff << std::endl;
    double count = 0;
    size_t idx = num_possible_values_max-1;
    for(;idx < num_possible_values_max;idx--){
      count += hist[idx];
      //std::cout << count << std::endl;
      if(count > cutoff) break;
    }

    return sc.smallest_score * weight + idx + 1;
  }

  int64_t custom_threshold_get2(const std::array<uint64_t,undefined_rank+1>& target_distribution,
    const double sensitivity) const {
    
    size_t sum_count = 0;
    for(uint8_t i = 0; i < undefined_rank; i++){
      sum_count += target_distribution[i];
    }
    std::array<double,undefined_rank> freq{};
    for(uint8_t i = 0; i < undefined_rank; i++){
      freq[i] = static_cast<double>(target_distribution[i]) / sum_count;
    }/*
    std::array<double,undefined_rank> freq2{};
    for(uint8_t i = 0; i < undefined_rank; i++){
      freq2[i] = static_cast<double>(1) / undefined_rank;
    }*/

    std::array<double,num_possible_values_max> base_hist{};
    for(size_t i = 0; i < undefined_rank; i++){
      for(size_t j = 0; j < undefined_rank; j++){
        base_hist[sc.score_matrix[i][j]-sc.smallest_score] += freq[i]; 
      }
    }
    std::array<double,num_possible_values_max> hist = base_hist;
    /*std::cout << "Base Distribution" << std::endl;
    for(size_t i = 0; i < base_hist.size(); i++){
      std::cout << base_hist[i] << std::endl;
    }*/

    for(uint8_t convolution_num = 1; convolution_num < weight; convolution_num++){
      std::array<double,num_possible_values_max> next_hist{};
      for(size_t i = 0; i < num_possible_values_single * convolution_num; i++){
        for(size_t j = 0; j < num_possible_values_single; j++){
          next_hist[i+j] += hist[i] * base_hist[j];
        }
      }
      hist = next_hist;
    }

    /*std::cout << "Distribution" << std::endl;
    for(size_t i = 0; i < hist.size(); i++){
      std::cout << hist[i] << std::endl;
    }*/

    const double cutoff = sensitivity;// / constexpr_pow(undefined_rank,weight); 
    //std::cout << "Cutoff " << cutoff << std::endl;
    double count = 0;
    size_t idx = num_possible_values_max-1;
    for(;idx < num_possible_values_max;idx--){
      count += hist[idx];
      //std::cout << count << std::endl;
      if(count > cutoff) break;
    }

    return sc.smallest_score * weight + idx + 1;
  }

  void get_next_qgram(std::array<uint8_t,4>& qgram)const{
    if(qgram[3] == undefined_rank-1){
      for(uint8_t i = 3; i < 4; i--){
        if(qgram[i] == undefined_rank-1) qgram[i] = 0;
        else{
          qgram[i]++;
          break;
        }
      }
    } else {
      qgram[3]++;
    }
  }

  template<const uint8_t qgram_length>
  const ScoreQgramcodePair2* const env_get(const uint64_t code) const {
    static const FullMatrix<ScoreClass,qgram_length> fm{};
    return fm.sorted_env_get(code);
  }

  template<const size_t sizeof_unit>
  int64_t context_sensitive_threshold_get(const std::vector<BytesUnit<sizeof_unit,3>>& hash_data,
  const GttlBitPacker<sizeof_unit,3>& packer, const double sensitivity) const {
    static constexpr const uint64_t seed = (1 << weight) - 1;
    static constexpr const SpacedSeedEncoder<ScoreClass,seed,3> spaced_seed_encoder{};
    static constexpr const auto num_of_primary_env = spaced_seed_encoder.num_of_primary_env_get();
    static constexpr const auto& subqgram_length_arr = spaced_seed_encoder.subqgram_length_arr_get();

    std::array<uint64_t,num_of_primary_env> modFactors{};
    std::array<uint64_t,num_of_primary_env> divFactors{};
    std::array<int64_t,num_of_primary_env> minArr{};

    for(size_t i = 0; i < num_of_primary_env; i++){
      modFactors[i] = constexpr_pow(undefined_rank,subqgram_length_arr[i]);
      uint64_t divFactor = 1;
      for(size_t j = i+1; j < num_of_primary_env; j++){
        divFactor *= constexpr_pow(undefined_rank,subqgram_length_arr[j]);
      }
      divFactors[i] = divFactor;
      minArr[i] = sc.smallest_score * subqgram_length_arr[i];
    }

    std::array<std::array<uint32_t,constexpr_pow(undefined_rank,3)>,3> code_count{};

    RunTimeClass rt{};
    const size_t total = hash_data.size();
    for(size_t i = 0; i < total; i++){
      const size_t code = hash_data[i].template decode_at<0>(packer);
      for(size_t j = 0; j < num_of_primary_env; j++){
        code_count[j][(code/divFactors[j])%modFactors[j]]++;
      }
    }

    std::array<std::array<double,constexpr_pow(undefined_rank,3)>,3> freq{};
    for(size_t i = 0; i < num_of_primary_env; i++){
      for(size_t j = 0; j < modFactors[i]; j++){
        freq[i][j] = static_cast<double>(code_count[i][j])/total;
      }
    }
    rt.show("Extracted from hash data");

    std::array<std::array<double,num_possible_values_max>,num_of_primary_env> hist{};
    constexpr_for<0,num_of_primary_env,1>([&] (auto k){
      for(size_t i = 0; i < modFactors[k]; i++){
        const auto& env = env_get<subqgram_length_arr[k]>(i);
        for(size_t j = 0; j < modFactors[k]; j++){
          hist[k][env[j].score-minArr[k]] += freq[k][i];
        }
      }
    });

    std::array<double,num_possible_values_max> convoluted = hist[0];
    size_t num = subqgram_length_arr[0];
    for(size_t k = 1; k < num_of_primary_env; k++){
      std::array<double,num_possible_values_max> temp_conv{}; 
      for(size_t i = 0; i < num_possible_values_single * num; i++){
        for(size_t j = 0; j < num_possible_values_single * subqgram_length_arr[k]; j++){
          temp_conv[i+j] += convoluted[i] * hist[k][j];
        }
      }
      convoluted = temp_conv;
      num += subqgram_length_arr[k];
    }
    rt.show("Constructed histogram");
    
    const double cutoff = sensitivity;// / constexpr_pow(undefined_rank,weight); 
    //std::cout << "Cutoff " << cutoff << std::endl;
    double count = 0;
    size_t idx = num_possible_values_max-1;
    for(;idx < num_possible_values_max;idx--){
      count += convoluted[idx];
      //std::cout << count << std::endl;
      if(count > cutoff) break;
    }

    return smallest_score + idx + 1;
  }
};

template<typename T>
struct FilterStats {
  T appears;
  size_t length;
  constexpr FilterStats(): appears(0),length(0){};
};

template<class ScoreClass, const size_t weight>
class ThresholdEvaluate {
  static constexpr const ScoreClass sc{};
  static constexpr const size_t undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;
  
  static constexpr const UnsortedQmer<char_spec,undefined_rank,weight> index_table{};
  const FullMatrix<ScoreClass,weight> m{};
  static constexpr const  size_t size = constexpr_pow(undefined_rank,weight);
  std::array<FilterStats<double>,size> stats;
  
  public:
  void eval(const std::array<size_t,undefined_rank+1>& target_distribution,const int64_t threshold){
    std::array<FilterStats<double>,size> _stats{};
    
    size_t sum_count = 0;
    for(uint8_t i = 0; i < undefined_rank; i++){
      sum_count += target_distribution[i];
    }
    std::array<double,undefined_rank> freq{};
    for(uint8_t i = 0; i < undefined_rank; i++){
      freq[i] = static_cast<double>(target_distribution[i]) / sum_count;
    }

    for(size_t qgram_code = 0; qgram_code < size; qgram_code++){
      const auto qgram_ref = index_table.qgram_get(qgram_code);
      double p = static_cast<double>(1);
      for(size_t i = 0; i < weight; i++){
        p *= freq[qgram_ref[i]];
      }
      _stats[qgram_code].appears = p;
      
      const auto env = m.sorted_env_get(qgram_code);
      size_t count = 0;
      for(size_t i = 0; i < env.size(); i++){
        if(env[i].score < threshold) break;
        count++;
      }
      _stats[qgram_code].length = count;
    }

    stats = _stats;
  }

  const std::array<FilterStats<double>,size>& stats_get() const {
    return stats;
  }
};
#endif