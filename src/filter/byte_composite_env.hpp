#ifndef BYTE_COMPOSITE_ENV_HPP
#define BYTE_COMPOSITE_ENV_HPP

#ifndef MAX_SUBQGRAM_LENGTH
#define MAX_SUBQGRAM_LENGTH 3
#endif

#include <list>
#include "filter/qgram_environment.hpp"
#include "utilities/constexpr_for.hpp"
#include "utilities/unused.hpp"
#include "filter/spaced_seeds.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "filter/env_matrix.hpp"

#ifndef CORRECTION_REGION_SIZE
#define CORRECTION_REGION_SIZE 20
#endif

#ifdef __SSSE3__
#include <xmmintrin.h>
#include <x86intrin.h>
#endif

#include <iostream>
#include <vector>
template<class ScoreClass, const size_t seed>
class BytesCompositeEnvironment {
  //ScoreClass
  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const auto alpha_size = alpha.size();
  
  //SpacedSeedEncoder
  static constexpr const SpacedSeedEncoder<ScoreClass,seed,MAX_SUBQGRAM_LENGTH> spaced_seed_encoder{};
  static constexpr const auto seed_bitset = spaced_seed_encoder.seed_bitset_get();
  static constexpr const size_t span = spaced_seed_encoder.span_get();
  static constexpr const auto weight = spaced_seed_encoder.weight_get();
  static constexpr const auto num_of_primary_env = spaced_seed_encoder.num_of_primary_env_get();
  static constexpr const auto qgram_length_arr = spaced_seed_encoder.subqgram_length_arr_get();
  static constexpr const auto env_threshold_arr = spaced_seed_encoder.env_threshold_arr_get();
  static constexpr const int8_t threshold = spaced_seed_encoder.threshold_get();
  static constexpr const int8_t max_possible_score = weight*sc.highest_score;
  
  //main_env
  std::array<std::vector<ScoreQgramcodePair>,num_of_primary_env> qgram_env_group{};

  //local score correction
  std::array<float,undefined_rank> background_correction;
  float local_score_correction;

  std::array<uint8_t,weight>
  reshuffle_no_simd(const std::array<uint8_t,weight>& reconstructed_qgram,
                    const std::array<uint8_t,weight>& permutation) const {
    std::array<uint8_t,weight> transformed_qgram{};
    for(uint8_t idx = 0; idx < weight; idx++)
    {
      transformed_qgram[permutation[idx]] = reconstructed_qgram[idx];
    }
    return transformed_qgram;
  }

  template<const uint8_t env_idx>
  std::vector<ScoreQgramcodePair> env_get_param(const uint16_t qgram_idx) const {
    static const QgramEnvironment<ScoreClass,qgram_length_arr[env_idx],env_threshold_arr[env_idx]> env{};
    return env.qgram_env(qgram_idx);
  }

  template<const uint8_t qgram_length,const int8_t threshold>
  std::vector<ScoreQgramcodePair> env_get(const uint16_t qgram_idx) const {
    static const QgramEnvironment<ScoreClass,qgram_length,threshold> env{};
    return env.qgram_env(qgram_idx);
  }
  
  template<const uint8_t env_idx>
  std::array<uint8_t,qgram_length_arr[env_idx]> unsorted_subqgram_get(uint16_t qgram_code) const {
    static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length_arr[env_idx]> unsorted_q{};
    return unsorted_q.qgram_get(qgram_code);
  }

  std::array<uint8_t,weight> reconstruct_qgram(std::array<uint16_t,num_of_primary_env> qgram_codes){
    std::array<uint8_t,weight> reconstructed_qgram{};
    uint8_t i = 0;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      constexpr const uint8_t qgram_length = qgram_length_arr[env_idx];
      const auto subqgram = unsorted_subqgram_get<env_idx>(qgram_codes[env_idx]);
      constexpr_for<0,qgram_length,1>([&] (auto qgram_idx)
      {
        assert(subqgram[qgram_idx] < undefined_rank);
        reconstructed_qgram[i] = subqgram[qgram_idx];
        i++;
      });
    });

    return reconstructed_qgram;
  }

  uint64_t qgram2code(std::array<uint8_t,weight> qgram) const {
    uint64_t code = 0;
    for(uint8_t idx = 0; idx < weight; idx++){
      code *= undefined_rank;
      code += qgram[idx];
    }
    return code;
  }
  
  template<const uint8_t sizeof_query_unit>
  void direct_add(const std::array<uint16_t,num_of_primary_env>& qgram_codes,
                  const std::array<uint8_t,weight>& permutation,const bool sorted,
                  const size_t seqnum,const size_t position,
                  std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                  const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<env_idx>(qgram_codes[env_idx]);
    });

    for(size_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      if(qgram_env_group[env_idx].empty()) return;
    }

    const int8_t corrected_threshold = threshold - static_cast<int8_t>(local_score_correction);

    const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};
    uint64_t code;

    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }
    if(score < corrected_threshold) return;
    //std::cout << "Start" << std::endl;
    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
    if(!sorted) reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    code = qgram2code(reconstructed_qgram);
    /*
    query_vec.emplace_back({query_packer,{
                            static_cast<uint64_t>(code),
                            static_cast<uint64_t>(seqnum),
                            static_cast<uint64_t>(position)}});*/
    const BytesUnit<sizeof_query_unit,3> bu{
              query_packer,
                                    {static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}
            };
    query_vec.push_back(bu);
    for(;;){
      score -= qgram_env_group[loopid][wheel[loopid]].score;
      wheel[loopid]++;
      if(wheel[loopid] >= qgram_env_group[loopid].size()){
        wheel[loopid] = 0;
        score += qgram_env_group[loopid][0].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][0].code;

        if(loopid == 0) break;
        loopid--;
      }else{
        score += qgram_env_group[loopid][wheel[loopid]].score;
        if(score >= corrected_threshold){
          //std::cout << (int) score << '\t' << (int) corrected_threshold << std::endl;
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          
          auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
          if(!sorted) reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          code = qgram2code(reconstructed_qgram);
          //std::cout << (int) code << std::endl;
          const BytesUnit<sizeof_query_unit,3> bu{
              query_packer,
                                    {static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}
            };
          query_vec.push_back(bu);
          /*query_vec.emplace_back({query_packer,{
                                  static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}});*/
          //std::cout << (int) constructed_env.size() << std::endl;
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=qgram_env_group[last_loop].size()-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }
        if(loopid < last_loop) loopid = last_loop;
      }
    }
  }

  void score_correction_set(const char* seq, const size_t seq_len, size_t position){
    local_score_correction = 0;
    for(uint8_t i = span-1; i < span; i-- && position++){
      if(!seed_bitset[i]) continue;

      size_t start,end;
      start = (position < CORRECTION_REGION_SIZE) ? 0 : (position - CORRECTION_REGION_SIZE);
      end = (position >= seq_len - CORRECTION_REGION_SIZE - 1) ? (seq_len - 1) : (position + CORRECTION_REGION_SIZE);
      const size_t length = end - start; 
      //std::cout << (int) start << '\t' << (int) position << '\t' << (int) *(seq+position) << '\t' << (int) end << std::endl;
      float sum = 0;
      for(size_t seq_idx = start; seq_idx <= end; seq_idx++){
        if(seq_idx == position) continue;

        sum += sc.score_matrix[seq[position]][seq[seq_idx]];
      }

      local_score_correction -= sum/length;
      local_score_correction += background_correction[seq[position]];
    }
  }

  public:
  BytesCompositeEnvironment(){};

  void set_background_data(const std::array<uint64_t,alpha_size+1>& target_distribution){
    background_correction_set(target_distribution);
  }

  void background_correction_set(const std::array<uint64_t,alpha_size+1>& target_distribution) {
    //calculate frequency of background data
    std::array<double,undefined_rank> frequency{};
    size_t sum_count = 0;
    for(uint8_t i = 0; i < alpha_size; i++){
      sum_count += target_distribution[i];
    }
    //std::cout << (int) sum_count << std::endl;
    for(uint8_t i = 0; i < undefined_rank; i++){
      frequency[i] = target_distribution[i]/sum_count;
    }

    //calculate background correction 
    for(uint8_t i = 0; i < undefined_rank; i++){
      double sum_char = 0;
      for(uint8_t j = 0; j < undefined_rank; j++){
        sum_char += frequency[j] * sc.score_matrix[i][j];
      }
      background_correction[i] = sum_char/undefined_rank;
    }
  }

  template<const uint8_t sizeof_query_unit>
  void process_seed(const char* seq,const size_t seqnum, const size_t seq_len,const size_t position,
                    std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                    const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    const auto encode_info = spaced_seed_encoder.encode(seq);
    score_correction_set(seq,seq_len,position);
    direct_add<sizeof_query_unit>(encode_info.codes,encode_info.permutation,
      encode_info.sorted,seqnum,position,query_vec,query_packer);
/*#ifdef __SSSE3__
    direct_add(encode_info.codes,encode_info.permutation,encode_info.sorted,seqnum,
              position,with_simd);
#else
    direct_add(encode_info.codes,encode_info.permutation,encode_info.sorted,seqnum,
              position,false);
#endif*/
  }

  template<const uint8_t sizeof_query_unit>
  void process_multiseq(const GttlMultiseq* query,std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
  const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    
    const size_t total_seq_num = query->sequences_number_get();
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++){
      const size_t seq_len = query->sequence_length_get(seqnum);
      if(seq_len >= span){
        const char* curr_seq = query->sequence_ptr_get(seqnum);
        for(size_t i = 0; i < seq_len - span + 1; i++){
          process_seed<sizeof_query_unit>(curr_seq+i,seqnum,seq_len,i,query_vec,query_packer);
        }
      }
    }
  }

  template<const uint8_t sizeof_query_unit>
  void dbsort(std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,const uint64_t hash_bits){
    if(sizeof_query_unit == 8){
      ska_lsb_radix_sort<size_t>(hash_bits,
                              reinterpret_cast<uint64_t *>
                              (query_vec.data()),
                              query_vec.size());
    } else {
      ska_large_lsb_small_radix_sort(sizeof_query_unit,hash_bits,
                                    reinterpret_cast<uint8_t *>(query_vec.data()),
                                    query_vec.size(),false);
    }
  }

  constexpr uint8_t span_get() const {
    return span;
  }

  constexpr uint8_t weight_get() const {
    return weight;
  }
};

template<class ScoreClass, const size_t seed>
class BytesCompositeEnvironment2 {
  //ScoreClass
  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const auto alpha_size = alpha.size();
  
  //SpacedSeedEncoder
  static constexpr const SpacedSeedEncoder<ScoreClass,seed,MAX_SUBQGRAM_LENGTH> spaced_seed_encoder{};
  static constexpr const auto seed_bitset = spaced_seed_encoder.seed_bitset_get();
  static constexpr const size_t span = spaced_seed_encoder.span_get();
  static constexpr const auto weight = spaced_seed_encoder.weight_get();
  static constexpr const auto num_of_primary_env = spaced_seed_encoder.num_of_primary_env_get();
  static constexpr const auto qgram_length_arr = spaced_seed_encoder.subqgram_length_arr_get();
  static constexpr const int8_t max_possible_score = weight*sc.highest_score;

  static constexpr const uint64_t unsorted_size = constexpr_pow(undefined_rank,weight);

  int8_t threshold;
  
  //main_env
  std::array<std::vector<ScoreQgramcodePair2>,num_of_primary_env> qgram_env_group;

  //local score correction
  std::array<float,undefined_rank> background_correction;
  float local_score_correction;

  std::array<uint8_t,weight>
  reshuffle_no_simd(const std::array<uint8_t,weight>& reconstructed_qgram,
                    const std::array<uint8_t,weight>& permutation) const {
    std::array<uint8_t,weight> transformed_qgram{};
    for(uint8_t idx = 0; idx < weight; idx++)
    {
      transformed_qgram[permutation[idx]] = reconstructed_qgram[idx];
    }
    return transformed_qgram;
  }

  template<const uint8_t qgram_length>
  std::vector<ScoreQgramcodePair2> env_get_param(const uint16_t qgram_idx) const {
    static const EnvMatrix2<ScoreClass,qgram_length> env{};
    return env.sorted_env_get(qgram_idx);
  }
  
  template<const uint8_t env_idx>
  std::array<uint8_t,qgram_length_arr[env_idx]> unsorted_subqgram_get(uint16_t qgram_code) const {
    static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length_arr[env_idx]> unsorted_q{};
    return unsorted_q.qgram_get(qgram_code);
  }

  std::array<uint8_t,weight> reconstruct_qgram(std::array<uint16_t,num_of_primary_env> qgram_codes){
    std::array<uint8_t,weight> reconstructed_qgram{};
    uint8_t i = 0;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      constexpr const uint8_t qgram_length = qgram_length_arr[env_idx];
      const auto subqgram = unsorted_subqgram_get<env_idx>(qgram_codes[env_idx]);
      constexpr_for<0,qgram_length,1>([&] (auto qgram_idx)
      {
        assert(subqgram[qgram_idx] < undefined_rank);
        reconstructed_qgram[i] = subqgram[qgram_idx];
        i++;
      });
    });

    return reconstructed_qgram;
  }

  uint64_t qgram2code(std::array<uint8_t,weight> qgram) const {
    uint64_t code = 0;
    for(uint8_t idx = 0; idx < weight; idx++){
      code *= undefined_rank;
      code += qgram[idx];
    }
    return code;
  }
  
  template<const uint8_t sizeof_query_unit>
  void direct_add(const std::array<uint16_t,num_of_primary_env>& qgram_codes,
                  const std::array<uint8_t,weight>& permutation,const bool sorted,
                  const size_t seqnum,const size_t position,
                  std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                  const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    //std::cout << "Init env: " << seqnum << '\t' << position << std::endl;
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<qgram_length_arr[env_idx]>(qgram_codes[env_idx]);
    });

    const int8_t corrected_threshold = threshold - static_cast<int8_t>(local_score_correction);
    const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};
    uint64_t code;
    //std::cout << "Iterating: " << seqnum << '\t' << position << std::endl;
    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }
    if(score < corrected_threshold) return;
    //std::cout << "Start" << std::endl;
    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
    if(!sorted) reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    code = qgram2code(reconstructed_qgram);
    /*
    query_vec.emplace_back({query_packer,{
                            static_cast<uint64_t>(code),
                            static_cast<uint64_t>(seqnum),
                            static_cast<uint64_t>(position)}});*/
    const BytesUnit<sizeof_query_unit,3> bu{
              query_packer,
                                    {static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}
            };
    query_vec.push_back(bu);
    
    for(;;){
      score -= qgram_env_group[loopid][wheel[loopid]].score;
      wheel[loopid]++;
      if(wheel[loopid] >= qgram_env_group[loopid].size()){
        wheel[loopid] = 0;
        score += qgram_env_group[loopid][0].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][0].code;

        if(loopid == 0) break;
        loopid--;
      }else{
        score += qgram_env_group[loopid][wheel[loopid]].score;
        if(score >= corrected_threshold){
          //std::cout << (int) score << '\t' << (int) corrected_threshold << std::endl;
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          
          auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
          if(!sorted) reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          code = qgram2code(reconstructed_qgram);
          //std::cout << (int) code << std::endl;
          const BytesUnit<sizeof_query_unit,3> bu{
              query_packer,
                                    {static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}
            };
          query_vec.push_back(bu);
          /*query_vec.emplace_back({query_packer,{
                                  static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)}});*/
          //std::cout << (int) constructed_env.size() << std::endl;
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=qgram_env_group[last_loop].size()-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }
        if(loopid < last_loop) loopid = last_loop;
      }
    }
  }

  void score_correction_set(const char* seq, const size_t seq_len, size_t position){
    local_score_correction = 0;
    for(uint8_t i = span-1; i < span; i-- && position++){
      if(!seed_bitset[i]) continue;

      size_t start,end;
      start = (position < CORRECTION_REGION_SIZE) ? 0 : (position - CORRECTION_REGION_SIZE);
      end = (position >= seq_len - CORRECTION_REGION_SIZE - 1) ? (seq_len - 1) : (position + CORRECTION_REGION_SIZE);
      const size_t length = end - start; 
      //std::cout << (int) start << '\t' << (int) position << '\t' << (int) *(seq+position) << '\t' << (int) end << std::endl;
      float sum = 0;
      for(size_t seq_idx = start; seq_idx <= end; seq_idx++){
        if(seq_idx == position) continue;

        sum += sc.score_matrix[seq[position]][seq[seq_idx]];
      }

      local_score_correction -= sum/length;
      local_score_correction += background_correction[seq[position]];
    }
  }

  public:
  BytesCompositeEnvironment2(){};

  void set_background_data(const std::array<uint64_t,alpha_size+1>& target_distribution,const double sensitivity){
    background_correction_set(target_distribution);
    const BGDistribution<ScoreClass,weight> distribution{};
    threshold = distribution.custom_threshold_get(sensitivity);
  }

  void background_correction_set(const std::array<uint64_t,alpha_size+1>& target_distribution) {
    //calculate frequency of background data
    std::array<double,undefined_rank> frequency{};
    size_t sum_count = 0;
    for(uint8_t i = 0; i < alpha_size; i++){
      sum_count += target_distribution[i];
    }
    //std::cout << (int) sum_count << std::endl;
    for(uint8_t i = 0; i < undefined_rank; i++){
      frequency[i] = target_distribution[i]/sum_count;
    }

    //calculate background correction 
    for(uint8_t i = 0; i < undefined_rank; i++){
      double sum_char = 0;
      for(uint8_t j = 0; j < undefined_rank; j++){
        sum_char += frequency[j] * sc.score_matrix[i][j];
      }
      background_correction[i] = sum_char/undefined_rank;
    }
  }

  template<const uint8_t sizeof_query_unit>
  void process_seed(const char* seq,const size_t seqnum, const size_t seq_len,const size_t position,
                    std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                    const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    const auto encode_info = spaced_seed_encoder.encode(seq);
    score_correction_set(seq,seq_len,position);
    direct_add<sizeof_query_unit>(encode_info.codes,encode_info.permutation,
      encode_info.sorted,seqnum,position,query_vec,query_packer);
  }

  template<const uint8_t sizeof_query_unit>
  void process_multiseq(const GttlMultiseq* query,std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
  const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    std::cout << "Curr threshold: " << (int)threshold << std::endl;
    const size_t total_seq_num = query->sequences_number_get();
    
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++){
      const size_t seq_len = query->sequence_length_get(seqnum);
      if(seq_len >= span){
        const char* curr_seq = query->sequence_ptr_get(seqnum);
        for(size_t i = 0; i < seq_len - span + 1; i++){
          process_seed<sizeof_query_unit>(curr_seq+i,seqnum,seq_len,i,query_vec,query_packer);
        }
      }
    }

    std::cout << "Query elem count: " << query_vec.size() << std::endl;
  }

  template<const uint8_t sizeof_query_unit>
  void dbsort(std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,const uint64_t hash_bits){
    if(sizeof_query_unit == 8){
      ska_lsb_radix_sort<size_t>(hash_bits,
                              reinterpret_cast<uint64_t *>
                              (query_vec.data()),
                              query_vec.size());
    } else {
      ska_large_lsb_small_radix_sort(sizeof_query_unit,hash_bits,
                                    reinterpret_cast<uint8_t *>(query_vec.data()),
                                    query_vec.size(),false);
    }
  }

  constexpr uint8_t span_get() const {
    return span;
  }

  constexpr uint8_t weight_get() const {
    return weight;
  }
};

#endif