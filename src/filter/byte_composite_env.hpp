#ifndef BYTE_COMPOSITE_ENV_HPP
#define BYTE_COMPOSITE_ENV_HPP

#ifndef MAX_SUBQGRAM_LENGTH
#define MAX_SUBQGRAM_LENGTH 3
#endif

#ifdef __SSSE3__
#include <xmmintrin.h>
#include <x86intrin.h>
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

#include <thread>
#include <iostream>
#include <vector>

/*
Difference: non predetermined threshold
*/
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

  int8_t threshold;
  std::array<size_t,num_of_primary_env> env_size;
  
  //main_env
  //std::array<std::vector<ScoreQgramcodePair2>,num_of_primary_env> qgram_env_group;

  //local score correction
  std::array<double,undefined_rank> background_correction;
  //double local_score_correction;

  //state of qgram and permutation
  static constexpr const size_t simd_vec_len = 16; // based on __mm_shuffle_epi8 instruction

  //mutex
  std::mutex mut;

  //stats
  //static constexpr const size_t stats_size = constexpr_pow(undefined_rank,2);

  void invert_permutation(std::array<uint8_t,simd_vec_len>& perm) const {
    uint8_t temp_p[weight];
    for(size_t i = 0; i < weight; i++){
      temp_p[perm[i]] = i;
    }
    for(size_t i = 0; i < weight; i++){
      perm[i] = temp_p[i];
    }
  }

  void reshuffle_no_simd(std::array<uint8_t,simd_vec_len>& qgram, std::array<uint8_t,simd_vec_len>& perm) const {
    std::array<uint8_t,simd_vec_len> transformed_qgram{};
    for(uint8_t idx = 0; idx < weight; idx++)
    {
      transformed_qgram[idx] = qgram[perm[idx]];
    }
    qgram = transformed_qgram;
  }

  void reshuffle_with_simd(std::array<uint8_t,simd_vec_len>& qgram, std::array<uint8_t,simd_vec_len>& perm){
#if defined __SSSE3__
    __m128i loaded_qgram = _mm_loadu_si128((__m128i*)(qgram.data()));
    __m128i loaded_perm = _mm_loadu_si128((__m128i*)(perm.data()));
    _mm_storeu_si128((__m128i*)(qgram.data()),_mm_shuffle_epi8(loaded_qgram,loaded_perm));
#endif
  }

  template<const uint8_t qgram_length>
  const ScoreQgramcodePair2* const env_get_param(const uint16_t qgram_idx) const {
    static const EnvMatrix2<ScoreClass,qgram_length> env{};
    return env.sorted_env_get(qgram_idx);
  }
  
  template<const uint8_t env_idx>
  const uint8_t* const unsorted_subqgram_get(uint32_t qgram_code) const {
    static constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length_arr[env_idx]> unsorted_q{};
    return unsorted_q.qgram_get(qgram_code);
  }

  void reconstruct_qgram(std::array<uint16_t,num_of_primary_env> qgram_codes,uint8_t* qgram) const {
    uint8_t i = 0;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      constexpr const uint8_t qgram_length = qgram_length_arr[env_idx];
      const uint8_t* const subqgram = unsorted_subqgram_get<env_idx>(qgram_codes[env_idx]);
      constexpr_for<0,qgram_length,1>([&] (auto qgram_idx)
      {
        assert(subqgram[qgram_idx] < undefined_rank);
        qgram[i] = subqgram[qgram_idx];
        i++;
      });
    });
  }

  uint64_t qgram2code(const uint8_t* _qgram) const {
    uint64_t code = 0;
    for(uint8_t idx = 0; idx < weight; idx++){
      code *= undefined_rank;
      code += _qgram[idx];
    }
    return code;
  }
  
  template<const uint8_t sizeof_query_unit>
  void direct_add(std::array<uint8_t,simd_vec_len>& qgram, std::array<uint8_t,simd_vec_len>& perm,
                  const std::array<size_t,num_of_primary_env>& qgram_codes,
                  const size_t seqnum,const size_t position,const bool sorted, 
                  const bool with_simd,const double corrected_threshold,
                  std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                  const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    //std::cout << "Init env: " << seqnum << '\t' << position << std::endl;
    //Init all primary environments
    std::array<const ScoreQgramcodePair2*,num_of_primary_env> qgram_env_group;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<qgram_length_arr[env_idx]>(qgram_codes[env_idx]);
    });
    constexpr const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};
    
    //std::cout << "Iterating: " << seqnum << '\t' << position << std::endl;
    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }
    if(score < corrected_threshold) return;
    //std::cout << "Start" << std::endl;

    invert_permutation(perm);
    reconstruct_qgram(unsorted_qgram_codes,qgram.data());
    if(!sorted and !with_simd) reshuffle_no_simd(qgram,perm);
    else if(with_simd) reshuffle_with_simd(qgram,perm);
    uint64_t code = qgram2code(qgram.data());
    //std::cout << code << '\t' << (int) score << std::endl;
    //mut.lock();
    query_vec.emplace_back(query_packer,std::array<uint64_t,3>{
                            static_cast<uint64_t>(code),
                            static_cast<uint64_t>(seqnum),
                            static_cast<uint64_t>(position)});
    //mut.unlock();
    
    for(;;){
      score -= qgram_env_group[loopid][wheel[loopid]].score;
      wheel[loopid]++;
      if(wheel[loopid] >= env_size[loopid]){
        wheel[loopid] = 0;
        score += qgram_env_group[loopid][0].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][0].code;

        if(loopid == 0) break;
        loopid--;
      }else{
        score += qgram_env_group[loopid][wheel[loopid]].score;
        if(score >= corrected_threshold){
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          //std::cout << wheel[0] << std::endl;
          reconstruct_qgram(unsorted_qgram_codes,qgram.data());
          if(!sorted and !with_simd) reshuffle_no_simd(qgram,perm);
          else if(with_simd) reshuffle_with_simd(qgram,perm);
          code = qgram2code(qgram.data());
          //std::cout << code << '\t' << (int) score << std::endl;
          //mut.lock();
          query_vec.emplace_back(query_packer,std::array<uint64_t,3>{
                                  static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)});
          //mut.unlock();
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=env_size[last_loop]-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }
        if(loopid < last_loop) loopid = last_loop;
      }
    }
  }

  template<const uint8_t qgram_length>
  const ScoreQgramcodePair2* const env_get_mmseqs(const uint16_t qgram_idx) const {
    static const FullMatrix<ScoreClass,qgram_length> env{};
    return env.sorted_env_get(qgram_idx);
  }

  template<const uint8_t sizeof_query_unit>
  void direct_add_mmseqs(std::array<uint8_t,simd_vec_len>& qgram,
                        const std::array<size_t,num_of_primary_env>& qgram_codes,
                        const size_t seqnum,const size_t position,const double corrected_threshold,
                        std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                        const GttlBitPacker<sizeof_query_unit,3>& query_packer){
                        //,std::array<FilterStats<size_t>,stats_size>& stats){
    //std::cout << "Init env: " << seqnum << '\t' << position << std::endl;
    //Init all primary environmentss
    std::array<const ScoreQgramcodePair2*,num_of_primary_env> qgram_env_group;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_mmseqs<qgram_length_arr[env_idx]>(qgram_codes[env_idx]);
    });
    constexpr const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};
    
    //std::cout << "Iterating: " << seqnum << '\t' << position << std::endl;
    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }

    //stats[unsorted_qgram_codes[0]].appears += 1;
    //stats[unsorted_qgram_codes[1]].appears += 1;
    if(score < corrected_threshold) return;

    //size_t length = 1,count=1;

    reconstruct_qgram(unsorted_qgram_codes,qgram.data());
    uint64_t code = qgram2code(qgram.data());
    //std::cout << code << '\t' << (int) score << std::endl;

    //mut.lock();
    query_vec.emplace_back(query_packer,std::array<uint64_t,3>{
                            static_cast<uint64_t>(code),
                            static_cast<uint64_t>(seqnum),
                            static_cast<uint64_t>(position)});
    //mut.unlock();
    
    for(;;){
      score -= qgram_env_group[loopid][wheel[loopid]].score;
      wheel[loopid]++;
      if(wheel[loopid] >= env_size[loopid]){
        wheel[loopid] = 0;
        score += qgram_env_group[loopid][0].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][0].code;
        if(loopid == 0) break;
        loopid--;
      }else{
        score += qgram_env_group[loopid][wheel[loopid]].score;
        if(score >= corrected_threshold){
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          reconstruct_qgram(unsorted_qgram_codes,qgram.data());
          code = qgram2code(qgram.data());
          //std::cout << code << '\t' << (int) score << std::endl;
          //mut.lock();
          query_vec.emplace_back(query_packer,std::array<uint64_t,3>{
                                  static_cast<uint64_t>(code),
                                  static_cast<uint64_t>(seqnum),
                                  static_cast<uint64_t>(position)});
          //mut.unlock();
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=env_size[last_loop]-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }
        if(loopid < last_loop) loopid = last_loop;
      }
    }
  }

  double score_correction_set(const char* seq, const size_t seq_len, size_t position){
    double local_score_correction = 0;
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
    return local_score_correction;
  }

  public:
  BytesCompositeEnvironment2(){
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      env_size[i] = constexpr_pow(undefined_rank,qgram_length_arr[i]);
    }
  };

  void set_background_data(const std::array<uint64_t,alpha_size+1>& target_distribution,const double sensitivity){
    background_correction_set(target_distribution);
    const BGDistribution<ScoreClass,weight> distribution{};
    threshold = distribution.custom_threshold_get(sensitivity);
  }

  void set_background_data(const std::array<uint64_t,alpha_size+1>& target_distribution,const int64_t _threshold){
    background_correction_set(target_distribution);
    threshold = _threshold;
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
  void process_seed(std::array<uint8_t,simd_vec_len>& qgram, std::array<uint8_t,simd_vec_len>& perm,
                    const char* seq,const size_t seqnum, const size_t seq_len,const size_t position, 
                    GTTL_UNUSED const bool with_simd,const bool correct,const double correct_ratio,
                    std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                    const GttlBitPacker<sizeof_query_unit,3>& query_packer) {
    bool sorted;
    const auto qgram_codes = spaced_seed_encoder.encode(seq,qgram.data(),perm.data(),sorted);
    const double local_score_correction = score_correction_set(seq,seq_len,position);
    const double corrected_threshold = threshold - local_score_correction * correct * correct_ratio;
#if defined __SSSE3__
    direct_add<sizeof_query_unit>(qgram,perm,qgram_codes,seqnum,position,sorted,with_simd,corrected_threshold,query_vec,query_packer);
#else
    direct_add<sizeof_query_unit>(qgram,perm,qgram_codes,seqnum,position,sorted,false,corrected_threshold,query_vec,query_packer);
#endif
  }

  template<const uint8_t sizeof_query_unit>
  void process_seed_mmseqs(std::array<uint8_t,simd_vec_len>& qgram,
                          const char* seq,const size_t seqnum, const size_t seq_len,const size_t position, 
                          const bool correct,const double correct_ratio, 
                          std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
                          const GttlBitPacker<sizeof_query_unit,3>& query_packer){
                          //,std::array<FilterStats<size_t>,stats_size>& stats) {
    const auto qgram_codes = spaced_seed_encoder.encode_unsorted(seq);
    const double local_score_correction = score_correction_set(seq,seq_len,position);
    const double corrected_threshold = threshold - local_score_correction * correct * correct_ratio;
    direct_add_mmseqs<sizeof_query_unit>(qgram,qgram_codes,seqnum,position,corrected_threshold,query_vec,query_packer);
  }

  template<const uint8_t sizeof_query_unit>
  void process_multiseq(std::array<uint8_t,simd_vec_len>& qgram, std::array<uint8_t,simd_vec_len>& perm, 
    const GttlMultiseq* query,std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
    const GttlBitPacker<sizeof_query_unit,3>& query_packer, const bool mmseqs, const bool with_simd, 
    const bool correct,const double correct_ratio,const size_t start, const size_t end) {
    
    //std::array<FilterStats<size_t>,stats_size> stats{};
    //size_t old_size = 0;
    for(size_t seqnum = start; seqnum < end; seqnum++){
      const size_t seq_len = query->sequence_length_get(seqnum);
      if(seq_len >= span){
        const char* curr_seq = query->sequence_ptr_get(seqnum);
        if(!mmseqs){
          for(size_t i = 0; i < seq_len - span + 1; i++){
            process_seed<sizeof_query_unit>(qgram,perm,curr_seq+i,seqnum,seq_len,i,with_simd,correct,correct_ratio,query_vec,query_packer);
            /*std::cout << i << '\t' << query_vec.size() - old_size << std::endl;
            old_size = query_vec.size();*/
          }
        } else {
          for(size_t i = 0; i < seq_len - span + 1; i++){
            process_seed_mmseqs<sizeof_query_unit>(qgram,curr_seq+i,seqnum,seq_len,i,correct,correct_ratio,query_vec,query_packer);
            /*std::cout << i << '\t' << query_vec.size() - old_size << std::endl;
            old_size = query_vec.size();*/
          }
        }
      }
    }
    /*double sum = 0;
    size_t total_count = 0;
    for(size_t i = 0; i < stats.size(); i++){
      total_count += stats[i].appears;
      std::cout << (int) i << '\t' << stats[i].appears << '\t' << (int) stats[i].length << std::endl;
    }
    for(size_t i = 0; i < stats.size(); i++){
      std::cout << (int) i << '\t' << static_cast<double>(stats[i].appears)/total_count << std::endl;
      sum += static_cast<double>(stats[i].appears)/total_count * stats[i].length;
    }
    std::cout << sum << std::endl;*/
  }

  template<const uint8_t sizeof_query_unit>
  void process_pthread(const GttlMultiseq* query,std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
    const GttlBitPacker<sizeof_query_unit,3>& query_packer, const bool mmseqs, const bool with_simd, 
    const bool correct,const double correct_ratio,const size_t num_threads){
    
    std::cout << "Curr threshold: " << (int)threshold << std::endl;
    const size_t total_seq_num = query->sequences_number_get();
    
    
    const size_t total_thread_num = (total_seq_num < num_threads) ? total_seq_num : num_threads;

    const std::array<uint8_t,simd_vec_len> ref_qgram{};
    std::array<uint8_t,simd_vec_len> ref_perm{};
    for(size_t i = weight; i < simd_vec_len; i++){
      ref_perm[i] = 0x80;
    }

    std::vector<std::vector<BytesUnit<sizeof_query_unit,3>>> threads_vec{total_thread_num,query_vec};

    std::vector<std::array<uint8_t,simd_vec_len>> threads_qgram{total_thread_num,ref_qgram};
    std::vector<std::array<uint8_t,simd_vec_len>> threads_perm{total_thread_num,ref_perm};

    std::vector<std::thread> threads;
    std::vector<size_t> thread_work_idx{total_thread_num+1,0};
      
    thread_work_idx[0] = 0;
    const size_t max_thread_load = total_seq_num/total_thread_num + bool(total_seq_num % total_thread_num);
    size_t overcount = (total_thread_num-total_seq_num%total_thread_num)*(total_seq_num>total_thread_num)%total_thread_num;

    size_t curr_load = max_thread_load;
    for(size_t i = 1; i < total_thread_num; i++){
      if(overcount){
        overcount--;
        curr_load--;
      }
      assert(curr_load < total_seq_num); 
      thread_work_idx[i] = curr_load;
      curr_load += max_thread_load;
    }

    thread_work_idx[total_thread_num] = total_seq_num;

    for(uint8_t thread_idx = 0; thread_idx < total_thread_num; thread_idx++){
      const size_t start_idx = thread_work_idx[thread_idx];
      const size_t end_idx = thread_work_idx[thread_idx+1];

      threads.push_back(std::thread(&BytesCompositeEnvironment2::process_multiseq<sizeof_query_unit>,
                                    this,std::ref(threads_qgram[thread_idx]),std::ref(threads_perm[thread_idx]),
                                    query,std::ref(threads_vec[thread_idx]),std::ref(query_packer),mmseqs,with_simd,correct,
                                    correct_ratio,start_idx,end_idx));
    }
    for(uint8_t thread_idx = 0; thread_idx < total_thread_num; thread_idx++){
      threads[thread_idx].join();
      query_vec.insert(query_vec.end(),
                      std::make_move_iterator(threads_vec[thread_idx].begin()),
                      std::make_move_iterator(threads_vec[thread_idx].end()));
    }
  }

  template<const uint8_t sizeof_query_unit>
  void sort(std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,const uint64_t hash_bits){
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