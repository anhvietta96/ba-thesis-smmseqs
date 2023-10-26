#ifndef COMPOSITE_ENV_HPP
#define COMPOSITE_ENV_HPP

#ifndef MAX_SUBQGRAM_LENGTH
#define MAX_SUBQGRAM_LENGTH 3
#endif

#include <list>
#include "filter/qgram_environment.hpp"
#include "utilities/constexpr_for.hpp"
#include "utilities/unused.hpp"
#include "filter/spaced_seeds.hpp"
#include <thread>
#include <omp.h>

#ifndef NUM_THREADS
#define NUM_THREADS 8
#endif

#define OUTPUT_QGRAM
#define VECTOR

#ifndef CORRECTION_REGION_SIZE
#define CORRECTION_REGION_SIZE 20
#endif

#ifdef __SSSE3__
#include <xmmintrin.h>
#include <x86intrin.h>
#endif

#include "utilities/runtime_class.hpp"

#include <iostream>
#include <vector>

//code_check
//ifndef
//blosum header
//AVX2

struct NumEnvElem {
  size_t position;
  int8_t score;
  size_t code;
  
  public:
  NumEnvElem(const size_t _position,const int8_t _score, const size_t _code): 
    position(_position),score(_score),code(_code){};
};

template<const uint8_t qgram_length>
struct LocalEnvElem {
  size_t position;
  int8_t score;
  std::array<uint8_t,qgram_length> qgram;

  LocalEnvElem(const size_t _position, const int8_t _score,
    const std::array<uint8_t,qgram_length>& _qgram): position(_position),score(_score),qgram(_qgram){};

  void set_qgram(const uint8_t* ref_qgram)
  {
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = ref_qgram[i];
    }
  }

  NumEnvElem to_num_type(const size_t undefined_rank) const {
    size_t code = 0;
    for(uint8_t i = 0; i < qgram_length; i++){
      code *= undefined_rank;
      code += qgram[i];
    }
    return NumEnvElem(position,score,code);
  }
};

struct ScorePosition {
  uint16_t start = 0;
  uint16_t end = 0;

  ScorePosition(): start{0},end{0}{};

  ScorePosition(const uint16_t& _start,const uint16_t& _end): start{_start},end{_end}{};

  void operator= (const ScorePosition& other){
    start = other.start;
    end = other.end;
  }
};

template<const uint8_t num_of_primary_env>
struct TempEnvElem {
  int8_t score = 0;
  std::array<ScorePosition,num_of_primary_env> positions{};

  TempEnvElem(const int8_t& _score,const std::array<ScorePosition,num_of_primary_env>& _sp): score{_score},positions{_sp}{};
};

struct EnvInfo {
  int8_t max_score;
  std::vector<ScorePosition> threshold_arr;

  EnvInfo(): max_score{0}, threshold_arr{}{};
};

template<class ScoreClass, const size_t seed>
class CompositeEnvironment {
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
  std::array<uint64_t,undefined_rank> background_correction;
  float local_score_correction;

#ifdef OUTPUT_QGRAM
  using OutputType = LocalEnvElem<weight>;
#else
  using OutputType = NumEnvElem;
#endif
  
  std::vector<OutputType> constructed_env{};
#ifdef VECTOR
  std::vector<OutputType> temp_env{};
  std::array<size_t,max_possible_score-threshold+1> indexer{};
#else
  std::list<OutputType> temp_env{};
  typename std::list<OutputType>::iterator indexer[max_possible_score-threshold+1];
#endif
  
  std::mutex env_mutex;

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
#ifdef __SSSE3__
  void reshuffle_with_simd(const std::array<uint8_t,weight>& permutation,const size_t position){
    size_t start_idx = constructed_env.size()-1;
    
    while(start_idx < constructed_env.size())
    {
      if(constructed_env[start_idx].position != position) break;
      start_idx--;
    }

    const size_t num_of_qgrams = constructed_env.size()-1-start_idx;
    if(num_of_qgrams == 0) return;

    constexpr const uint8_t qgram_per_simd_vector = 16 / weight;

    const size_t num_of_simd_operations = (num_of_qgrams) ? ((num_of_qgrams-1) /
                                           qgram_per_simd_vector + 1) : 0;
    size_t qgram_pos = start_idx+1;
    size_t curr_qgram = 0;

    //std::cout <<  (int) num_of_qgrams << '\t' << (int) qgram_pos << std::endl;
    
    //invert permutation
    std::array<uint8_t,weight> inv_permutation{};
    for(uint8_t i = 0; i < weight; i++)
    {
      inv_permutation[permutation[i]] = i;
    }

    //create __m128i for transformation
    uint8_t transformation[16]{};
    for(uint8_t i = 0; i < qgram_per_simd_vector; i++){
      for(uint8_t j = 0; j < weight; j++){
        transformation[i * weight + j] = inv_permutation[j] + i * weight;
      }
    }
    for(uint8_t i = 15; transformation[i] == 0; i--){
      transformation[i] = 0x80;
    }

    __m128i *transform_m128i = (__m128i *) transformation;
    __m128i loaded_transformation = _mm_load_si128(transform_m128i);
    
    for (size_t simd_op_no = 0; simd_op_no < num_of_simd_operations; simd_op_no++)
    {
      uint8_t buf[16]{};
      uint8_t filled_qgram = 0;
      while(curr_qgram < num_of_qgrams and filled_qgram < qgram_per_simd_vector)
      {
        //fill buffer
        std::array<uint8_t,weight> qgram = constructed_env[qgram_pos].qgram;
        for(uint8_t char_idx = 0; char_idx < weight; char_idx++)
        {
          assert(qgram[char_idx] < undefined_rank);
          buf[filled_qgram*weight+char_idx] = qgram[char_idx];
        }
        filled_qgram++;
        curr_qgram++;
        qgram_pos++;
      }
      std::cout << "Buffer\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)buf[j] << '\t';
      }
      std::cout << std::endl;
      __m128i *vector = (__m128i *) buf;
      __m128i loaded_vector = _mm_load_si128(vector);

      //apply on buffer
      __m128i loaded_result = _mm_shuffle_epi8(loaded_vector,loaded_transformation);
      _mm_store_si128(vector,loaded_result);

      uint8_t* result = (uint8_t*) vector;

      //save back in env
      for(size_t saved_qgram = filled_qgram; saved_qgram <= qgram_pos; saved_qgram--)
      {
        size_t curr_saved_qgram = qgram_pos - saved_qgram;
        
        //assert(curr_saved_qgram < constructed_env.size());

        uint8_t extracted_from_result[weight]{};
        for(uint8_t j = 0; j < weight; j++)
        {
          extracted_from_result[j] = result[(filled_qgram-saved_qgram)*weight+j];
        }
        
        constructed_env[curr_saved_qgram].set_qgram(extracted_from_result);
        /*std::cout << (int) curr_saved_qgram << std::endl;
        for(uint8_t i = 0; i < weight; i++){
          std::cout << (int)constructed_env[curr_saved_qgram].qgram[i] << '\t';
        }
        std::cout << std::endl;*/
      }
      std::cout << "Result\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)result[j] << '\t';
      }
      std::cout << std::endl;
    }
  }
#endif

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

  void env_group_get(
    const std::array<uint16_t,num_of_primary_env> sorted_code) {

    //std::array<std::vector<ScoreQgramcodePair>,num_of_primary_env> env_group{};
    
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get<qgram_length_arr[env_idx],
                                          env_threshold_arr[env_idx]>(
                                          sorted_code[env_idx]);
    });
    
    //return env_group;
  }

  std::array<EnvInfo,num_of_primary_env> get_env_info(const 
            std::array<uint16_t,num_of_primary_env>& encoded_qgram_arr) {
    std::array<EnvInfo,num_of_primary_env> env_info{};
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx){
      /*
      qgram_env_group[env_idx] = env_get<qgram_length_arr[env_idx],
                                          env_threshold_arr[env_idx]>(
                                          encoded_qgram_arr[env_idx]);
                                          */
      qgram_env_group[env_idx] = env_get_param<env_idx>(encoded_qgram_arr[env_idx]);
    });

    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      if(qgram_env_group[env_idx].empty()){
        return env_info;
      }
    }

    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      const auto qgram_env = qgram_env_group[env_idx];
      env_info[env_idx].max_score = qgram_env[0].score;

      int8_t curr_score = qgram_env[0].score;
      uint16_t start = 0;
      
      for(uint16_t idx = 1; idx < qgram_env.size(); idx++)
      {
        const auto elem = qgram_env[idx];
        if(elem.score < curr_score)
        {
          env_info[env_idx].threshold_arr.push_back(ScorePosition(start,idx-1));
          for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.score -1); j++)
          {
            env_info[env_idx].threshold_arr.push_back(ScorePosition(__UINT16_MAX__,__UINT16_MAX__));  
          }
          start = idx;
          curr_score = elem.score;
        }
      }
      env_info[env_idx].threshold_arr.push_back(ScorePosition(start,qgram_env.size()-1));
    });
    /*for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      const auto qgram_env = qgram_env_group[env_idx];
      env_info[env_idx].max_score = qgram_env[0].score;
      
      int8_t curr_score = qgram_env[0].score;
      uint16_t start = 0;
      
      for(uint16_t idx = 1; idx < qgram_env.size(); idx++){
        const auto elem = qgram_env[idx];
        if(elem.score < curr_score){
          env_info[env_idx].threshold_arr.push_back(ScorePosition(start,idx-1));
          for(uint16_t j = 0; j < static_cast<uint16_t>(curr_score - elem.score -1); j++){
            env_info[env_idx].threshold_arr.push_back(ScorePosition(__UINT16_MAX__,__UINT16_MAX__));  
          }
          start = idx;
          curr_score = elem.score;
        }
      }
      env_info[env_idx].threshold_arr.push_back(ScorePosition(start,qgram_env.size()-1));
    }*/
    return env_info;
  }

  std::vector<TempEnvElem<num_of_primary_env>> create_temp_env(std::array<EnvInfo,num_of_primary_env> env_info){
    std::vector<TempEnvElem<num_of_primary_env>> temp_env{};
    
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      if(env_info[env_idx].threshold_arr.empty()){
        return temp_env;
      }
    }
    size_t wheel[num_of_primary_env]{};
    int8_t loopid = num_of_primary_env - 1;
    
    int8_t max_score = 0;
    std::array<ScorePosition,num_of_primary_env> sp_arr{};
    for(size_t i = 0; i < num_of_primary_env;i++){
      max_score += env_info[i].max_score;
      sp_arr[i] = env_info[i].threshold_arr[0];
    }
    if(max_score >= threshold){
      temp_env.push_back(TempEnvElem<num_of_primary_env>(max_score,sp_arr));
    }
    else {
      return temp_env;
    }

    bool entry_valid;

    for (;;){
      wheel[loopid]++;
      if(wheel[loopid] == env_info[loopid].threshold_arr.size()){
        wheel[loopid] = 0;
        if(loopid == 0) break;
        loopid--;
      } 
      else {
        entry_valid = true;
        int8_t score = max_score;
        for(size_t i = 0; i < num_of_primary_env;i++){
          score -= wheel[i];
          if(env_info[i].threshold_arr[wheel[i]].start == __UINT16_MAX__) entry_valid = false;
          if(score < threshold){
            entry_valid = false;
            break;
          }
          sp_arr[i] = env_info[i].threshold_arr[wheel[i]];
          //sp_arr[i].end = env_info[i].threshold_arr[wheel[i]].end;
        }
        
        if(entry_valid){
          //insertion sort
          uint16_t idx = 0;
          while(idx < temp_env.size())
          {
            if(temp_env[idx].score < score)
            {
              break;
            }
            idx++;
          }
          TempEnvElem<num_of_primary_env> elem(score,sp_arr);
          if(idx == temp_env.size())
          {
            temp_env.push_back(elem);
          }
          else
          {
            temp_env.insert(temp_env.begin()+idx,elem);
          }
        }

        if(loopid < num_of_primary_env-1) loopid = num_of_primary_env-1;
      }
    }
/*
    for(size_t i = 0; i < temp_env.size(); i++){
      std::cout << (int) temp_env[i].score << '\t';
      for(uint8_t j = 0; j < num_of_primary_env; j++){
        std::cout << (int) temp_env[i].positions[j].start << '\t' << (int) temp_env[i].positions[j].end << '\t';
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
*/
    return temp_env;
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

  size_t get_composite_code(std::array<uint8_t,weight> qgram){
    size_t code = 0;
    for(uint8_t i = 0; i < weight; i++){
      code *= undefined_rank;
      code += qgram[i];
    }
    return code;
  }

  void add_in_curr_env(const std::vector<TempEnvElem<num_of_primary_env>> temp_env,
                      const std::array<uint8_t,weight> permutation, const bool sorted,
                      const size_t position, const bool with_simd){
    if(temp_env.empty()) return;

    std::array<uint16_t,num_of_primary_env> extracted_qgram_codes{};

    for(auto const& elem : temp_env)
    {
      const int8_t score = elem.score;
      const auto index_group = elem.positions;
      
      size_t wheel[num_of_primary_env]{};
      int8_t loopid = num_of_primary_env - 1;
      size_t starting_points[num_of_primary_env]{};

      //init wheels (starting point on qgram environments)
      for(uint8_t i = 0; i < num_of_primary_env; i++){
        starting_points[i] = index_group[i].start;
        wheel[i] = index_group[i].start;
      }

      //1st iteration
      for(uint8_t i = 0; i < num_of_primary_env; i++){
        extracted_qgram_codes[i] = qgram_env_group[i][wheel[i]].code;
      }
      auto reconstructed_qgram = reconstruct_qgram(extracted_qgram_codes);
      if(!sorted and !with_simd){
        reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
      }
      const OutputType env_elem{position,score,reconstructed_qgram};
      constructed_env.push_back(env_elem);

      for (;;)
      {
        wheel[loopid]++;
        if(wheel[loopid] > index_group[loopid].end)
        {
          wheel[loopid] = starting_points[loopid];
          if(loopid == 0) break;
          loopid--;
        }
        else {
          for(uint8_t i = 0; i < num_of_primary_env; i++){
            extracted_qgram_codes[i] = qgram_env_group[i][wheel[i]].code;
          }

          auto reconstructed_qgram = reconstruct_qgram(extracted_qgram_codes);
          if(!sorted and !with_simd){
            reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          }
          const OutputType env_elem{position,score,reconstructed_qgram};
          constructed_env.push_back(env_elem);

          if(loopid < num_of_primary_env-1) loopid = num_of_primary_env-1;
        }
      }
    }
#ifdef __SSSE3__
    if(!sorted and with_simd){
      reshuffle_with_simd(permutation,position);
    }
#endif
  }

  void direct_add(const std::array<uint16_t,num_of_primary_env>& qgram_codes,
                  const std::array<uint8_t,weight>& permutation,const bool sorted,
                  const size_t position, GTTL_UNUSED const bool with_simd) {
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<env_idx>(qgram_codes[env_idx]);
    });
    std::fill(std::begin(indexer),std::end(indexer),0);
    temp_env.clear();

    const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};

    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }
    if(score < threshold) return;

    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
    if(!sorted and !with_simd){
      reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    }
    const OutputType env_elem{position,score,reconstructed_qgram};
    constructed_env.push_back(env_elem);

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
        if(score >= threshold){
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          
          auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
          if(!sorted and !with_simd){
            reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          }

          const OutputType env_elem{position,score,reconstructed_qgram};
          
          //insertion_sort
          temp_env.insert(std::next(temp_env.begin(),indexer[max_possible_score-score]),env_elem);
          for(uint8_t i = max_possible_score-score; i < indexer.size(); i++){
            indexer[i]++;
          }
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=qgram_env_group[last_loop].size()-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }

        if(loopid < last_loop) loopid = last_loop;
      }
    }

    //constructed_env.insert(constructed_env.end(),std::begin(temp_env),std::end(temp_env));
    for(const auto elem : temp_env){
      constructed_env.push_back(elem);
    }
#ifdef __SSSE3__
    if(!sorted and with_simd){
      reshuffle_with_simd(permutation,position);
    }
#endif
  }

  void direct_add_list(const std::array<uint16_t,num_of_primary_env> qgram_codes,
                  const std::array<uint8_t,weight> permutation,const bool& sorted,
                  const size_t& position, GTTL_UNUSED const bool with_simd) {
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<env_idx>(qgram_codes[env_idx]);
    });
    temp_env.clear();
    for(uint8_t i = 0; i < max_possible_score-threshold+1; i++){
      indexer[i] = temp_env.end();
    }

    const auto last_loop = num_of_primary_env-1;
    //Init wheels and loopid
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env - 1;

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};

    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score += qgram_env_group[env_idx][wheel[env_idx]].score;
      unsorted_qgram_codes[env_idx] = qgram_env_group[env_idx][wheel[env_idx]].code;
    }
    if(score < threshold) return;

    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
    if(!sorted and !with_simd){
      reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    }
    const OutputType env_elem{position,score,reconstructed_qgram};
    temp_env.insert(temp_env.begin(),env_elem);
    indexer[max_possible_score-score]=temp_env.begin();

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
        if(score >= threshold){
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
          
          auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
          if(!sorted and !with_simd){
            reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          }

          const OutputType env_elem{position,score,reconstructed_qgram};
          const uint8_t indexer_idx = max_possible_score-score;

          //insertion_sort
          if(indexer[indexer_idx] == temp_env.end()){
            for(auto it = temp_env.begin();; it++){
              if(it == temp_env.end() or (*it).score < score){
                temp_env.insert(it,env_elem);
                indexer[indexer_idx] = std::prev(it);
                break;
              }
            }
          }
          else{
            temp_env.insert(std::next(indexer[indexer_idx]),env_elem);
            indexer[indexer_idx]++;
          }
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=qgram_env_group[last_loop].size()-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }

        if(loopid < last_loop) loopid = last_loop;
      }
    }

    //constructed_env.insert(constructed_env.end(),std::begin(temp_env),std::end(temp_env));
    for(const auto elem : temp_env){
      constructed_env.push_back(elem);
    }
#ifdef __SSSE3__
    if(!sorted and with_simd){
      reshuffle_with_simd(permutation,position);
    }
#endif
  }

  void single_thread_run(const size_t& start_idx,const size_t& end_idx,
                        const std::array<uint8_t,weight>& permutation,bool sorted,
                        const size_t& position, bool with_simd){
    size_t wheel[num_of_primary_env]{};
    uint8_t loopid = num_of_primary_env-1;
    const uint8_t last_loop = loopid;
    
    std::array<size_t,num_of_primary_env> start_arr{};
    start_arr[0] = start_idx;

    std::array<size_t,num_of_primary_env> end_arr{};
    end_arr[0] = end_idx;
    for(uint8_t env_idx = 1; env_idx < num_of_primary_env; env_idx++){
      end_arr[env_idx] = qgram_env_group[env_idx].size();
    }

    int8_t score = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes{};

    for(uint8_t i = 0; i < num_of_primary_env; i++){
      score += qgram_env_group[i][start_arr[i]].score;
      unsorted_qgram_codes[i] = qgram_env_group[i][start_arr[i]].code;
    }
    //1st loop
    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
    if(!sorted and !with_simd){
      reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    }
    const OutputType env_elem{position,score,reconstructed_qgram};
    if(true){
      std::lock_guard<std::mutex> guard{env_mutex};
#ifdef VECTOR      
      temp_env.insert(std::next(temp_env.begin(),indexer[max_possible_score-score]),env_elem);
      for(uint8_t i = max_possible_score-score; i < indexer.size(); i++){
        indexer[i]++;
      }
#endif
#ifdef LIST
      if(temp_env.empty()){
        temp_env.insert(temp_env.begin(),env_elem);
        indexer[max_possible_score-score]=temp_env.begin();
      } else {
        const uint8_t indexer_idx = max_possible_score-score;
        
        if(indexer[indexer_idx] == temp_env.end()){
          for(auto it = temp_env.begin();; it++){
            if(it == temp_env.end() or (*it).score < score){
              temp_env.insert(it,env_elem);
              indexer[indexer_idx] = std::prev(it);
              break;
            }
          }
        }
        else{
          temp_env.insert(std::next(indexer[indexer_idx]),env_elem);
          indexer[indexer_idx]++;
        }
      }
#endif
    }
    
    for(;;){
      score -= qgram_env_group[loopid][wheel[loopid]].score;
      wheel[loopid]++;
      
      if(wheel[loopid] >= end_arr[loopid]){
        wheel[loopid] = start_arr[loopid];
        score += qgram_env_group[loopid][start_arr[loopid]].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][start_arr[loopid]].code;
        if(loopid == 0) break;
        loopid--;
      }else{
        score += qgram_env_group[loopid][wheel[loopid]].score;
        unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
        
        if(score >= threshold){
          reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
          if(!sorted and !with_simd){
            reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
          }

          const OutputType env_elem{position,score,reconstructed_qgram};

          std::lock_guard<std::mutex> guard{env_mutex};
#ifdef VECTOR          
          temp_env.insert(std::next(temp_env.begin(),indexer[max_possible_score-score]),env_elem);
          for(uint8_t i = max_possible_score-score; i < indexer.size(); i++){
            indexer[i]++;
          }
#else
          const uint8_t indexer_idx = max_possible_score-score;
          
          if(indexer[indexer_idx] == temp_env.end()){
            for(auto it = temp_env.begin();; it++){
              if(it == temp_env.end() or (*it).score < score){
                temp_env.insert(it,env_elem);
                indexer[indexer_idx] = std::prev(it);
                break;
              }
            }
          }
          else{
            temp_env.insert(std::next(indexer[indexer_idx]),env_elem);
            indexer[indexer_idx]++;
          }
#endif
        }else{
          score -= qgram_env_group[last_loop][wheel[last_loop]].score;
          wheel[last_loop]=end_arr[loopid]-1;
          score += qgram_env_group[last_loop][wheel[last_loop]].score;
          unsorted_qgram_codes[last_loop] = qgram_env_group[last_loop][wheel[last_loop]].code;
        }
        if(loopid < last_loop) loopid = last_loop;
      }
    }
  }

  void direct_add_multithreaded(const std::array<uint16_t,num_of_primary_env> qgram_codes,
                  const std::array<uint8_t,weight>& permutation,bool sorted,
                  const size_t& position, bool with_simd) {
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<env_idx>(qgram_codes[env_idx]);
    });

    for(size_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      if(qgram_env_group[env_idx].empty()) return;
    }
    
    temp_env.clear();
#ifdef VECTOR
    std::fill(std::begin(indexer),std::end(indexer),0);
#else
    std::fill(std::begin(indexer),std::end(indexer),temp_env.end());
#endif
    
    //std::cout << (int)qgram_env_group[0][0].score << std::endl;
    int8_t max_rest_score = 0;
    for(uint8_t env_idx = 1; env_idx < num_of_primary_env; env_idx++){
      max_rest_score += qgram_env_group[env_idx][0].score;
      //std::cout << (int)qgram_env_group[env_idx][0].score << std::endl;
    }

    if(max_rest_score + qgram_env_group[0][0].score < threshold) return;

    if(num_of_primary_env >= 2){
      //Init threads
      std::array<std::thread,NUM_THREADS> threads;
      std::array<size_t,NUM_THREADS+1> thread_work_idx{};
      size_t end_in_first = 0;
      for(;end_in_first < qgram_env_group[0].size() and qgram_env_group[0][end_in_first].score >= threshold - max_rest_score; end_in_first++);
      
      const uint8_t total_thread_num = (end_in_first < NUM_THREADS) ? end_in_first : NUM_THREADS;
      
      thread_work_idx[total_thread_num] = end_in_first;
      thread_work_idx[0] = 0;

      const size_t max_thread_load = end_in_first/total_thread_num + bool(end_in_first % total_thread_num);
      size_t overcount = (total_thread_num-end_in_first%total_thread_num)*(end_in_first>total_thread_num)%total_thread_num;

      size_t curr_load = max_thread_load;
      for(size_t i = 1; i < total_thread_num; i++){
        if(overcount){
          overcount--;
          curr_load--;
        }
        assert(curr_load < end_in_first); 
        thread_work_idx[i] = curr_load;
        curr_load += max_thread_load;
      }

      for(uint8_t thread_idx = 0; thread_idx < total_thread_num; thread_idx++){
        const size_t start_idx = thread_work_idx[thread_idx];
        const size_t end_idx = thread_work_idx[thread_idx+1];

        threads[thread_idx] = std::thread(&CompositeEnvironment::single_thread_run,
                                          this,start_idx,end_idx,
                                          permutation,sorted,position,with_simd);
      }
      for(uint8_t thread_idx = total_thread_num-1; thread_idx < total_thread_num; thread_idx--){
        threads[thread_idx].join();
      }
    }

    //constructed_env.insert(constructed_env.end(),std::begin(temp_env),std::end(temp_env));
    for(const auto elem : temp_env){
      constructed_env.push_back(elem);
    }
#ifdef __SSSE3__
    if(!sorted and with_simd){
      reshuffle_with_simd(permutation,position);
    }
#endif
  }

  void direct_add_openmp(const std::array<uint16_t,num_of_primary_env> qgram_codes,
                  const std::array<uint8_t,weight> permutation,const bool& sorted,
                  const size_t& position, GTTL_UNUSED const bool with_simd) {
    //Init all primary environments
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      qgram_env_group[env_idx] = env_get_param<env_idx>(qgram_codes[env_idx]);
    });

    //Reset temp_env & indexer
    temp_env.clear();
#ifdef VECTOR
    std::fill(std::begin(indexer),std::end(indexer),0);
#else
    std::fill(std::begin(indexer),std::end(indexer),temp_env.end());
#endif
    const int8_t corrected_threshold = threshold - static_cast<int8_t>(local_score_correction);
    //std::cout << (int)(corrected_threshold < threshold) << '\t' << (int)corrected_threshold << '\t' << (int)threshold << '\t' << (int) local_score_correction << std::endl;

    int8_t score_1 = 0;
    std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes_1{};

    const auto last_loop = num_of_primary_env-1;
    //1st loop
    for(uint8_t env_idx = 0; env_idx < num_of_primary_env; env_idx++){
      score_1 += qgram_env_group[env_idx][0].score;
      unsorted_qgram_codes_1[env_idx] = qgram_env_group[env_idx][0].code;
    }
    const int8_t max_res_score = score_1-qgram_env_group[0][0].score;
    const int8_t threshold_1 = corrected_threshold-max_res_score;
    if(score_1 < threshold) return;
    size_t cutoff_idx = 0;
    for(; cutoff_idx < qgram_env_group[0].size();cutoff_idx++){
      if(qgram_env_group[0][cutoff_idx].score < threshold_1) break;
    }

    auto reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes_1);
    if(!sorted and !with_simd){
      reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
    }
    const OutputType env_elem{position,score_1,reconstructed_qgram};
    constructed_env.push_back(env_elem);

    #pragma omp parallel for
    for(size_t i = 0; i < cutoff_idx; i++){
      const auto elem_in_first = qgram_env_group[0][i];
      const int8_t score_in_first = elem_in_first.score;
      if(score_in_first < threshold_1) continue;
      
      const int8_t local_threshold = corrected_threshold - score_in_first;
      
      //Init wheels and loopid
      size_t wheel[num_of_primary_env]{};
      uint8_t loopid = num_of_primary_env - 1;

      std::array<uint16_t,num_of_primary_env> unsorted_qgram_codes = unsorted_qgram_codes_1;
      unsorted_qgram_codes[0] = elem_in_first.code;
      int8_t score = max_res_score;

      for(;;){
        score -= qgram_env_group[loopid][wheel[loopid]].score;
        wheel[loopid]++;
        
        if(wheel[loopid] >= qgram_env_group[loopid].size()){
          wheel[loopid] = 0;
          score += qgram_env_group[loopid][0].score;
          unsorted_qgram_codes[loopid] = qgram_env_group[loopid][0].code;

          if(loopid == 1) break;
          loopid--;
        }else{
          score += qgram_env_group[loopid][wheel[loopid]].score;
          if(score >= local_threshold){
            unsorted_qgram_codes[loopid] = qgram_env_group[loopid][wheel[loopid]].code;
            const int8_t true_score = score + score_in_first;
            
            reconstructed_qgram = reconstruct_qgram(unsorted_qgram_codes);
            if(!sorted and !with_simd){
              reconstructed_qgram = reshuffle_no_simd(reconstructed_qgram,permutation);
            }

            const OutputType env_elem{position,true_score,reconstructed_qgram};
            #pragma omp critical
            if(true){
              
              temp_env.insert(std::next(temp_env.begin(),indexer[max_possible_score-true_score]),env_elem);
              for(uint8_t i = max_possible_score-true_score; i < indexer.size(); i++){
                indexer[i]++;
              }
            
            }
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

    //constructed_env.insert(constructed_env.end(),std::begin(temp_env),std::end(temp_env));
    for(const auto elem : temp_env){
      constructed_env.push_back(elem);
    }
#ifdef __SSSE3__
    if(!sorted and with_simd){
      reshuffle_with_simd(permutation,position);
    }
#endif
  }

  void score_correction_set(const char* seq, const size_t seq_len, size_t position){
    local_score_correction = 0;
    for(uint8_t i = span-1; i < span; i-- && position++){
      if(!seed_bitset[i]) continue;

      size_t start,end;
      start = (position < CORRECTION_REGION_SIZE) ? 0 : (position - CORRECTION_REGION_SIZE);
      end = (position >= seq_len - CORRECTION_REGION_SIZE) ? (seq_len - 1) : (position + CORRECTION_REGION_SIZE);
      const size_t length = end - start; 
      
      float sum = 0;
      for(size_t seq_idx = start; seq_idx <= end; seq_idx++){
        if(seq_idx == position) continue;

        sum += sc.score_matrix[alpha.char_to_rank(seq[position])][alpha.char_to_rank(seq[seq_idx])];
      }

      local_score_correction -= sum/length;
      local_score_correction += background_correction[alpha.char_to_rank(seq[position])];
    }
  }

  public:
  CompositeEnvironment(){};

  void background_correction_set(const std::array<size_t,alpha_size+1>& target_distribution) {
    //calculate frequency of background data
    std::array<float,undefined_rank> frequency;
    size_t sum_count = 0;
    for(uint8_t i = 0; i < alpha_size; i++){
      sum_count += target_distribution[i];
    }
    for(uint8_t i = 0; i < undefined_rank; i++){
      frequency[i] = target_distribution[i]/sum_count;
    }

    //calculate background correction 
    for(uint8_t i = 0; i < undefined_rank; i++){
      float sum_char = 0;
      for(uint8_t j = 0; j < undefined_rank; j++){
        sum_char += frequency[j] * sc.score_matrix[i][j];
      }
      background_correction[i] = sum_char/undefined_rank;
    }
  }

  void process_seed(const char* seq,const size_t seq_len,const size_t position,
                    GTTL_UNUSED const bool with_simd) {
    //RunTimeClass rt{};
    //preprocess seed
    const auto encode_info = spaced_seed_encoder.encode(seq);
    //rt.show("Finished encode spaced seed");
    score_correction_set(seq,seq_len,position);
    /*for(size_t i = 0; i < num_of_primary_env; i++){
      std::cout << (int)encode_info.codes[i] << '\t';
    }
    std::cout << std::endl;*/
#ifdef __SSSE3__
    direct_add_multithreaded(encode_info.codes,encode_info.permutation,encode_info.sorted,
              position,with_simd);
#else
    direct_add_multithreaded(encode_info.codes,encode_info.permutation,encode_info.sorted,
              position,false);
#endif
/*
    //threshold_arr
    const auto env_info = get_env_info(encode_info.codes);
    //rt.show("Finished get env info");
    
    //temp_env
    const auto temp_env = create_temp_env(env_info);
    //rt.show("Finished creating temp env");
    //add into curr env
#ifdef __SSSE3__
    add_in_curr_env(temp_env,encode_info.permutation,encode_info.sorted,
                    position,with_simd);
#else
    add_in_curr_env(temp_env,encode_info.permutation,encode_info.sorted,
                    position,false);
#endif
    //rt.show("Finished adding in curr env");*/
    
  }

  constexpr size_t num_chars() const {
    return undefined_rank;
  } 

  size_t size() const
  {
    return constructed_env.size();
  }

  OutputType elem_get(const size_t idx) const
  {
    return constructed_env[idx];
  }

  void reset() {
    constructed_env.clear();
  }

  constexpr uint8_t span_get() const {
    return span;
  }

  constexpr uint8_t weight_get() const {
    return weight;
  }
};

#endif
