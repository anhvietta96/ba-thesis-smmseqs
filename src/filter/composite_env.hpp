#include "filter/qgram_environment.hpp"
#include "utilities/constexpr_for.hpp"
//#include "filter/multiset_code.hpp"

#include <xmmintrin.h>
#include <x86intrin.h>
#include "utilities/runtime_class.hpp"

#include <iostream>
#include <tuple>

static constexpr const int8_t predef_threshold_arr[5] = { 22, 27, 35, 40, 50};

template<const uint8_t num_of_primary_env,const uint8_t weight>
constexpr std::array<uint8_t,num_of_primary_env> create_qgram_length_arr(){
  std::array<uint8_t,num_of_primary_env> qgram_length_arr{};
  for(uint8_t i = 0; i < num_of_primary_env; i++){
    qgram_length_arr[i] = 3;
  }

  constexpr const uint8_t overcount = (num_of_primary_env * 3) - weight;
  for(uint8_t i = 0; i < overcount; i++){
    qgram_length_arr[i]--;
  }

  return qgram_length_arr;
}

template<const uint8_t num_of_primary_env,const uint8_t weight>
constexpr std::array<int8_t,num_of_primary_env> create_env_threshold_arr(const std::array<uint8_t,num_of_primary_env>& qgram_length_arr,const int8_t& threshold,const int8_t& highest_score){
  std::array<int8_t,num_of_primary_env> env_threshold_arr{};
  for(uint8_t i = 0; i < num_of_primary_env; i++){
    env_threshold_arr[i] = threshold - (weight - qgram_length_arr[i]) * highest_score;
  }
  return env_threshold_arr;
}

template<class ScoreClass,const int8_t threshold,const uint8_t weight>
class PrimaryEnvGroup {
  static constexpr const ScoreClass sc{};
  static constexpr const uint8_t num_of_primary_env = weight % 3 ? (weight / 3 + 1) : 
                                                      weight / 3;
  static constexpr const std::array<uint8_t,num_of_primary_env> qgram_length_arr = create_qgram_length_arr<num_of_primary_env,weight>();
  static constexpr const std::array<int8_t,num_of_primary_env> env_threshold_arr = create_env_threshold_arr<num_of_primary_env,weight>(qgram_length_arr,threshold,sc.highest_score);
/*
  constexpr std::array<uint8_t,num_of_primary_env> create_qgram_length_arr(){
    std::array<uint8_t,num_of_primary_env> qgram_length_arr{};
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      qgram_length_arr[i] = 3;
    }

    constexpr const uint8_t decrement = (num_of_primary_env * 3) - weight;
    for(uint8_t i = 0; i < decrement; i--){
      qgram_length_arr[i]--;
    }

    return qgram_length_arr;
  }

  constexpr std::array<int8_t,num_of_primary_env> create_env_threshold_arr(){
    std::array<int8_t,num_of_primary_env> env_threshold_arr{};
    for(uint8_t i = 0; i < num_of_primary_env; i++){
      env_threshold_arr[i] = threshold - (weight - qgram_length_arr) * sc.highest_score;
    }
    return env_threshold_arr;
  }
*/
  public:
  std::array<std::vector<ScoreQgramcodePair>,num_of_primary_env> env_group_get(
    const std::array<uint16_t,num_of_primary_env> sorted_code) const {

    std::array<std::vector<ScoreQgramcodePair>,num_of_primary_env> env_group{};
    
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      const QgramEnvironment<ScoreClass,qgram_length_arr[env_idx],
                    env_threshold_arr[env_idx]> env{};
      env_group[env_idx] = env.qgram_env(sorted_code[env_idx]);
    });
    
    return env_group;
  }

  constexpr uint8_t num_of_primary_get() const {
    return num_of_primary_env;
  }

  constexpr std::array<uint8_t,num_of_primary_env> qgram_length_arr_get() const {
    return qgram_length_arr;
  }
};

template<const uint8_t qgram_length>
struct LocalEnvElem {
  size_t position;
  int8_t score;
  uint8_t qgram[qgram_length];

  LocalEnvElem(const size_t _position, const int8_t _score,
    const std::array<uint8_t,qgram_length> _qgram){
    position = _position;
    score = _score;
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = _qgram[i];
    }
  };

  void set_qgram(const uint8_t* ref_qgram)
  {
    for(uint8_t i = 0; i < qgram_length; i++)
    {
      qgram[i] = ref_qgram[i];
    }
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

template<class ScoreClass, const uint8_t weight>
class CompositeEnvironment {
  //ScoreClass
  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  
  //PrimaryEnvGroup
  static constexpr const int8_t threshold = (weight > 8 || weight < 4) ? predef_threshold_arr[4] : predef_threshold_arr[weight-4];
  static constexpr const PrimaryEnvGroup<ScoreClass,threshold,weight> primary_env{};
  static constexpr const auto num_of_primary_env = primary_env.num_of_primary_get();
  static constexpr const auto qgram_length_arr = primary_env.qgram_length_arr_get();

  //Sorted/UnsortedQmer
  /*
  static constexpr const SortedQmer<char_spec,undefined_rank,2> sorted_q_2{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,2> unsorted_q_2{};
  static constexpr const SortedQmer<char_spec,undefined_rank,3> sorted_q_3{};
  static constexpr const UnsortedQmer<char_spec,undefined_rank,3> unsorted_q_3{};
  */

  //main_env
  std::vector<LocalEnvElem<weight>> constructed_env{};

  std::array<uint8_t,weight> process_qgram(std::array<uint8_t,weight> reconstructed_qgram,
  const std::array<uint8_t,weight> permutation) const
  {
    std::array<uint8_t,weight> transformed_qgram{};
    for(uint8_t idx = 0; idx < weight; idx++)
    {
      transformed_qgram[permutation[idx]] = reconstructed_qgram[idx];
    }
    return transformed_qgram;
  }

  void reshuffle_with_simd(const std::array<uint8_t,weight> permutation)
  {
    size_t start_idx = constructed_env.size()-1;
    const auto position = constructed_env[start_idx].position;
    
    while(start_idx < constructed_env.size())
    {
      if(constructed_env[start_idx].position != position) break;
      start_idx--;
    }

    const size_t num_of_qgrams = constructed_env.size()-1-start_idx;
    if(num_of_qgrams == 0) return;

    constexpr const uint8_t qgram_per_simd_vector = 16 / weight;

    const size_t num_of_simd_operations = (num_of_qgrams) ? ((num_of_qgrams-1) / qgram_per_simd_vector + 1) : 0;
    size_t qgram_pos = start_idx+1;
    size_t curr_qgram = 0;
    //std::cout << (int)position << '\t'<<(int) i << '\t' << (int) num_of_qgrams << '\t' << (int) num_of_simd_operations << '\t' << (int) qgram_pos << std::endl;
    
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
    /*
    std::cout << "Permutation\t";
    for(uint8_t j = 0; j < 16; j++)
    {
      std::cout << (int)transformation[j] << '\t';
    }
    std::cout << std::endl;
*/
    __m128i *transform_m128i = (__m128i *) transformation;
    __m128i loaded_transformation = _mm_load_si128(transform_m128i);
    
    for (size_t simd_op_no = 0; simd_op_no < num_of_simd_operations; simd_op_no++)
    {
      uint8_t buf[16]{};
      uint8_t filled_qgram = 0;
      while(curr_qgram < num_of_qgrams and filled_qgram < qgram_per_simd_vector)
      {
        //fill buffer
        auto qgram = constructed_env[qgram_pos].qgram;
        for(uint8_t char_idx = 0; char_idx < weight; char_idx++)
        {
          buf[filled_qgram*weight+char_idx] = qgram[char_idx];
        }
        //std::cout << (int)qgram_pos << std::endl;
        filled_qgram++;
        curr_qgram++;
        qgram_pos++;
      }/*
      std::cout << "Buffer\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)buf[j] << '\t';
      }
      std::cout << std::endl;*/

      __m128i *vector = (__m128i *) buf;
      __m128i loaded_vector = _mm_load_si128(vector);

      //apply on buffer
      __m128i loaded_result = _mm_shuffle_epi8(loaded_vector,loaded_transformation);
      _mm_store_si128(vector,loaded_result);

      uint8_t* result = (uint8_t*) vector;
      //std::cout << (int) filled_qgram << std::endl;
      //save back in env
      for(uint8_t saved_qgram = filled_qgram; saved_qgram > 0; saved_qgram--)
      {
        uint8_t curr_saved_qgram = qgram_pos - saved_qgram;
        //std::cout << '2' << '\t' << (int) curr_saved_qgram << std::endl;
        assert(curr_saved_qgram < constructed_env.size() and curr_saved_qgram >= 0);

        uint8_t extracted_from_result[weight]{};
        for(uint8_t j = 0; j < weight; j++)
        {
          extracted_from_result[j] = result[(filled_qgram-saved_qgram)*weight+j];
        }
        
        constructed_env[curr_saved_qgram].set_qgram(extracted_from_result);
      }/*
      std::cout << "Result\t";
      for(uint8_t j = 0; j < 16; j++)
      {
        std::cout << (int)result[j] << '\t';
      }
      std::cout << std::endl;*/
    }
  }

  std::array<uint16_t,num_of_primary_env> 
  extract_code(const size_t unsorted_qgram_code) const {
    std::array<uint16_t,num_of_primary_env> encoded_qgram_arr{};
    
    size_t reserve_code = unsorted_qgram_code;
    uint8_t back = weight;

    constexpr_for<0,num_of_primary_env-1,1>([&] (auto env_idx)
    {
      back -= qgram_length_arr[env_idx];
      const uint16_t unsorted_code = reserve_code / power_1(undefined_rank,back);
      constexpr const SortedQmer<char_spec,undefined_rank,qgram_length_arr[env_idx]> sorted_q{};
      encoded_qgram_arr[env_idx] = sorted_q.sorted_code_get(unsorted_code);
      reserve_code = reserve_code % power_1(undefined_rank,back);
    });

    constexpr const SortedQmer<char_spec,undefined_rank,qgram_length_arr[num_of_primary_env-1]> sorted_q{};
    encoded_qgram_arr[num_of_primary_env-1] = sorted_q.sorted_code_get(reserve_code);
    
    return encoded_qgram_arr;
  }

  std::array<EnvInfo,num_of_primary_env> get_env_info(
      const std::array<uint16_t,num_of_primary_env> encoded_qgram_arr) const {
    std::array<EnvInfo,num_of_primary_env> env_info{};
    const auto qgram_env_group = primary_env.env_group_get(encoded_qgram_arr);
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      env_info[env_idx].max_score = qgram_env_group[env_idx][0].score;
      const auto qgram_env = qgram_env_group[env_idx];

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
    return env_info;
  }

  std::vector<TempEnvElem<num_of_primary_env>> create_temp_env(std::array<EnvInfo,num_of_primary_env> env_info){
    std::vector<TempEnvElem<num_of_primary_env>> temp_env{};
    
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

  std::array<uint8_t,weight> reconstruct_qgram(std::array<uint16_t,num_of_primary_env> qgram_codes){
    std::array<uint8_t,weight> reconstructed_qgram{};
    uint8_t i = 0;
    constexpr_for<0,num_of_primary_env,1>([&] (auto env_idx)
    {
      constexpr const uint8_t qgram_length = qgram_length_arr[env_idx];
      constexpr const UnsortedQmer<char_spec,undefined_rank,qgram_length> unsorted_q{};
      const auto subqgram = unsorted_q.qgram_get(qgram_codes[env_idx]);
      constexpr_for<0,qgram_length,1>([&] (auto qgram_idx)
      {
        reconstructed_qgram[i] = subqgram[qgram_idx];
        i++;
      });
    });

    return reconstructed_qgram;
  }

  void add_in_curr_env(const std::array<uint16_t,num_of_primary_env> encoded_qgram_arr,
                      const std::vector<TempEnvElem<num_of_primary_env>> temp_env,
                      const std::array<uint8_t,weight> permutation, const bool sorted,
                      const size_t position, const bool with_simd){
    
    if(temp_env.empty()) return;
    const auto qgram_env_group = primary_env.env_group_get(encoded_qgram_arr);

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
      if(!sorted and !with_simd) reconstructed_qgram = process_qgram(reconstructed_qgram,permutation);
      const LocalEnvElem<weight> env_elem{position,score,reconstructed_qgram};
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
          if(!sorted and !with_simd) reconstructed_qgram = process_qgram(reconstructed_qgram,permutation);
          const LocalEnvElem<weight> env_elem{position,score,reconstructed_qgram};
          constructed_env.push_back(env_elem);

          if(loopid < num_of_primary_env-1) loopid = num_of_primary_env-1;
        }
      }
    }

    if(!sorted and with_simd){
      reshuffle_with_simd(permutation);
    }
  }

  public:
  void process_seed(const size_t unsorted_qgram_code,const std::array<uint8_t,weight> permutation, const bool sorted,
                    const size_t position, const bool with_simd) {
    RunTimeClass rt{};
    //extract code
    const auto encoded_qgram_arr = extract_code(unsorted_qgram_code);
    rt.show("Finished extracting code");
    //threshold_arr
    const auto env_info = get_env_info(encoded_qgram_arr);
    rt.show("Finished get env info");
    //temp_env
    const auto temp_env = create_temp_env(env_info);
    rt.show("Finished creating temp env");
    //add into curr env
    add_in_curr_env(encoded_qgram_arr,temp_env,permutation,sorted,position,with_simd);
    rt.show("Finished adding in curr env");
  }

  size_t size() const
  {
    return constructed_env.size();
  }

  LocalEnvElem<weight> elem_get(const size_t idx) const
  {
    return constructed_env[idx];
  }
};