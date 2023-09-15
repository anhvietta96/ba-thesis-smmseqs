#ifndef MMSEQS2_IT_HPP
#define MMSEQS2_IT_HPP

#include "sequences/literate_multiseq.hpp"
#include "filter/composite_env.hpp"
#include "filter/InvIntHash.hpp"

//Problem: data structure (associative array?)
//vector basetype in hash unit
//some bugs in hash func

struct MMseqs2Hit {
  size_t target_seq_num;
  size_t target_seq_pos;
  size_t query_seq_num;
  size_t query_seq_pos;
  size_t code;
  int8_t score;

  MMseqs2Hit(const size_t _target_seq_num,const size_t _target_seq_pos,
              const size_t _query_seq_num,const size_t _query_seq_pos,
              const size_t _code,const int8_t _score):
              target_seq_num(_target_seq_num),target_seq_pos(_target_seq_pos),
              query_seq_num(_query_seq_num),query_seq_pos(_query_seq_pos),
              code(_code),score(_score){};
};

template<class ScoreClass,template<class,size_t> class HashFunc,const size_t seed>
class MMseqs2 {
  CompositeEnvironment<ScoreClass,seed> env_constructor{};
  Multiseq_Hash<ScoreClass,HashFunc,seed> multiseq_hash{};
  std::vector<MMseqs2Hit> hits_data;
  
  void create_target_data(const GttlMultiseq* target){
    multiseq_hash.hash(target);
  }

  void iterate_over_query(const GttlMultiseq* query,const bool with_simd){
    const size_t total_seq_num = query->sequences_number_get();
    const auto span = env_constructor.span_get();
    const auto num_chars = env_constructor.num_chars();
    for(size_t seqnum = 0; seqnum < total_seq_num; seqnum++)
    {
      const char* curr_seq = query->sequence_ptr_get(seqnum);
      const size_t seq_len = query->sequence_length_get(seqnum);
      if(seq_len >= span)
      {
        for(size_t i = 0; i < seq_len - span + 1; i++)
        {
          env_constructor.process_seed(curr_seq+i,seq_len,i+1,with_simd);
        }
        for(size_t i = 0; i < env_constructor.size(); i++)
        {
          const auto elem = env_constructor.elem_get(i);
          const auto num_elem = elem.to_num_type(num_chars);
          find_hits(num_elem,seqnum);
        }
      }
    }
  }

  void find_hits(const NumEnvElem& query_hit, const size_t query_seq_num){
    const auto packer = multiseq_hash.packer_get();
    for(const auto& target_elem: multiseq_hash.bytes_unit_vec_get()){
      const auto target_elem_code = target_elem.template decode_at<2>(*packer);
      if(query_hit.code == target_elem_code){
        const auto target_seq_num = target_elem.template decode_at<0>(*packer);
        const auto target_seq_pos = target_elem.template decode_at<1>(*packer);
        const auto query_seq_pos = query_hit.position;
        const auto score = query_hit.score;
        hits_data.push_back(MMseqs2Hit(target_seq_num,target_seq_pos,
                                      query_seq_num,query_seq_pos,
                                      query_hit.code,score));
      }
    }
  }

  public:
  MMseqs2(const GttlMultiseq* query, const GttlMultiseq* target,
          const bool with_simd){
    constexpr const ScoreClass sc{};
    constexpr const auto char_spec = sc.character_spec;
    constexpr const auto undefined_rank = sc.num_of_chars;
    const LiterateMultiseq<char_spec,undefined_rank> target_literate{*target};
    const auto rank_dist = target_literate.rank_dist_get();

    env_constructor.background_correction_set(rank_dist);
    
    create_target_data(target);

    iterate_over_query(query,with_simd);
  };

  size_t size() const {
    return hits_data.size();
  }

  const MMseqs2Hit& hit_get(const size_t idx) const {
    return hits_data[idx];
  }
};
#endif