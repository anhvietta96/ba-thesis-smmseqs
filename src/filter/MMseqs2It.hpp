#ifndef MMSEQS2_IT_HPP
#define MMSEQS2_IT_HPP

#include "sequences/literate_multiseq.hpp"
#include "filter/byte_composite_env.hpp"
#include "filter/InvIntHash.hpp"
#include "utilities/unused.hpp"
#include "utilities/runtime_class.hpp"
#include "filter/distribution.hpp"

#define TIME

template<class ScoreClass,template<class,size_t> class HashFunc,const size_t seed>
class MMseqs2 {
  static constexpr const uint8_t min_unit_size = 8;
  static constexpr const uint8_t max_unit_size = 9;

  static constexpr const Multiseq_Hash<ScoreClass,HashFunc,seed> multiseq_hash{};

  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  
  template<const uint8_t sizeof_query_unit, const uint8_t sizeof_target_unit,
    const uint8_t sizeof_match_unit>
  void find_hits(const std::vector<BytesUnit<sizeof_query_unit,3>>& query_vec,
    const std::vector<BytesUnit<sizeof_target_unit,3>>& target_vec,
    std::vector<BytesUnit<sizeof_match_unit,4>>& match_vec,
    const GttlBitPacker<sizeof_query_unit,3>& query_packer, 
    const GttlBitPacker<sizeof_target_unit,3>& target_packer,
    const GttlBitPacker<sizeof_match_unit,4>& match_packer) const {
    size_t target_idx = 0, query_idx = 0;
    while(target_idx < target_vec.size() and query_idx < query_vec.size()){
      const auto target_code = target_vec[target_idx].template decode_at<0>(target_packer);
      const auto query_code = query_vec[query_idx].template decode_at<0>(query_packer);
      
      if(target_code < query_code){
        //skip ahead on target
        target_idx++;
        while(target_idx < target_vec.size()){
          if(target_vec[target_idx].template decode_at<0>(target_packer) >= query_code) break;
          target_idx++;
        }
      
      } else if(target_code > query_code){
        //skip ahead on query
        query_idx++;
        while(query_idx < query_vec.size()){
          if(query_vec[query_idx].template decode_at<0>(query_packer) >= target_code) break;
          query_idx++;
        }
      
      } else {
        //collect matches
        size_t target_end,query_end;
        for(target_end = target_idx+1; target_end < target_vec.size(); target_end++){
          if(target_vec[target_end].template decode_at<0>(target_packer) != target_code) break; 
        }
        for(query_end = query_idx+1; query_end < query_vec.size(); query_end++){
          if(query_vec[query_end].template decode_at<0>(query_packer) != query_code) break; 
        }

        for(size_t query_tmp_idx = query_idx; query_tmp_idx < query_end; query_tmp_idx++){
          const auto query_seqnum = query_vec[query_tmp_idx].template decode_at<1>(query_packer);
          const auto query_seqpos = query_vec[query_tmp_idx].template decode_at<2>(query_packer);
          
          for(size_t target_tmp_idx = target_idx; target_tmp_idx < target_end; target_tmp_idx++){
            const auto target_seqnum = target_vec[target_tmp_idx].template decode_at<1>(target_packer);
            const auto target_seqpos = target_vec[target_tmp_idx].template decode_at<2>(target_packer);
            
            match_vec.emplace_back(match_packer,std::array<uint64_t,4>{static_cast<uint64_t>(target_seqnum),
                                                                      static_cast<uint64_t>(query_seqnum),
                                                                      static_cast<uint64_t>(target_seqpos),
                                                                      static_cast<uint64_t>(query_seqpos)});
          }
        }
        target_idx = target_end;
        query_idx = query_end;
      }
    }
  }

  template<const uint8_t sizeof_query_unit, const uint8_t sizeof_target_unit,
    const uint8_t sizeof_match_unit>
  void query_all_vs_all(GttlMultiseq* query, GttlMultiseq* target,const double sensitivity,
                        GTTL_UNUSED const bool with_simd,const bool short_header,const bool show, 
                        const bool mmseqs, const bool ctxsens,const bool correct,const double correct_ratio,
                        const size_t num_threads) const {
    const size_t target_seq_len_bits = target->sequences_length_bits_get();
    const size_t target_seq_num_bits = target->sequences_number_bits_get();
    const size_t query_seq_len_bits = query->sequences_length_bits_get();
    const size_t query_seq_num_bits = query->sequences_number_bits_get();
    const size_t hashbits = multiseq_hash.hashbits_get();

    LiterateMultiseq<char_spec,undefined_rank> literate_target{target};
    literate_target.perform_sequence_encoding();
    LiterateMultiseq<char_spec,undefined_rank> literate_query{query};
    literate_query.perform_sequence_encoding();
    
    std::vector<BytesUnit<sizeof_query_unit,3>> query_hash_data;
    std::vector<BytesUnit<sizeof_target_unit,3>> target_hash_data;
    std::vector<BytesUnit<sizeof_match_unit,4>> matches;
    const GttlBitPacker<sizeof_query_unit,3> query_packer{
                                    {static_cast<int>(hashbits),
                                    static_cast<int>(query_seq_num_bits),
                                    static_cast<int>(query_seq_len_bits)}};
    const GttlBitPacker<sizeof_target_unit,3> target_packer{
                                    {static_cast<int>(hashbits),
                                    static_cast<int>(target_seq_num_bits),
                                    static_cast<int>(target_seq_len_bits)}};
    const GttlBitPacker<sizeof_match_unit,4> match_packer{{
                                        {static_cast<int>(target_seq_num_bits),
                                        static_cast<int>(query_seq_num_bits),
                                        static_cast<int>(target_seq_len_bits),
                                        static_cast<int>(query_seq_len_bits)}}};
#ifdef TIME    
    RunTimeClass rt{};
    multiseq_hash.template hash<sizeof_target_unit>(target,target_hash_data,target_packer);
    rt.show("Target hashing finished");
    multiseq_hash.template sort<sizeof_target_unit>(target_hash_data,hashbits);
    rt.show("Target sorting finished");
    constexpr const uint8_t weight = multiseq_hash.weight_get();
    const BGDistribution<Blosum62,weight> distribution{};
    int64_t threshold;
    if(ctxsens){
      threshold = distribution.template context_sensitive_threshold_get<sizeof_target_unit>(target_hash_data,target_packer,sensitivity);
    } else {
      threshold = distribution.custom_threshold_get2(literate_target.rank_dist_get(),sensitivity);
    }
    std::cout << "Threshold " << (int) threshold << std::endl;
    rt.show("Evaluated threshold");
    const BytesCompositeEnvironment2<ScoreClass,seed> env_constructor{literate_target.rank_dist_get(),threshold};
    //env_constructor.set_background_data();
    rt.show("Prepare Query Input");
    env_constructor.template process_openmp<sizeof_query_unit>(query,query_hash_data,query_packer,mmseqs,with_simd,correct,correct_ratio,num_threads);
    rt.show("Constructed Cartesian products");
    std::cout << "Total generated qgrams: " << query_hash_data.size() << std::endl;
    /*for(size_t i = 0; i < query_hash_data.size(); i++){
      std::cout << query_hash_data[i].template decode_at<0>(query_packer) << '\t' << query_hash_data[i].template decode_at<2>(query_packer) << std::endl;
    }*/
    env_constructor.template sort<sizeof_query_unit>(query_hash_data,hashbits);
    rt.show("Query sorting finished");
    find_hits<sizeof_query_unit,sizeof_target_unit,sizeof_match_unit>(query_hash_data,target_hash_data,matches,query_packer,target_packer,match_packer);
    rt.show("Merging query - target completed");
    hit_sort<sizeof_match_unit>(matches,target_seq_num_bits+query_seq_num_bits);
    rt.show("Match data sorted");
    std::cout << "Matches found: " << matches.size() << std::endl;
#else
    multiseq_hash.template hash<sizeof_target_unit>(target,target_hash_data,target_packer);
    
    multiseq_hash.template sort<sizeof_target_unit>(target_hash_data,hashbits);
    
    constexpr const uint8_t weight = multiseq_hash.weight_get();
    const BGDistribution<Blosum62,weight> distribution{};
    //const auto threshold = distribution.template context_sensitive_threshold_get<sizeof_target_unit>(target_hash_data,target_packer,sensitivity);
    const auto threshold = distribution.custom_threshold_get2(literate_target.rank_dist_get(),sensitivity);
    //std::cout << "Threshold " << (int) threshold << std::endl;
    const BytesCompositeEnvironment2<ScoreClass,seed> env_constructor{literate_target.rank_dist_get(),threshold};
    //env_constructor.set_background_data(literate_target.rank_dist_get(),threshold);
    
    env_constructor.template process_pthread<sizeof_query_unit>(query,query_hash_data,query_packer,mmseqs,with_simd,correct,correct_ratio,num_threads);
    
    std::cout << "Total generated qgrams: " << query_hash_data.size() << std::endl;
    env_constructor.template sort<sizeof_query_unit>(query_hash_data,hashbits);
    
    find_hits<sizeof_query_unit,sizeof_target_unit,sizeof_match_unit>(query_hash_data,target_hash_data,matches,query_packer,target_packer,match_packer);
    
    hit_sort<sizeof_match_unit>(matches,target_seq_num_bits+query_seq_num_bits);
    
    std::cout << "Matches found: " << matches.size() << std::endl;
#endif
    if(show){
      if(!short_header){
        std::cout << "#target_seq_num" << '\t' << "query_seq_num" << '\t' << 
        "target_seq_pos" << '\t' << "query_seq_pos" << '\t' << std::endl;
        for(size_t i = 0; i < matches.size(); i++){
          const auto hit = matches[i];
          std::cout << (int) hit.template decode_at<0>(match_packer) << '\t' 
          << (int) hit.template decode_at<1>(match_packer) << '\t' << 
          (int) hit.template decode_at<2>(match_packer) << '\t' 
          << (int) hit.template decode_at<3>(match_packer) <<  std::endl;
        }
      } else {
        std::cout << "#target_header" << '\t' << "query_header" << '\t' << 
        "target_seq_pos" << '\t' << "query_seq_pos" << '\t' << std::endl;
        for(size_t i = 0; i < matches.size(); i++){
          const auto hit = matches[i];
          const auto target_seq_num = hit.template decode_at<0>(match_packer);
          size_t target_sh_offset, target_sh_len;
          std::tie(target_sh_offset,target_sh_len) = target->short_header_get(target_seq_num);
          const std::string_view target_seq_header = target->header_get(target_seq_num);

          const auto query_seq_num = hit.template decode_at<1>(match_packer);
          size_t query_sh_offset, query_sh_len;
          std::tie(query_sh_offset,query_sh_len) = query->short_header_get(query_seq_num);
          const std::string_view query_seq_header = query->header_get(query_seq_num);

          std::cout << target_seq_header.substr(target_sh_offset,target_sh_len) << '\t' 
          << query_seq_header.substr(query_sh_offset,query_sh_len) << '\t' << 
          (int) hit.template decode_at<2>(match_packer) << '\t' 
          << (int) hit.template decode_at<3>(match_packer) <<  std::endl;
        }
      }
    }
  }

  uint8_t sizeof_unit_get(const size_t total_bits) const {
    if(total_bits <= 64) return 8;
    if(total_bits % 8 == 0) return total_bits/8;
    return total_bits/8+1;
  }

  template<const size_t sizeof_match_unit>
  void hit_sort(std::vector<BytesUnit<sizeof_match_unit,4>>& container, const size_t sort_bits) const {
    if constexpr (sizeof_match_unit == 8){
      ska_lsb_radix_sort<uint64_t>(sort_bits,
                              reinterpret_cast<uint64_t *>
                              (container.data()),
                              container.size());
    }
    else {
      ska_large_lsb_small_radix_sort(sizeof_match_unit,sort_bits,reinterpret_cast<uint8_t *>(
                                  container.data()),container.size(),
                                  false);
    }
  }

  public:
  MMseqs2(GttlMultiseq* query, GttlMultiseq* target,const double sensitivity,
          const bool with_simd,const bool short_header,const bool show, const bool mmseqs, const bool ctxsens,
          const bool correct, const double correct_ratio, const size_t num_threads){
    
    const size_t target_seq_len_bits = target->sequences_length_bits_get();
    const size_t target_seq_num_bits = target->sequences_number_bits_get();
    const size_t query_seq_len_bits = query->sequences_length_bits_get();
    const size_t query_seq_num_bits = query->sequences_number_bits_get();
    const size_t hashbits = multiseq_hash.hashbits_get();

    const auto sizeof_target_unit = sizeof_unit_get(hashbits+target_seq_num_bits+target_seq_len_bits);
    const auto sizeof_query_unit = sizeof_unit_get(hashbits+query_seq_num_bits+query_seq_len_bits);
    const auto sizeof_match_unit = sizeof_unit_get(query_seq_num_bits+query_seq_len_bits+target_seq_num_bits+target_seq_len_bits);

    //std::cout << (int) sizeof_target_unit << '\t' << (int) sizeof_query_unit << '\t' << (int) sizeof_match_unit << std::endl;
    /*std::cout << (int) hashbits+target_seq_num_bits+target_seq_len_bits << '\t' 
    << (int) hashbits+query_seq_num_bits+query_seq_len_bits << '\t' 
    << (int) query_seq_num_bits+query_seq_len_bits+target_seq_num_bits+target_seq_len_bits << std::endl;*/
    
    std::cout << "Hashbits: " << (int) hashbits << std::endl;
    std::cout << "Target_seqnum_bits: " << (int) target_seq_num_bits << std::endl;
    std::cout << "Target_seqpos_bits: " << (int) target_seq_len_bits << std::endl; 
    std::cout << "Query_seqnum_bits: " << (int) query_seq_num_bits << std::endl;
    std::cout << "Query_seqpos_bits: " << (int) query_seq_len_bits << std::endl;

    if(sizeof_match_unit > max_unit_size or sizeof_query_unit > max_unit_size
    or sizeof_target_unit > max_unit_size){
      std::cerr << " Over max allowed size" << std::endl;
      return;
    }

    constexpr_for<min_unit_size,max_unit_size+1,1>([&] (auto constexpr_target_size){
      constexpr_for<min_unit_size,max_unit_size+1,1>([&] (auto constexpr_query_size){
        constexpr_for<min_unit_size,max_unit_size+1,1>([&] (auto constexpr_match_size){
          if(constexpr_query_size == sizeof_query_unit and 
          constexpr_target_size == sizeof_target_unit and
          constexpr_match_size == sizeof_match_unit){
            query_all_vs_all<constexpr_query_size,constexpr_target_size,
            constexpr_match_size>(query,target,sensitivity,with_simd,short_header,show,mmseqs,ctxsens,correct,correct_ratio,num_threads);
          }
        });
      });
    });
  };
};
#endif