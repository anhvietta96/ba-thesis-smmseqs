#ifndef INVINT_HASH_HPP
#define INVINT_HASH_HPP

#include<bitset>
#include<vector>
#include "sequences/alphabet.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/unused.hpp"

#ifndef SEED_SIZE
#define SEED_SIZE 64
#endif

constexpr uint64_t get_span(const std::bitset<SEED_SIZE>& seed_bitset)
{
  uint64_t span = 0;
  constexpr_for<0,SEED_SIZE,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      span = idx;
    }
  });
  return span+1;
}

template<const uint64_t span>
constexpr uint64_t get_weight(const std::bitset<SEED_SIZE>& seed_bitset)
{
  uint64_t weight = 0;
  constexpr_for<0,span,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      weight++;
    }
  });
  return weight;
}

constexpr uint64_t count_block(const std::bitset<SEED_SIZE>& seed_bitset, const uint64_t span){
  //static_assert(seed_bitset[span-1]==1);
  uint64_t block_count = 1;
  bool in_block = true;
  for(uint64_t i = span-1; i < seed_bitset.size(); i--){
    block_count += (!in_block && seed_bitset[i]);
    in_block = seed_bitset[i];
  }
  return block_count;
}

struct HashBlock{
  uint64_t to_be_removed = 0;
  uint64_t to_be_added = 0;
  uint64_t rank_removed = 0;
  uint64_t rank_added = 0;
};

template<const uint64_t num_blocks,const uint64_t span, const uint64_t weight>
constexpr std::array<HashBlock,num_blocks> divide_seed(const std::bitset<SEED_SIZE>& seed_bitset){
  std::array<HashBlock,num_blocks> block_schematics{};
  uint64_t arr_idx = 0;
  uint64_t block_start = 0;
  bool in_block = true;
  uint64_t curr_rank = weight;

  for(uint64_t i = span-1; i < seed_bitset.size(); i--){
    if(in_block && !seed_bitset[i]){
      block_schematics[arr_idx].to_be_removed = block_start;
      block_schematics[arr_idx].to_be_added = i;
      block_schematics[arr_idx].to_be_removed = curr_rank+(i-block_start+1);
      block_schematics[arr_idx].to_be_removed = curr_rank;
      arr_idx++;
    }
    if(!in_block && seed_bitset[i]){
      block_start = i;
      in_block = true;
    }
    if(seed_bitset[i]) curr_rank--;
  }
  return block_schematics;
}

template<const uint64_t seed>
class SpacedSeedInterpreter {
  static constexpr const std::bitset<SEED_SIZE> seed_bitset{seed};
  static constexpr const uint64_t span = get_span(seed_bitset);
  static constexpr const uint64_t weight = get_weight<span>(seed_bitset);
  static constexpr const uint64_t padding = SEED_SIZE - span;
  static constexpr const auto num_blocks = count_block(seed_bitset,span);

  HashBlock block_schematics[num_blocks]{};

  public:
  constexpr SpacedSeedInterpreter(){
    uint64_t arr_idx = 0;
    uint64_t block_start = 0;
    bool in_block = true;
    uint64_t curr_rank = weight;

    for(uint64_t i = span-1; i < seed_bitset.size(); i--){
      if(in_block && !seed_bitset[i]){
        block_schematics[arr_idx].to_be_removed = block_start;
        block_schematics[arr_idx].to_be_added = span-1-i;
        block_schematics[arr_idx].rank_removed = curr_rank+(span-1-i-block_start-1);
        block_schematics[arr_idx].rank_added = curr_rank;
        arr_idx++;
        in_block = false;
      }
      if(!in_block && seed_bitset[i]){
        block_start = span-1-i;
        in_block = true;
      }
      if(seed_bitset[i]) curr_rank--;
    }
    block_schematics[arr_idx].to_be_removed = block_start;
    block_schematics[arr_idx].to_be_added = span;
    block_schematics[arr_idx].rank_removed = curr_rank+(span-block_start-1);
    block_schematics[arr_idx].rank_added = curr_rank;
  };

  constexpr uint64_t span_get() const {
    return span;
  }

  constexpr uint64_t weight_get() const {
    return weight;
  }

  constexpr const std::bitset<SEED_SIZE>& seed_bitset_get() const {
    return seed_bitset;
  }

  constexpr uint64_t num_blocks_get() const {
    return num_blocks;
  }

  constexpr const HashBlock* block_schematics_get() const {
    return block_schematics;
  }
};

template<class ScoreClass,const uint64_t seed>
class InvIntHashFunc {
  static constexpr const ScoreClass sc{};
  static constexpr const auto char_spec = sc.character_spec;
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};

  static constexpr const SpacedSeedInterpreter<seed> seed_interpreter{};
  static constexpr const auto block_schematics = seed_interpreter.block_schematics_get();
  static constexpr const auto num_blocks = seed_interpreter.num_blocks_get();
  static constexpr const auto span = seed_interpreter.span_get();
  static constexpr const auto weight = seed_interpreter.weight_get();
  static constexpr const auto seed_bitset = seed_interpreter.seed_bitset_get();

  std::array<uint64_t,weight> factor_table{};
  std::array<uint64_t,weight*undefined_rank> precomputed_table{};

  public:
  constexpr InvIntHashFunc(){
    uint64_t factor = 1;
    for(uint64_t i = 0; i < weight; i++){
      factor_table[i] = factor;
      factor *= undefined_rank;
    }

    for(uint64_t ch = 0; ch < undefined_rank; ch++){
      for(uint64_t rank = 0; rank < weight; rank++){
        precomputed_table[ch*weight+rank] = ch * factor_table[rank];
      }
    }
  }

  uint64_t first_hash_get(const char* seq, const uint64_t seq_len) const {
    uint64_t hashval = 0;
    if(seq_len < span) return 0;
    for(uint64_t i = span-1; i < span; i--){
      if(!seed_bitset[i]) continue;

      hashval *= undefined_rank;
      auto ch = seq[span-1-i];
      if(ch >= undefined_rank){
        std::cout << "Invalid char" << std::endl;
        ch = undefined_rank-1;
      }
      hashval += ch;
    }
    return hashval;
  }

  uint64_t next_hash_get(const char* seq, GTTL_UNUSED const uint64_t seq_len, uint64_t hashval) const {
    for(uint64_t block_num = 0; block_num < num_blocks; block_num++){
      auto ch_to_remove = seq[block_schematics[block_num].to_be_removed];
      auto rank_removed = block_schematics[block_num].rank_removed;
      if(ch_to_remove >= undefined_rank){
        ch_to_remove = undefined_rank-1;
        std::cout << "Invalid char" << std::endl;
      }
      assert(ch_to_remove < undefined_rank && rank_removed < weight);
      hashval -= precomputed_table[ch_to_remove*weight+rank_removed];
    }
    hashval *= undefined_rank;
    for(uint64_t block_num = 0; block_num < num_blocks; block_num++){
      auto ch_to_add = seq[block_schematics[block_num].to_be_added];
      if(ch_to_add >= undefined_rank){
        ch_to_add = undefined_rank-1;
        std::cout << "Invalid char" << std::endl;
      }
      auto rank_added = block_schematics[block_num].rank_added;
      assert(ch_to_add < undefined_rank && rank_added < weight);
      hashval += precomputed_table[ch_to_add*weight+rank_added];
    }
    return hashval;
  }

  constexpr uint64_t span_get() const {
    return span;
  }

  constexpr uint64_t weight_get() const {
    return weight;
  }

  constexpr uint64_t max_hashval_get() const {
    return constexpr_pow(undefined_rank,weight);
  }
};

/*
LiteratureMultiseq

*/

template<class ScoreClass,template<class,uint64_t> class HashFunc,const uint64_t seed>
class Multiseq_Hash {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const HashFunc<ScoreClass,seed> hash_func{};
  static constexpr const auto span = hash_func.span_get();
  static constexpr const auto weight = hash_func.weight_get();

  const uint64_t hashbits = gttl_required_bits(hash_func.max_hashval_get());
  
  static constexpr const int8_t num_bytes_hash = 8;
  std::vector<BytesUnit<num_bytes_hash,3>> bytes_unit_vec;
  GttlBitPacker<num_bytes_hash,3> *hashed_qgram_packer = nullptr;
  std::array<uint64_t,3> num_bits{};

  public:
  Multiseq_Hash(){
    num_bits[0] = hashbits;
  };

  uint64_t hashbits_get() const {
    return hashbits;
  }

  template<const uint64_t sizeof_unit>
  void dbhash(const GttlMultiseq *multiseq,std::vector<BytesUnit<sizeof_unit,3>>&
            target_hash_container,const GttlBitPacker<sizeof_unit,3>& target_packer) const {
    if(multiseq == nullptr) return;
    uint64_t hashval;
    const uint64_t total_seq_num = multiseq->sequences_number_get();
    for(uint64_t seq_num = 0; seq_num < total_seq_num; seq_num++){
      //std::cout << (int) seq_num << std::endl;
      const uint64_t seq_len = multiseq->sequence_length_get(seq_num);
      if(seq_len < span) continue;
      const char* seq = multiseq->sequence_ptr_get(seq_num);
      hashval = hash_func.first_hash_get(seq,seq_len);
      //std::cout << (int)0 << '\t' << (int) hashval << std::endl;
      const BytesUnit<sizeof_unit,3> bu_0{target_packer,
                                                {static_cast<uint64_t>(hashval),
                                                static_cast<uint64_t>(seq_num),
                                                static_cast<uint64_t>(0)}};
      target_hash_container.push_back(bu_0);
      /*target_hash_container.emplace_back(BytesUnit<sizeof_unit,3>{target_packer,
                                              {static_cast<uint64_t>(hashval),
                                              static_cast<uint64_t>(seq_num),
                                              static_cast<uint64_t>(0)}});*/
      for(uint64_t seq_pos = 1; seq_pos <= seq_len - span; seq_pos++){
        hashval = hash_func.next_hash_get(seq+seq_pos-1,seq_len,hashval);
        //std::cout << (int) seq_pos << '\t' << (int) hashval << std::endl;
        const BytesUnit<sizeof_unit,3> bu{target_packer,
                                                {static_cast<uint64_t>(hashval),
                                                static_cast<uint64_t>(seq_num),
                                                static_cast<uint64_t>(seq_pos)}};
        target_hash_container.push_back(bu);
        /*target_hash_container.emplace_back(BytesUnit<sizeof_unit,3>{target_packer,
                                              {static_cast<uint64_t>(hashval),
                                              static_cast<uint64_t>(seq_num),
                                              static_cast<uint64_t>(seq_pos)}});*/
      }
    }
  };

  template<const uint64_t sizeof_unit>
  void dbsort(std::vector<BytesUnit<sizeof_unit,3>>& target_hash_container,const uint64_t hash_bits) const {
    if constexpr (sizeof_unit == 8){
      ska_lsb_radix_sort<uint64_t>(hash_bits,
                              reinterpret_cast<uint64_t *>
                              (target_hash_container.data()),
                              target_hash_container.size());
    }
    else {
      ska_large_lsb_small_radix_sort(sizeof_unit,hash_bits,reinterpret_cast<uint8_t *>(
                                  target_hash_container.data()),target_hash_container.size(),
                                  false);
    }
  }

  void hash(const GttlMultiseq *multiseq){
    assert(multiseq != nullptr);
    //LiterateMultiseq<char_spec,undefined_rank> literate_multiseq{*multiseq};
    //literate_multiseq.perform_sequence_encoding();
    
    const uint64_t total_seq_num = multiseq->sequences_number_get();
    const uint64_t seq_len_bits = multiseq->sequences_length_bits_get();
    const uint64_t seq_num_bits = multiseq->sequences_number_bits_get();
    //std::cout << (int) seq_len_bits << '\t' << (int) seq_num_bits << std::endl;
    //const uint64_t hash_bits = gttl_required_bits<uint64_t>(constexpr_pow(undefined_rank,weight));
    const uint64_t hash_bits = 8*num_bytes_hash-seq_len_bits-seq_num_bits;
    num_bits[0] = hash_bits;
    num_bits[1] = seq_num_bits;
    num_bits[2] = seq_len_bits;

    hashed_qgram_packer = new GttlBitPacker<num_bytes_hash,3>({static_cast<int>(hash_bits),
                                                              static_cast<int>(seq_num_bits),
                                                              static_cast<int>(seq_len_bits)});
    uint64_t hashval;
    
    for(uint64_t seq_num = 0; seq_num < total_seq_num; seq_num++){
      const uint64_t seq_len = multiseq->sequence_length_get(seq_num);
      const char* seq = multiseq->sequence_ptr_get(seq_num);
      hashval = hash_func.first_hash_get(seq,seq_len);
      /*const BytesUnit<num_bytes_hash,3> bu_0{*hashed_qgram_packer,
                                            {static_cast<uint64_t>(hashval),
                                            static_cast<uint64_t>(seq_num),
                                            static_cast<uint64_t>(0)}};
      bytes_unit_vec.push_back(bu_0);*/
      bytes_unit_vec.emplace_back(BytesUnit<num_bytes_hash,3>{
                                     *hashed_qgram_packer,
                                    {static_cast<uint64_t>(hashval),
                                    static_cast<uint64_t>(seq_num),
                                    static_cast<uint64_t>(0)}});
      for(uint64_t seq_pos = 1; seq_pos <= seq_len - span; seq_pos++){
        hashval = hash_func.next_hash_get(seq+seq_pos-1,seq_len,hashval);
        /*const BytesUnit<num_bytes_hash,3> bu{*hashed_qgram_packer,
                                            {static_cast<uint64_t>(hashval),
                                            static_cast<uint64_t>(seq_num),
                                            static_cast<uint64_t>(seq_pos)}};
        bytes_unit_vec.push_back(bu);*/
        bytes_unit_vec.emplace_back(BytesUnit<num_bytes_hash,3>{
                                    *hashed_qgram_packer,
                                    {static_cast<uint64_t>(hashval),
                                    static_cast<uint64_t>(seq_num),
                                    static_cast<uint64_t>(seq_pos)}});
      }
    }
    ska_lsb_radix_sort<uint64_t>(hash_bits,
                              reinterpret_cast<uint64_t *>
                              (bytes_unit_vec.data()),
                              bytes_unit_vec.size());
  };

  ~Multiseq_Hash(){
    if(hashed_qgram_packer) delete hashed_qgram_packer;
  }

  uint64_t size()const{
    return bytes_unit_vec.size();
  }

  const std::array<uint64_t,3>& num_bits_get() const {
    return num_bits;
  }

  const std::vector<BytesUnit<num_bytes_hash,3>>& bytes_unit_vec_get() const {
    return bytes_unit_vec;
  }

  const BytesUnit<num_bytes_hash,3>& bytes_unit_get(const uint64_t idx) const {
    return bytes_unit_vec[idx];
  }

  const GttlBitPacker<num_bytes_hash,3> *packer_get() const {
    return hashed_qgram_packer;
  }
};

#endif