#ifndef INVINT_HASH_HPP
#define INVINT_HASH_HPP

#include <bitset>
#include "sequences/alphabet.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/mathsupport.hpp"

#ifndef SEED_SIZE
#define SEED_SIZE 64
#endif

constexpr size_t get_span(const std::bitset<SEED_SIZE> seed_bitset)
{
  size_t span = 0;
  constexpr_for<0,SEED_SIZE,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      span = idx;
    }
  });
  return span+1;
}

template<const size_t span>
constexpr size_t get_weight(const std::bitset<SEED_SIZE>& seed_bitset)
{
  size_t weight = 0;
  constexpr_for<0,span,1>([&] (auto idx)
  {
    if(seed_bitset[idx] == 1)
    {
      weight++;
    }
  });
  return weight;
}

constexpr size_t count_block(const std::bitset<SEED_SIZE>& seed_bitset, const size_t span){
  //static_assert(seed_bitset[span-1]==1);
  size_t block_count = 1;
  bool in_block = true;
  for(size_t i = span-1; i < seed_bitset.size(); i--){
    block_count += (!in_block && seed_bitset[i]);
    in_block = seed_bitset[i];
  }
  return block_count;
}

struct HashBlock{
  size_t to_be_removed = 0;
  size_t to_be_added = 0;
  size_t rank_removed = 0;
  size_t rank_added = 0;
};

template<const size_t num_blocks,const size_t span, const size_t weight>
constexpr std::array<HashBlock,num_blocks> divide_seed(const std::bitset<SEED_SIZE>& seed_bitset){
  std::array<HashBlock,num_blocks> block_schematics{};
  size_t arr_idx = 0;
  size_t block_start = 0;
  bool in_block = true;
  size_t curr_rank = weight;

  for(size_t i = span-1; i < seed_bitset.size(); i--){
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

template<const size_t seed>
class SpacedSeedInterpreter {
  static constexpr const std::bitset<SEED_SIZE> seed_bitset{seed};
  static constexpr const size_t span = get_span(seed_bitset);
  static constexpr const size_t weight = get_weight<span>(seed_bitset);
  static constexpr const size_t padding = SEED_SIZE - span;
  static constexpr const auto num_blocks = count_block(seed_bitset,span);

  HashBlock block_schematics[num_blocks]{};

  public:
  constexpr SpacedSeedInterpreter(){
    size_t arr_idx = 0;
    size_t block_start = 0;
    bool in_block = true;
    size_t curr_rank = weight;

    for(size_t i = span-1; i < seed_bitset.size(); i--){
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

  constexpr size_t span_get() const {
    return span;
  }

  constexpr size_t weight_get() const {
    return weight;
  }

  constexpr const std::bitset<SEED_SIZE>& seed_bitset_get() const {
    return seed_bitset;
  }

  constexpr size_t num_blocks_get() const {
    return num_blocks;
  }

  constexpr const HashBlock* block_schematics_get() const {
    return block_schematics;
  }
};

template<class ScoreClass,const size_t seed>
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

  std::array<size_t,weight> factor_table{};
  std::array<size_t,weight*undefined_rank> precomputed_table{};

  public:
  constexpr InvIntHashFunc(){
    size_t factor = 1;
    for(size_t i = 0; i < weight; i++){
      factor_table[i] = factor;
      factor *= undefined_rank;
    }

    for(size_t ch = 0; ch < undefined_rank; ch++){
      for(size_t rank = 0; rank < weight; rank++){
        precomputed_table[ch*weight+rank] = ch * factor_table[rank];
      }
    }
  }

  size_t first_hash_get(const char* seq, const size_t seq_len) const {
    size_t hashval = 0;
    if(seq_len < span) return 0;
    for(size_t i = span-1; i < span; i--){
      if(!seed_bitset[i]) continue;

      hashval *= undefined_rank;
      hashval += alpha.char_to_rank(seq[span-1-i]);
    }
    return hashval;
  }

  size_t next_hash_get(const char* seq, const size_t seq_len, size_t hashval) const {
    /*for(size_t i = 0; i < num_blocks; i++){
      const auto block_0 = block_schematics[i];
      std::cout << (int) block_0.to_be_added << '\t' << (int) block_0.to_be_removed << '\t' 
      << (int) block_0.rank_added << '\t' << (int) block_0.rank_removed << '\t' << std::endl;
    }*/
    
    for(size_t block_num = 0; block_num < num_blocks; block_num++){
      auto ch_to_remove = alpha.char_to_rank(seq[block_schematics[block_num].to_be_removed]);
      auto rank_removed = block_schematics[block_num].rank_removed;
      /*std::cout << (int) ch_to_remove << '\t' << (int) rank_removed << '\t' 
      << (int) precomputed_table[ch_to_remove*weight+rank_removed] << std::endl;*/
      assert(ch_to_remove < undefined_rank && rank_removed < weight);
      hashval -= precomputed_table[ch_to_remove*weight+rank_removed];
    }
    hashval *= undefined_rank;
    for(size_t block_num = 0; block_num < num_blocks; block_num++){
      auto ch_to_add = alpha.char_to_rank(seq[block_schematics[block_num].to_be_added]);
      auto rank_added = block_schematics[block_num].rank_added;
      /*std::cout << (int) ch_to_add << '\t' << (int) rank_added << '\t' 
      << (int) precomputed_table[ch_to_add*weight+rank_added] << std::endl;*/
      assert(ch_to_add < undefined_rank && rank_added < weight);
      hashval += precomputed_table[ch_to_add*weight+rank_added];
    }
    return hashval;
  }

  constexpr size_t span_get() const {
    return span;
  }

  constexpr size_t weight_get() const {
    return weight;
  }
};

template<class ScoreClass,template<class,size_t> class HashFunc,const size_t seed>
class Multiseq_Hash {
  static constexpr const ScoreClass sc{};
  static constexpr const auto undefined_rank = sc.num_of_chars;
  static constexpr const auto char_spec = sc.character_spec;

  static constexpr const HashFunc<ScoreClass,seed> hash_func{};
  static constexpr const auto span = hash_func.span_get();
  static constexpr const auto weight = hash_func.weight_get();

  std::vector<BytesUnit<9,3>> bytes_unit_vec;
  GttlBitPacker<9,3> *hashed_qgram_packer = nullptr;

  public:
  Multiseq_Hash(){};

  void hash(const GttlMultiseq *multiseq){
    assert(multiseq != nullptr);
    //LiterateMultiseq<char_spec,undefined_rank> literate_multiseq{*multiseq};
    //literate_multiseq.perform_sequence_encoding();
    
    const size_t total_seq_num = multiseq->sequences_number_get();
    const size_t seq_len_bits = multiseq->sequences_length_bits_get();
    const size_t seq_num_bits = multiseq->sequences_number_bits_get();
    //std::cout << (int) seq_len_bits << '\t' << (int) seq_num_bits << std::endl;
    //const size_t hash_bits = gttl_required_bits<size_t>(constexpr_pow(undefined_rank,weight));
    const size_t hash_bits = 71-seq_len_bits-seq_num_bits;
    std::cout << "Query" << std::endl;
    hashed_qgram_packer = new GttlBitPacker<9,3>({static_cast<uint64_t>(seq_num_bits),
                                                  static_cast<uint64_t>(seq_len_bits),
                                                  static_cast<uint64_t>(hash_bits)});
    std::cout << "Query" << std::endl;
    size_t hashval;
    
    for(size_t seq_num = 0; seq_num < total_seq_num; seq_num++){
      const size_t seq_len = multiseq->sequence_length_get(seq_num);
      const auto seq = multiseq->sequence_ptr_get(seq_num);
      hashval = hash_func.first_hash_get(seq,seq_len);
      //std::cout << (int)0 << '\t' << (int) hashval << std::endl;
      const BytesUnit<9,3> bu_0{*hashed_qgram_packer,{static_cast<uint64_t>(seq_num),
                                                      static_cast<uint64_t>(0),
                                                      static_cast<uint64_t>(hashval)}};
      bytes_unit_vec.push_back(bu_0);
      for(size_t seq_pos = 1; seq_pos <= seq_len - span; seq_pos++){
        hashval = hash_func.next_hash_get(seq+seq_pos-1,seq_len,hashval);
        //std::cout << (int) seq_pos << '\t' << (int) hashval << std::endl;
        const BytesUnit<9,3> bu{*hashed_qgram_packer,{static_cast<uint64_t>(seq_num),
                                                      static_cast<uint64_t>(seq_pos),
                                                      static_cast<uint64_t>(hashval)}};
        bytes_unit_vec.push_back(bu);
      }
    }
  };

  ~Multiseq_Hash(){
    delete hashed_qgram_packer;
  }

  size_t size()const{
    return bytes_unit_vec.size();
  }

  const std::vector<BytesUnit<9,3>>& bytes_unit_vec_get() const {
    return bytes_unit_vec;
  }

  const BytesUnit<9,3>& bytes_unit_get(const size_t idx) const {
    return bytes_unit_vec[idx];
  }

  const GttlBitPacker<9,3> *packer_get() const {
    return hashed_qgram_packer;
  }
};

#endif