#include "sequences/alphabet.hpp"
#include <cassert>

size_t binom(size_t a,size_t b)
{
  if(2*b > a && a != b)
  {
      return binom(a,a-b);
  }
  if(b == 0 || a == b)
  {
    return 1;
  }
  return binom(a-1,b-1) + binom(a-1,b);
}

constexpr const size_t calc_max_num_seq(size_t seq_len, size_t alphabet_size)
{
  size_t count = 0;
  for(size_t idx = 1; idx <= alphabet_size; idx++)
  {
    if(idx <= seq_len)
    {
      count += binom(seq_len-1,idx-1)*binom(alphabet_size,idx);
    }
  }
  return count;
}

template <const char* char_spec,uint8_t undefined_rank,const size_t seq_len>
constexpr const char*** create_map()
{
  const GttlAlphabet<char_spec,undefined_rank> alpha{};
  constexpr const size_t total_seq_num = calc_max_num_seq(seq_len,alpha.size());
  char map[total_seq_num][seq_len];
  char temp_seq[seq_len] = { 0 };
  size_t code = 0;
  while(temp_seq[0] != alpha.size() - 1)
  {
    std::copy(temp_seq,temp_seq+seq_len,map[code]);
    if(temp_seq[seq_len-1] < alpha.size() - 1)
    {
      temp_seq[seq_len-1]++;
    }
    else
    {
      size_t idx_to_fix = seq_len-1;
      char value_to_fix = 0;
      for(size_t idx = seq_len-2; idx > 0; idx--)
      {
        if(temp_seq[idx] < alpha.size()-1)
        {
          idx_to_fix = idx;
          value_to_fix = temp_seq[idx]+1;
          break;
        }
      }
      for(size_t idx = idx_to_fix; idx < seq_len; idx++)
      {
        temp_seq[idx] = value_to_fix;
      }
    }
  }
  return &map;
}

template<const char* char_spec,uint8_t undefined_rank,const size_t seq_len>
class SortedQmer
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const char*** map = create_map<char_spec,undefined_rank,seq_len>();
  static constexpr const size_t size = calc_max_num_seq(seq_len,alpha.size());

  public:
  size_t get_size()
  {
    return this->size;
  }

  const char*** get_map()
  {
    return this->map;
  }
  /*char* sort(const char* seq)
  {
    
  }

  size_t encode(const char* seq)
  {

  }

  size_t decode(size_t code)
  {

  }*/
};