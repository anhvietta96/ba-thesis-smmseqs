#include "sequences/alphabet.hpp"
#include "utilities/constexpr_for.hpp"
#include <cassert>
#include <iostream>

constexpr size_t get_span(const char *seed, size_t index)
{
    return !(seed[index] - '\0') ? index : get_span(seed,index+1);
}

constexpr size_t get_weight(const char *seed, size_t index, size_t count)
{
    return !(seed[index] - '\0') ? count : !(seed[index] - '0') ? get_weight(seed,index+1,count) : get_weight(seed,index+1,count+1);
}

template<const char* char_spec,uint8_t undefined_rank,const char *seed>
class GttlSpacedSeed
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const size_t span = get_span(seed,0);
  static constexpr const size_t weight = get_weight(seed,0,0);

  public:
  constexpr GttlSpacedSeed(){};

  size_t encode(const char* seq_ptr, const size_t seq_len) const
  {
    //assert(seq_len == span);
    if(seq_len != span)
    {
        std::cout << seq_len << "      " << span << std::endl;
    }
    
    size_t code = 0;
    static constexpr const size_t alphabet_size = alpha.size();
    constexpr_for<0,span,1>([&] (auto idx)
    {
      if constexpr (seed[idx]-'0')
      {
        code *= alphabet_size;
        code += alpha.char_to_rank(seq_ptr[idx]);
      }
    });
    return code;
  }
};