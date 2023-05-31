#include "filter/sorted_q_mer.hpp"
#include "sequences/alphabet.hpp"

template<const char* char_spec,uint8_t undefined_rank,const size_t qgram_length>
class QgramProc
{
  private:
  static constexpr const GttlAlphabet<char_spec,undefined_rank> alpha{};
  static constexpr const SortedQmer<char_spec,undefined_rank,qgram_length> map;

  public:
  constexpr QgramProc(){};

  bool sort(const char* qgram_ptr,char* sorted_qgram, size_t* permutation) const
  {
    for(size_t idx = 0; idx < qgram_length; idx++)
    {
      sorted_qgram[idx] = qgram_ptr[idx];
      permutation[idx] = idx;
    }
    bool swapped = false;
    for(size_t pm = 0; pm < qgram_length; pm++)
    {
      for(size_t pl = pm; pl > 0 and alpha.char_to_rank(sorted_qgram[pl-1]) 
      > alpha.char_to_rank(sorted_qgram[pl]); pl--)
      {
        const uint8_t tmp_cc = sorted_qgram[pl-1];
        sorted_qgram[pl-1] = sorted_qgram[pl];
        sorted_qgram[pl] = tmp_cc;

        const size_t tmp_t = permutation[pl-1];
        permutation[pl-1] = permutation[pl];
        permutation[pl] = tmp_t;

        swapped = true;
      }
    }
    return swapped;
  }  

  size_t encode(const char* qgram) const
  {
    size_t start=0,end=ms_size-1;
    while(start!=end)
    {
      const size_t mid = (start+end)/2;
      const char comp = compare_qgram(map[mid],qgram,alpha);
      if(comp < 0)
      {
        end = mid;
      }
      else if(comp > 0)
      {
        start = mid;
      }
      else
      {
        return mid;
      }
    }
    return start;
  }
};