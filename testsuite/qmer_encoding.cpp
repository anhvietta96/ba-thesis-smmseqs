#include "filter/sorted_q_mer.hpp"

int main()
{
  static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
  constexpr const SortedQmer<nucleotides_upper_lower_ACTG,4,3> map_nuc_3{};

  /*static constexpr const char amino_acids[]
      = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
    constexpr const SortedQmer<amino_acids,20,2> map_prot_2{};*/

  constexpr const char test[] = "CCG";
  constexpr const size_t sz = sizeof(test)-1;
  char sorted[sz+1] = { 0 };
  size_t permutation[sz];
  bool sort = map_nuc_3.sort(test,sorted,permutation);
  std::cout << sort << std::endl;
  std::cout << sorted << std::endl;


  size_t code = map_nuc_3.encode(sorted);
  std::cout << code << std::endl;
}