#include <iostream>
#include "filter/sorted_q_mer.hpp"

int main(){
    static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
    SortedQmer<nucleotides_upper_lower_ACTG,4,2> map;
    const char** generated_map = *(map.get_map());
    const size_t size = map.get_size();
    for(size_t i = 0; i < size; i++)
    {
        std::cout << generated_map[i] << std::endl;
    }
}