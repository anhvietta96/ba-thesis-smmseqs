#include <iostream>
#include "filter/sorted_q_mer.hpp"

int main()
{
    static constexpr const char nucleotides_upper_lower_ACTG[] = "Aa|Cc|TtUu|Gg";
    static constexpr const size_t seq_len = 3;
    SortedQmer<nucleotides_upper_lower_ACTG,4,seq_len> map;
    
    char* generated_map = map.get_map();
    const size_t size = map.get_size();

    for(size_t i = 0; i < size*seq_len; i++)
    {
        if(i % seq_len == 0)
        {
            std::cout << '\n';
        }
        std::cout << generated_map[i];
    }
    std::cout << std::endl;
}