#include "gf2n.hpp"

int main()
{
    GF2n field(4, 0b1001);
    GF2n::u64 test_to_str = 0b010010101001110;
    std::cout << field.to_string(test_to_str) << std::endl;

    return 0;
}