#include "trie_poly.hpp"
#include <complex>

using complex = std::complex<double>;

int main()
{
    PolyTrie<double> f_R2({"alpha", "betta", "gamma"});
    f_R2.add_term({2, 3, 7}, -3.5);
    f_R2.add_term({1, 4, -4}, 5);
    f_R2.add_term({0, 0, 0}, -11.3);
    f_R2.print();
    return 0;
}