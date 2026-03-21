#include "trie_poly.hpp"
#include <complex>

using complex = std::complex<double>;

int main()
{
    PolyTrie<double> f1_R2({"alpha", "betta", "gamma"});
    f1_R2.add_term({2, 3, 7}, -3.5);
    f1_R2.add_term({1, 4, -4}, 5);
    f1_R2.add_term({0, 0, 0}, -11.3);

    PolyTrie<double> f2_R2({"alpha", "betta", "gamma"});
    f2_R2.add_term({2, 3, 7}, -3.5);
    f2_R2.add_term({1, 4, -4}, 5);
    f2_R2.add_term({0, 0, 0}, -11.3);
    f2_R2.add_term({3, 3, 3}, -25.0);

    f1_R2 += f2_R2;

    f1_R2.print();
    return 0;
}