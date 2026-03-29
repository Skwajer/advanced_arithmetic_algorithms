#include "trie_poly.hpp"
#include <utility>

int main()
{
    // (a) test: f = x^7 y^2 + x^3 y^2 - y + 1
    // F = {x y^2 - x, x - y^3}

    std::cout << "====== division test (a) ======" << std::endl;

    PolyTrie<double> f({"x","y"});
    f.add_term({7,2}, 1.0);
    f.add_term({3,2}, 1.0);
    f.add_term({0,1}, -1.0);
    f.add_term({0,0}, 1.0);

    PolyTrie<double> g1({"x","y"});
    g1.add_term({1,2}, 1.0);
    g1.add_term({1,0}, -1.0);

    PolyTrie<double> g2({"x","y"});
    g2.add_term({1,0}, 1.0);
    g2.add_term({0,3}, -1.0);

    std::vector<PolyTrie<double>> F;
    F.emplace_back(std::move(g1));
    F.emplace_back(std::move(g2));

    auto [Q, r] = f.divide(F, order::Lex{});

    std::cout << "f = ";
    f.print();

    for (size_t i = 0; i < Q.size(); ++i)
    {
        std::cout << "q" << i+1 << " = ";
        Q[i].print();
    }

    std::cout << "remainder r = ";
    r.print();
    return 0;
}