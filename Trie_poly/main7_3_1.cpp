#include "trie_poly.hpp"
#include <iostream>

int main() 
{
    using namespace order;
    // exersise 2.6.5 from Кокс book
    std::cout << "exersise 2.6.5 from Кокс book" << std::endl;
    //  a)
    std::cout << "a)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({2, 0, 1}, 4);
        f.add_term({0, 2, 0}, -7);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({1, 1, 2}, 1);  // xy
        g.add_term({1, 0, 4}, 3);  // 1

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto S = f.S_poly(g, Lex{});

        std::cout << "S(f, g) = "; S.print();
        std::cout << std::endl;

    }

    //  b)
    std::cout << "b)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({4, 1, 0}, 1);
        f.add_term({0, 0, 2}, -1);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({1, 0, 2}, 3);  // xy
        g.add_term({0, 1, 0}, -1);  // 1

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto S = f.S_poly(g, Lex{});

        std::cout << "S(f, g) = "; S.print();
        std::cout << std::endl;
    }

    //  c)
    std::cout << "c)" << std::endl;
    {
        PolyTrie<std::complex<double>> f({"x", "y", "z"});
        f.add_term({7, 2, 1}, std::complex<double>(1, 0));
        f.add_term({1, 1, 1}, std::complex<double>(0, 2));

        PolyTrie<std::complex<double>> g({"x", "y", "z"});
        g.add_term({7, 2, 1}, std::complex<double>(2, 0));  // xy
        g.add_term({0, 1, 0}, -1);  // 1

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto S = f.S_poly(g, Lex{});

        std::cout << "S(f, g) = "; S.print();
        std::cout << std::endl;
    }

    //  d)
    std::cout << "d)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({1, 1, 0}, 1);
        f.add_term({0, 0, 3}, 1);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({0, 0, 2}, 1);  // xy
        g.add_term({0, 0, 1}, -3);  // 1

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto S = f.S_poly(g, Lex{});

        std::cout << "S(f, g) = "; S.print();
        std::cout << std::endl;
    }

    return 0;
}