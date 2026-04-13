#include "trie_poly.hpp"
#include <iostream>

template<typename Comparator>
bool is_strictly_less(const std::vector<int>& a, const std::vector<int>& b, Comparator comp) 
{
    return comp(a, b) && !(comp(b, a));
}

int main() 
{
    using namespace order;
    // exersise 2.6.5 from Кокс book
    std::cout << "exercise 2.6.5 from Кокс book" << std::endl;
    //  a)
    std::cout << "a)" << std::endl;
    std::cout << "a)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({2, 0, 1}, 4);
        f.add_term({0, 2, 0}, -7);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({1, 1, 2}, 1);
        g.add_term({1, 0, 4}, 3);

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto lf = f.leading_monomial(Lex{});
        auto lg = g.leading_monomial(Lex{});
        
        std::vector<int> lcm(lf.size());
        for (size_t i = 0; i < lf.size(); ++i)
            lcm[i] = std::max(lf[i], lg[i]);

        auto S = f.S_poly(g, Lex{});
        auto multideg_S = S.leading_monomial(Lex{});

        std::cout << "S(f, g) = "; S.print();
        
        bool ok = is_strictly_less(multideg_S, lcm, Lex{});
        std::cout << " => " << (ok ? "multideg(S) < LCM" : "multideg(S) >= LCM") << std::endl;
        std::cout << std::endl;
    }

    //  b)
    std::cout << "b)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({4, 1, 0}, 1);
        f.add_term({0, 0, 2}, -1);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({1, 0, 2}, 3);
        g.add_term({0, 1, 0}, -1);

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto lf = f.leading_monomial(Lex{});
        auto lg = g.leading_monomial(Lex{});
        
        std::vector<int> lcm(lf.size());
        for (size_t i = 0; i < lf.size(); ++i)
            lcm[i] = std::max(lf[i], lg[i]);

        auto S = f.S_poly(g, Lex{});
        auto multideg_S = S.leading_monomial(Lex{});

        std::cout << "S(f, g) = "; S.print();
        
        bool ok = is_strictly_less(multideg_S, lcm, Lex{});
        std::cout << " => " << (ok ? "multideg(S) < LCM" : "multideg(S) >= LCM") << std::endl;
        std::cout << std::endl;
    }

    //  c)
    std::cout << "c)" << std::endl;
    {
        PolyTrie<std::complex<double>> f({"x", "y", "z"});
        f.add_term({7, 2, 1}, std::complex<double>(1, 0));
        f.add_term({1, 1, 1}, std::complex<double>(0, 2));

        PolyTrie<std::complex<double>> g({"x", "y", "z"});
        g.add_term({7, 2, 1}, std::complex<double>(2, 0));
        g.add_term({0, 1, 0}, -1);

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto lf = f.leading_monomial(Lex{});
        auto lg = g.leading_monomial(Lex{});
        
        std::vector<int> lcm(lf.size());
        for (size_t i = 0; i < lf.size(); ++i)
            lcm[i] = std::max(lf[i], lg[i]);

        auto S = f.S_poly(g, Lex{});
        auto multideg_S = S.leading_monomial(Lex{});

        std::cout << "S(f, g) = "; S.print();
        
        bool ok = is_strictly_less(multideg_S, lcm, Lex{});
        std::cout << " => " << (ok ? "multideg(S) < LCM" : "multideg(S) >= LCM") << std::endl;
        std::cout << std::endl;
    }

    //  d)
    std::cout << "d)" << std::endl;
    {
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({1, 1, 0}, 1);
        f.add_term({0, 0, 3}, 1);

        PolyTrie<double> g({"x", "y", "z"});
        g.add_term({0, 0, 2}, 1);
        g.add_term({0, 0, 1}, -3);

        std::cout << "f = "; f.print();
        std::cout << "g = "; g.print();

        auto lf = f.leading_monomial(Lex{});
        auto lg = g.leading_monomial(Lex{});
        
        std::vector<int> lcm(lf.size());
        for (size_t i = 0; i < lf.size(); ++i)
            lcm[i] = std::max(lf[i], lg[i]);

        auto S = f.S_poly(g, Lex{});
        auto multideg_S = S.leading_monomial(Lex{});

        std::cout << "S(f, g) = "; S.print();
        
        bool ok = is_strictly_less(multideg_S, lcm, Lex{});
        std::cout << " => " << (ok ? "multideg(S) < LCM" : "multideg(S) >= LCM") << std::endl;
        std::cout << std::endl;
    }

    return 0;
}