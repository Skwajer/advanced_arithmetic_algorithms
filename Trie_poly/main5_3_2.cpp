#include "trie_poly.hpp"
#include <utility>

int main()
{
    // (a) test: f = x^7 y^2 + x^3 y^2 - y + 1
    // F = {x y^2 - x, x - y^3}
    std::cout << "========== DIVISION TASK 1 ==========" << std::endl;
    {
        std::cout << "====== division test (a) (lex) f---g1---g2---r ======" << std::endl;

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

        std::cout << "f = ";
        f.print();

        std::cout << "g1 = ";
        g1.print();
        std::cout << "g2 = ";
        g2.print();
        std::cout << std::endl;

        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g1));
        F.emplace_back(std::move(g2));

        auto [Q, r] = f.divide(F, order::Lex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    {
        std::cout << "====== division test (a) (grlex) f---g1---g2---r ======" << std::endl;

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

        std::cout << "f = ";
        f.print();

        std::cout << "g1 = ";
        g1.print();
        std::cout << "g2 = ";
        g2.print();
        std::cout << std::endl;

        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g1));
        F.emplace_back(std::move(g2));

        auto [Q, r] = f.divide(F, order::GrLex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    {
        std::cout << "====== division test (b) (lex) f---g2---g1---r ======" << std::endl;

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

        std::cout << "f = ";
        f.print();

        std::cout << "g2 = ";
        g2.print();
        std::cout << "g1 = ";
        g1.print();
        std::cout << std::endl;

        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g2));
        F.emplace_back(std::move(g1));

        auto [Q, r] = f.divide(F, order::Lex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    {
        std::cout << "====== division test (b) (grlex) f---g2---g1---r ======" << std::endl;

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

        std::cout << "f = ";
        f.print();

        std::cout << "g2 = ";
        g2.print();
        std::cout << "g1 = ";
        g1.print();
        std::cout << std::endl;

        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g2));
        F.emplace_back(std::move(g1));

        auto [Q, r] = f.divide(F, order::GrLex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }
    std::cout << "========== DIVISION TASK 2 ==========" << std::endl;
    {
        std::cout << "====== division test (a) (lex) f---g1---g2---g3---r ======" << std::endl;

        PolyTrie<double> f({"x","y", "z"});
        f.add_term({1,2, 2}, 1.0);
        f.add_term({1,1, 0}, 1.0);
        f.add_term({0,1, 1}, -1.0);

        PolyTrie<double> g1({"x","y", "z"});
        g1.add_term({1,0, 0}, 1.0);
        g1.add_term({0,2, 0}, -1.0);

        PolyTrie<double> g2({"x","y", "z"});
        g2.add_term({0,1,0}, 1.0);
        g2.add_term({0,0,3}, -1.0);

        PolyTrie<double> g3({"x","y", "z"});
        g3.add_term({0,0,2}, 1.0);
        g3.add_term({0,0,0}, -1.0);

        std::cout << "f = ";
        f.print();

        std::cout << "g1 = ";
        g1.print();
        std::cout << "g2 = ";
        g2.print();
        std::cout << "g3 = ";
        g3.print();
        std::cout << std::endl;



        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g1));
        F.emplace_back(std::move(g2));
        F.emplace_back(std::move(g3));

        auto [Q, r] = f.divide(F, order::Lex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    {
        std::cout << "====== division test (b) (lex) f---g2---g3---g1---r ======" << std::endl;

        PolyTrie<double> f({"x","y", "z"});
        f.add_term({1,2, 2}, 1.0);
        f.add_term({1,1, 0}, 1.0);
        f.add_term({0,1, 1}, -1.0);

        PolyTrie<double> g1({"x","y", "z"});
        g1.add_term({1,0, 0}, 1.0);
        g1.add_term({0,2, 0}, -1.0);

        PolyTrie<double> g2({"x","y", "z"});
        g2.add_term({0,1,0}, 1.0);
        g2.add_term({0,0,3}, -1.0);

        PolyTrie<double> g3({"x","y", "z"});
        g3.add_term({0,0,2}, 1.0);
        g3.add_term({0,0,0}, -1.0);

        std::cout << "f = ";
        f.print();

        std::cout << "g1 = ";
        g1.print();
        std::cout << "g2 = ";
        g2.print();
        std::cout << "g3 = ";
        g3.print();
        std::cout << std::endl;



        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g2));
        F.emplace_back(std::move(g3));
        F.emplace_back(std::move(g1));

        auto [Q, r] = f.divide(F, order::Lex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    {
        std::cout << "====== division test (b) (lex) f---g3---g1---g2---r ======" << std::endl;

        PolyTrie<double> f({"x","y", "z"});
        f.add_term({1,2, 2}, 1.0);
        f.add_term({1,1, 0}, 1.0);
        f.add_term({0,1, 1}, -1.0);

        PolyTrie<double> g1({"x","y", "z"});
        g1.add_term({1,0, 0}, 1.0);
        g1.add_term({0,2, 0}, -1.0);

        PolyTrie<double> g2({"x","y", "z"});
        g2.add_term({0,1,0}, 1.0);
        g2.add_term({0,0,3}, -1.0);

        PolyTrie<double> g3({"x","y", "z"});
        g3.add_term({0,0,2}, 1.0);
        g3.add_term({0,0,0}, -1.0);

        std::cout << "f = ";
        f.print();

        std::cout << "g1 = ";
        g1.print();
        std::cout << "g2 = ";
        g2.print();
        std::cout << "g3 = ";
        g3.print();
        std::cout << std::endl;



        std::vector<PolyTrie<double>> F;
        F.emplace_back(std::move(g3));
        F.emplace_back(std::move(g1));
        F.emplace_back(std::move(g2));

        auto [Q, r] = f.divide(F, order::Lex{});


        for (size_t i = 0; i < Q.size(); ++i)
        {
            std::cout << "q" << i+1 << " = ";
            Q[i].print();
        }

        std::cout << "remainder r = ";
        r.print();
    }

    return 0;
}