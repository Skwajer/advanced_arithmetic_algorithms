#include "trie_poly.hpp"
#include <iostream>

int main()
{
    auto print_monomial = [](const std::vector<int>& degs,
                             const std::vector<std::string>& vars)
    {
        bool first = true;
        for (size_t i = 0; i < degs.size(); ++i)
        {
            if (degs[i] != 0)
            {
                if (!first) std::cout << " * ";
                std::cout << vars[i];
                if (degs[i] != 1)
                    std::cout << "^" << degs[i];
                first = false;
            }
        }
        if (first) std::cout << "1";
    };

    auto print_term = [&](const std::vector<int>& degs,
                          double coeff,
                          const std::vector<std::string>& vars)
    {
        std::cout << coeff << " * ";
        print_monomial(degs, vars);
    };


    {
        std::cout << "=============== TASK 1 ===============" << std::endl << std::endl;
        std::cout << "x > y > z"<< std::endl;
        // a part
        std::cout << "====== (a) ======" << std::endl;
        PolyTrie<double> f_a({"x", "y", "z"});
        f_a.add_term({1, 0, 0}, 2.0);
        f_a.add_term({0, 1, 0}, 3.0);
        f_a.add_term({0 ,0 ,1}, 1.0);
        f_a.add_term({2, 0, 0}, 1.0);
        f_a.add_term({0, 0, 2}, -1.0);
        f_a.add_term({3, 0, 0}, 1.0);

        std::vector<std::string> vars_a = {"x","y","z"};

        auto multideg_f_a_lex = f_a.find_multideg(order::Lex{});
        auto multideg_f_a_grlex = f_a.find_multideg(order::GrLex{});
        auto multideg_f_a_grevlex = f_a.find_multideg(order::GrevLex{});

        auto lm_f_a_lex = f_a.leading_monomial(order::Lex{});
        auto lc_f_a_lex = f_a.leading_coeff(order::Lex{});
        auto lt_f_a_lex = f_a.leading_term(order::Lex{});

        auto lm_f_a_grlex = f_a.leading_monomial(order::GrLex{});
        auto lc_f_a_grlex = f_a.leading_coeff(order::GrLex{});
        auto lt_f_a_grlex = f_a.leading_term(order::GrLex{});

        auto lm_f_a_grevlex = f_a.leading_monomial(order::GrevLex{});
        auto lc_f_a_grevlex = f_a.leading_coeff(order::GrevLex{});
        auto lt_f_a_grevlex = f_a.leading_term(order::GrevLex{});

        std::cout << "multideg(f_a) (lex) = ";
        for (auto d : multideg_f_a_lex) std::cout << d << " ";
        std::cout << "\nmultideg(f_a) (grlex) = ";
        for (auto d : multideg_f_a_grlex) std::cout << d << " ";
        std::cout << "\nmultideg(f_a) (grevlex) = ";
        for (auto d : multideg_f_a_grevlex) std::cout << d << " ";
        std::cout << "\n";

        std::cout << "lm(f_a) (lex) = ";
        print_monomial(lm_f_a_lex, vars_a);
        std::cout << "\nlc(f_a) (lex) = " << lc_f_a_lex;
        std::cout << "\nlt(f_a) (lex) = ";
        print_term(lt_f_a_lex.first, lt_f_a_lex.second, vars_a);
        std::cout << "\n";

        std::cout << "lm(f_a) (grlex) = ";
        print_monomial(lm_f_a_grlex, vars_a);
        std::cout << "\nlc(f_a) (grlex) = " << lc_f_a_grlex;
        std::cout << "\nlt(f_a) (grlex) = ";
        print_term(lt_f_a_grlex.first, lt_f_a_grlex.second, vars_a);
        std::cout << "\n";

        std::cout << "lm(f_a) (grevlex) = ";
        print_monomial(lm_f_a_grevlex, vars_a);
        std::cout << "\nlc(f_a) (grevlex) = " << lc_f_a_grevlex;
        std::cout << "\nlt(f_a) (grevlex) = ";
        print_term(lt_f_a_grevlex.first, lt_f_a_grevlex.second, vars_a);
        std::cout << "\n\n";


        // b part
        std::cout << "====== (b) ======" << std::endl;
        PolyTrie<double> f_b({"x", "y", "z"});
        f_b.add_term({2, 8, 0}, 2.0);
        f_b.add_term({5, 1, 4}, -3.0);
        f_b.add_term({1 ,1,3}, 1.0);
        f_b.add_term({1, 4, 0}, -1.0);

        std::vector<std::string> vars_b = {"x","y","z"};

        auto multideg_f_b_lex = f_b.find_multideg(order::Lex{});
        auto multideg_f_b_grlex = f_b.find_multideg(order::GrLex{});
        auto multideg_f_b_grevlex = f_b.find_multideg(order::GrevLex{});

        auto lm_f_b_lex = f_b.leading_monomial(order::Lex{});
        auto lc_f_b_lex = f_b.leading_coeff(order::Lex{});
        auto lt_f_b_lex = f_b.leading_term(order::Lex{});

        auto lm_f_b_grlex = f_b.leading_monomial(order::GrLex{});
        auto lc_f_b_grlex = f_b.leading_coeff(order::GrLex{});
        auto lt_f_b_grlex = f_b.leading_term(order::GrLex{});

        auto lm_f_b_grevlex = f_b.leading_monomial(order::GrevLex{});
        auto lc_f_b_grevlex = f_b.leading_coeff(order::GrevLex{});
        auto lt_f_b_grevlex = f_b.leading_term(order::GrevLex{});

        std::cout << "multideg(f_b) (lex) = ";
        for (auto d : multideg_f_b_lex) std::cout << d << " ";
        std::cout << "\nmultideg(f_b) (grlex) = ";
        for (auto d : multideg_f_b_grlex) std::cout << d << " ";
        std::cout << "\nmultideg(f_b) (grevlex) = ";
        for (auto d : multideg_f_b_grevlex) std::cout << d << " ";
        std::cout << "\n";

        std::cout << "lm(f_b) (lex) = ";
        print_monomial(lm_f_b_lex, vars_b);
        std::cout << "\nlc(f_b) (lex) = " << lc_f_b_lex;
        std::cout << "\nlt(f_b) (lex) = ";
        print_term(lt_f_b_lex.first, lt_f_b_lex.second, vars_b);
        std::cout << "\n";

        std::cout << "lm(f_b) (grlex) = ";
        print_monomial(lm_f_b_grlex, vars_b);
        std::cout << "\nlc(f_b) (grlex) = " << lc_f_b_grlex;
        std::cout << "\nlt(f_b) (grlex) = ";
        print_term(lt_f_b_grlex.first, lt_f_b_grlex.second, vars_b);
        std::cout << "\n";

        std::cout << "lm(f_b) (grevlex) = ";
        print_monomial(lm_f_b_grevlex, vars_b);
        std::cout << "\nlc(f_b) (grevlex) = " << lc_f_b_grevlex;
        std::cout << "\nlt(f_b) (grevlex) = ";
        print_term(lt_f_b_grevlex.first, lt_f_b_grevlex.second, vars_b);
        std::cout << "\n";

        std::cout << std::endl << std::endl;
    } 
        std::cout << "=============== TASK 3 ===============" << std::endl << std::endl;
        std::cout << "z > y > x"<< std::endl;
        std::cout << "====== (a) ======" << std::endl;
    
    PolyTrie<double> f_a({"z", "y", "x"});
    f_a.add_term({1, 0, 0}, 2.0);
    f_a.add_term({0, 1, 0}, 3.0);
    f_a.add_term({0 ,0 ,1}, 1.0);
    f_a.add_term({2, 0, 0}, 1.0);
    f_a.add_term({0, 0, 2}, -1.0);
    f_a.add_term({3, 0, 0}, 1.0);

    std::vector<std::string> vars_a = {"z","y","x"};

    auto multideg_f_a_lex = f_a.find_multideg(order::Lex{});
    auto multideg_f_a_grlex = f_a.find_multideg(order::GrLex{});
    auto multideg_f_a_grevlex = f_a.find_multideg(order::GrevLex{});

    auto lm_f_a_lex = f_a.leading_monomial(order::Lex{});
    auto lc_f_a_lex = f_a.leading_coeff(order::Lex{});
    auto lt_f_a_lex = f_a.leading_term(order::Lex{});

    auto lm_f_a_grlex = f_a.leading_monomial(order::GrLex{});
    auto lc_f_a_grlex = f_a.leading_coeff(order::GrLex{});
    auto lt_f_a_grlex = f_a.leading_term(order::GrLex{});

    auto lm_f_a_grevlex = f_a.leading_monomial(order::GrevLex{});
    auto lc_f_a_grevlex = f_a.leading_coeff(order::GrevLex{});
    auto lt_f_a_grevlex = f_a.leading_term(order::GrevLex{});

    std::cout << "multideg(f_a) (lex) = ";
    for (auto d : multideg_f_a_lex) std::cout << d << " ";
    std::cout << "\nmultideg(f_a) (grlex) = ";
    for (auto d : multideg_f_a_grlex) std::cout << d << " ";
    std::cout << "\nmultideg(f_a) (grevlex) = ";
    for (auto d : multideg_f_a_grevlex) std::cout << d << " ";
    std::cout << "\n";

    std::cout << "lm(f_a) (lex) = ";
    print_monomial(lm_f_a_lex, vars_a);
    std::cout << "\nlc(f_a) (lex) = " << lc_f_a_lex;
    std::cout << "\nlt(f_a) (lex) = ";
    print_term(lt_f_a_lex.first, lt_f_a_lex.second, vars_a);
    std::cout << "\n";

    std::cout << "lm(f_a) (grlex) = ";
    print_monomial(lm_f_a_grlex, vars_a);
    std::cout << "\nlc(f_a) (grlex) = " << lc_f_a_grlex;
    std::cout << "\nlt(f_a) (grlex) = ";
    print_term(lt_f_a_grlex.first, lt_f_a_grlex.second, vars_a);
    std::cout << "\n";

    std::cout << "lm(f_a) (grevlex) = ";
    print_monomial(lm_f_a_grevlex, vars_a);
    std::cout << "\nlc(f_a) (grevlex) = " << lc_f_a_grevlex;
    std::cout << "\nlt(f_a) (grevlex) = ";
    print_term(lt_f_a_grevlex.first, lt_f_a_grevlex.second, vars_a);
    std::cout << "\n\n";


    // b part
    std::cout << "====== (b) ======" << std::endl;
    PolyTrie<double> f_b({"z", "y", "x"});
    f_b.add_term({2, 8, 0}, 2.0);
    f_b.add_term({5, 1, 4}, -3.0);
    f_b.add_term({1 ,1,3}, 1.0);
    f_b.add_term({1, 4, 0}, -1.0);

    std::vector<std::string> vars_b = {"z","y","x"};

    auto multideg_f_b_lex = f_b.find_multideg(order::Lex{});
    auto multideg_f_b_grlex = f_b.find_multideg(order::GrLex{});
    auto multideg_f_b_grevlex = f_b.find_multideg(order::GrevLex{});

    auto lm_f_b_lex = f_b.leading_monomial(order::Lex{});
    auto lc_f_b_lex = f_b.leading_coeff(order::Lex{});
    auto lt_f_b_lex = f_b.leading_term(order::Lex{});

    auto lm_f_b_grlex = f_b.leading_monomial(order::GrLex{});
    auto lc_f_b_grlex = f_b.leading_coeff(order::GrLex{});
    auto lt_f_b_grlex = f_b.leading_term(order::GrLex{});

    auto lm_f_b_grevlex = f_b.leading_monomial(order::GrevLex{});
    auto lc_f_b_grevlex = f_b.leading_coeff(order::GrevLex{});
    auto lt_f_b_grevlex = f_b.leading_term(order::GrevLex{});

    std::cout << "multideg(f_b) (lex) = ";
    for (auto d : multideg_f_b_lex) std::cout << d << " ";
    std::cout << "\nmultideg(f_b) (grlex) = ";
    for (auto d : multideg_f_b_grlex) std::cout << d << " ";
    std::cout << "\nmultideg(f_b) (grevlex) = ";
    for (auto d : multideg_f_b_grevlex) std::cout << d << " ";
    std::cout << "\n";

    std::cout << "lm(f_b) (lex) = ";
    print_monomial(lm_f_b_lex, vars_b);
    std::cout << "\nlc(f_b) (lex) = " << lc_f_b_lex;
    std::cout << "\nlt(f_b) (lex) = ";
    print_term(lt_f_b_lex.first, lt_f_b_lex.second, vars_b);
    std::cout << "\n";

    std::cout << "lm(f_b) (grlex) = ";
    print_monomial(lm_f_b_grlex, vars_b);
    std::cout << "\nlc(f_b) (grlex) = " << lc_f_b_grlex;
    std::cout << "\nlt(f_b) (grlex) = ";
    print_term(lt_f_b_grlex.first, lt_f_b_grlex.second, vars_b);
    std::cout << "\n";

    std::cout << "lm(f_b) (grevlex) = ";
    print_monomial(lm_f_b_grevlex, vars_b);
    std::cout << "\nlc(f_b) (grevlex) = " << lc_f_b_grevlex;
    std::cout << "\nlt(f_b) (grevlex) = ";
    print_term(lt_f_b_grevlex.first, lt_f_b_grevlex.second, vars_b);
    std::cout << "\n";

    std::cout << std::endl << std::endl; 


    return 0;
}