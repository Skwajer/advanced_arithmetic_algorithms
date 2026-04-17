#include "trie_poly.hpp"
#include <iostream>
#include <vector>

int main() 
{
    using namespace order;

    std::cout << "EXERCISE 2.7.1 FROM KOKS BOOK" << std::endl << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y"});
        f1.add_term({3, 0}, 1);
        f1.add_term({1, 1}, -2);

        PolyTrie<double> f2({"x", "y"});
        f2.add_term({2, 1}, 1);
        f2.add_term({0, 2}, -2);
        f2.add_term({1, 0}, 1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        auto Grobner_basis =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        std::cout << "GRLEX order:" << std::endl;
        for (size_t i = 0; i < Grobner_basis.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";  
            Grobner_basis[i].print();
        }

        std::cout << std::endl;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "EXERCISE 2.7.2 FROM KOKS BOOK" << std::endl << std::endl;
    {
        std::cout << "a)" << "LEX order:" << std::endl;
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y"});
        f1.add_term({2, 1}, 1);
        f1.add_term({0, 0}, -1);

        PolyTrie<double> f2({"x", "y"});
        f2.add_term({1, 2}, 1);
        f2.add_term({1, 0}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        auto Grobner_basis_lex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, Lex{});
        for (size_t i = 0; i < Grobner_basis_lex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_lex[i].print();
        }

        std::cout << std::endl;

        std::cout << "a)" << "GRLEX order:" << std::endl;
        auto Grobner_basis_grlex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        for (size_t i = 0; i < Grobner_basis_grlex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_grlex[i].print();
        }

        std::cout << std::endl;
    }

    {
        std::cout << "b)" << "LEX order:" << std::endl;
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y"});
        f1.add_term({2, 0}, 1);
        f1.add_term({0, 1}, 1);

        PolyTrie<double> f2({"x", "y"});
        f2.add_term({4, 0}, 1);
        f2.add_term({2, 1}, 2);
        f2.add_term({0, 2}, 1);
        f2.add_term({0, 0}, 3);
        
        basis.push_back(f1);
        basis.push_back(f2);
        auto Grobner_basis_lex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, Lex{});
        for (size_t i = 0; i < Grobner_basis_lex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_lex[i].print();
        }

        std::cout << std::endl;

        std::cout << "b)" << "GRLEX order:" << std::endl;
        auto Grobner_basis_grlex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        for (size_t i = 0; i < Grobner_basis_grlex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_grlex[i].print();
        }

        std::cout << std::endl;
    }

    {
        std::cout << "c)" << "LEX order:" << std::endl;
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({1, 0, 0}, 1);
        f1.add_term({0, 0, 4}, -1);

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({0, 1, 0}, 1);
        f2.add_term({0, 0, 5}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        auto Grobner_basis_lex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, Lex{});
        for (size_t i = 0; i < Grobner_basis_lex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_lex[i].print();
        }

        std::cout << std::endl;

        std::cout << "c)" << "GRLEX order:" << std::endl;
        auto Grobner_basis_grlex =
            f1.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        for (size_t i = 0; i < Grobner_basis_grlex.size(); i++)
        {
            std::cout << "f" << i+1 << " = ";
            Grobner_basis_grlex[i].print();
        }

        std::cout << std::endl;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "EXERCISE 2.8.1 FROM KOKS BOOK" << std::endl << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({1, 3, 0}, 1);
        f.add_term({0, 0, 2}, -1);
        f.add_term({0, 5, 0}, 1);
        f.add_term({0, 0, 3}, -1);
        std::cout << "f = "; 
        f.print();

        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({3, 0, 0}, -1);
        f1.add_term({0, 1, 0}, 1);

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({2, 1, 0}, 1);
        f2.add_term({0, 0, 1}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        auto Grobner_basis_lex =
            f.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        std::cout << "LEX order:" << std::endl;
        std::cout << "I = \n    <" << std::endl;
        for (size_t i = 0; i < Grobner_basis_lex.size(); i++)
        {
            std::cout << "     f" << i+1 << " = ";  
            Grobner_basis_lex[i].print();
        }
        std::cout << "    >" << std::endl;

        auto [quotions, r] =  f.divide(Grobner_basis_lex, GrLex{});
        std:: cout << "f " << (r.isZero() ? "doesn't belong to I": "belong to I") << std::endl;

        auto Grobner_basis_grlex =
            f.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        std::cout << "\nGRLEX order:" << std::endl;
        std::cout << "I = \n    <" << std::endl;
        for (size_t i = 0; i < Grobner_basis_grlex.size(); i++)
        {
            std::cout << "     f" << i+1 << " = ";  
            Grobner_basis_grlex[i].print();
        }
        std::cout << "    >" << std::endl;

        auto [quotions_grlex, r_grlex] =  f.divide(Grobner_basis_grlex, GrLex{});
        std:: cout << "f " << (r.isZero() ? "doesn't belong to I": "belong to I") << std::endl;

        std::cout << std::endl;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "EXERCISE 2.8.2 FROM KOKS BOOK" << std::endl << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f({"x", "y", "z"});
        f.add_term({3, 0, 1}, 1);
        f.add_term({0, 2, 0}, -2);
        std::cout << "f = "; 
        f.print();

        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({1, 0, 1}, 1);
        f1.add_term({0, 1, 0}, -1);

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({1, 1, 0}, 1);
        f2.add_term({0, 0, 2}, 2);

        PolyTrie<double> f3({"x", "y", "z"});
        f3.add_term({0, 1, 0}, 1);
        f3.add_term({0, 0, 1}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        basis.push_back(f3);
        auto Grobner_basis_lex =
            f.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        std::cout << "LEX order:" << std::endl;
        std::cout << "I = \n    <" << std::endl;
        for (size_t i = 0; i < Grobner_basis_lex.size(); i++)
        {
            std::cout << "     f" << i+1 << " = ";  
            Grobner_basis_lex[i].print();
        }
        std::cout << "    >" << std::endl;

        auto [quotions, r] =  f.divide(Grobner_basis_lex, GrLex{});
        std:: cout << "f " << (r.isZero() ? "doesn't belong to I": "belong to I") << std::endl;

        auto Grobner_basis_grlex =
            f.PolyTrie<double>::build_GrobnerBasis(basis, GrLex{});
        std::cout << "\nGRLEX order:" << std::endl;
        std::cout << "I = \n    <" << std::endl;
        for (size_t i = 0; i < Grobner_basis_grlex.size(); i++)
        {
            std::cout << "     f" << i+1 << " = ";  
            Grobner_basis_grlex[i].print();
        }
        std::cout << "    >" << std::endl;

        auto [quotions_grlex, r_grlex] =  f.divide(Grobner_basis_grlex, GrLex{});
        std:: cout << "f " << (r.isZero() ? "doesn't belong to I": "belong to I") << std::endl;

        std::cout << std::endl;
    }
    
}