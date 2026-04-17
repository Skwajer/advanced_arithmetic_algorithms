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
    
}