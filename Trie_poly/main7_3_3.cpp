#include "trie_poly.hpp"
#include <iostream>
#include <vector>

int main() 
{
    using namespace order;

    std::cout << "exercise 2.6.8 from Кокс book" << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({0, 1, 0}, 1); // y
        f1.add_term({2, 0, 0}, -1); // -x^2

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({0, 0, 1}, 1); // -x^2
        f2.add_term({3, 0, 0}, -1); // -x^2
        
        basis.push_back(f1);
        basis.push_back(f2);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {x^2, y^2} is Groebner basis? " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "exercise 2.6.9 from Кокс book" << std::endl << std::endl;
    std::cout << "a)" << std::endl;    
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({2, 0, 0}, 1);
        f1.add_term({0, 1, 0}, -1); 

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({3, 0, 0}, 1); 
        f2.add_term({0, 0, 1}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, GrLex{});
        std::cout << "Basis {x^2 - y, x^3 - z} is Groebner basis? (GrLex) " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "b)" << std::endl;  
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({2, 0, 0}, 1);
        f1.add_term({0, 1, 0}, -1); 

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({3, 0, 0}, 1); 
        f2.add_term({0, 0, 1}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, InvLex{});
        std::cout << "Basis {x^2 - y, x^3 - z} is Groebner basis? (InvLex) " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "c)" << std::endl;  
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({1, 2, 0}, 1);
        f1.add_term({1, 0, 1}, -1);
        f1.add_term({0, 1, 0}, 1);

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({1, 1, 0}, 1); 
        f2.add_term({0, 0, 2}, -1);

        PolyTrie<double> f3({"x", "y", "z"});
        f3.add_term({1, 0, 0}, 1); 
        f3.add_term({0, 1, 4}, -1);
        
        basis.push_back(f1);
        basis.push_back(f2);
        basis.push_back(f3);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {x^2 - xy + y, xy - z^2, x - yz^4} is Groebner basis? (GrLex) " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "exercise 21.21 from Von zur Gathen Modern Computer Algebra:" << std::endl << std::endl;
    std::cout << "i)" << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y"});
        f1.add_term({1, 0}, 1);
        f1.add_term({0, 1}, 1); 

        PolyTrie<double> f2({"x", "y"});
        f2.add_term({0, 2}, 1); 
        f2.add_term({0, 0}, -1); 
        
        basis.push_back(f1);
        basis.push_back(f2);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {x + y, y^2 - 1} is Groebner basis? " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "ii)" << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"y", "x"});
        f1.add_term({1, 0}, 1);
        f1.add_term({0, 1}, 1); 

        PolyTrie<double> f2({"y", "x"});
        f2.add_term({2, 0}, 1); 
        f2.add_term({0, 0}, -1); 
        
        basis.push_back(f1);
        basis.push_back(f2);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {y + x, y^2 - 1} is Groebner basis? " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }

    std::cout << "iii)" << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y"});
        f1.add_term({2, 0}, 1);
        f1.add_term({0, 2}, 1);
        f1.add_term({0, 0}, -1);

        PolyTrie<double> f2({"x", "y"});
        f2.add_term({1, 1}, 1); 
        f2.add_term({0, 0}, -1);

        PolyTrie<double> f3({"x", "y",});
        f3.add_term({1, 0}, 1); 
        f3.add_term({0, 3}, 1);
        f3.add_term({0, 1}, -1);

        basis.push_back(f1);
        basis.push_back(f2);
        basis.push_back(f3);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {x^2 + y^2 - 1, xy - 1, x + y^3 - y} is Groebner basis? " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }
    std::cout << "iv)" << std::endl;
    {
        std::vector<PolyTrie<double>> basis;
        
        PolyTrie<double> f1({"x", "y", "z"});
        f1.add_term({1, 1, 1}, 1);
        f1.add_term({0, 0, 0}, -1);

        PolyTrie<double> f2({"x", "y", "z"});
        f2.add_term({1, 0, 0}, 1); 
        f2.add_term({0, 1, 0}, -1);

        PolyTrie<double> f3({"x", "y", "z"});
        f3.add_term({0, 2, 1}, 1); 
        f3.add_term({0, 0, 0}, -1);

        basis.push_back(f1);
        basis.push_back(f2);
        basis.push_back(f3);
        
        bool isGB = f1.PolyTrie<double>::isGroebnerBasis(basis, Lex{});
        std::cout << "Basis {xyz - 1, x - y, y^2z - 1} is Groebner basis? " << (isGB ? "Yes" : "No") << std::endl;
        std::cout << std::endl;  
    }
    

}