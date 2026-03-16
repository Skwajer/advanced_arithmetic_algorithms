#include <iostream>
#include <iterator>
#include "../Numerical_solvers/SImple_iteration.hpp"

int main()
{
        try {
        std::cout << "=== Part (a): x^3 + 3x^2 - 1 = 0 ===\n";
        auto eq1 = [](double x) { return x*x*x + 3*x*x - 1; };
        
        double root1 = simpleIteration(0.5, eq1);
        std::cout << "root 1 (≈0.53): " << root1 << "\n\n";
        
        double root2 = simpleIteration(-2.9, eq1);
        std::cout << "root 2 (≈-2.88): " << root2 << "\n\n";
        
        double root3 = simpleIteration(-0.6, eq1);
        std::cout << "root 3 (≈-0.65): " << root3 << "\n\n";
        
        std::cout << "=== Part (b): x^4 - x^3 = 0 ===\n";
        auto eq2 = [](double x) { return x*x*x*x - x*x*x; };
        
        double root4 = simpleIteration(0.5, eq2);
        std::cout << "root x=0: " << root4 << "\n\n";
        
        double root5 = simpleIteration(1.5, eq2);
        std::cout << "root x=1: " << root5 << "\n\n";
        
        std::cout << "=== Part (c): x^2 - 3x + 2 = 0 ===\n";
        auto eq3 = [](double x) { return x*x - 3*x + 2; };
        
        double root6 = simpleIteration(0.5, eq3);
        std::cout << "root x=1: " << root6 << "\n\n";
        
        double root7 = simpleIteration(2.5, eq3);
        std::cout << "root x=2: " << root7 << "\n\n";

        std::cout << "RESULTS" << std::endl;
        std::cout << "=== Part (a): x^3 + 3x^2 - 1 = 0 ===\n";
        std::cout << "root 1 (≈0.53): " << root1 << "\n\n";
        std::cout << "root 2 (≈-2.88): " << root2 << "\n\n";
        std::cout << "root 3 (≈-0.65): " << root3 << "\n\n";

        std::cout << "=== Part (b): x^4 - x^3 = 0 ===\n";
        
        std::cout << "root x=0: " << root4 << "\n\n";
        
        std::cout << "root x=1: " << root5 << "\n\n";

        std::cout << "=== Part (c): x^2 - 3x + 2 = 0 ===\n";
        
        std::cout << "root x=1: " << root6 << "\n\n";
        
        std::cout << "root x=2: " << root7 << "\n\n";
        


        
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        
    }
    
    return 0;
}
