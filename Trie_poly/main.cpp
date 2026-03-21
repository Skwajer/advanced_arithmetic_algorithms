#include "trie_poly.hpp"
#include <complex>
#include <vector>

using complex = std::complex<double>;

int main()
{
    PolyTrie<double> f1_R3({"alpha", "betta", "gamma"});
    f1_R3.add_term({2, 3, 7}, -3.5);
    f1_R3.add_term({1, 4, -4}, 5);
    f1_R3.add_term({0, 0, 0}, -11.3);

    PolyTrie<double> f2_R3({"alpha", "betta", "gamma"});
    f2_R3.add_term({2, 3, 7}, -3.5);
    f2_R3.add_term({1, 4, -4}, 5);
    f2_R3.add_term({0, 0, 0}, -11.3);
    f2_R3.add_term({3, 3, 3}, -25.0);

    f1_R3 += f2_R3;

    f1_R3.print();
    std::vector<double> point_R3 = {5, 7, 2};
    std::cout << "value of f1_R3 in point{5, 7, 2} = " << f1_R3.evaluate(point_R3) << std::endl;

    PolyTrie<double> P({"x", "y"});
    
    // P(x,y) = 3*x^2*y + 2*x*y^2 - 5*x + 4*y - 7
    P.add_term({2, 1}, 3.0);  
    P.add_term({1, 2}, 2.0);
    P.add_term({1, 0}, -5.0); 
    P.add_term({0, 1}, 4.0);  
    P.add_term({0, 0}, -7.0);  
    
    std::cout << "P(x,y) = ";
    P.print();
    std::cout << std::endl;
    
    // Тест 1: точка (x=1, y=1)
    // 3*1^2*1 + 2*1*1^2 - 5*1 + 4*1 - 7 = 3 + 2 - 5 + 4 - 7 = -3
    std::vector<double> point1 = {1.0, 1.0};
    double val1 = P.evaluate(point1);
    std::cout << "P(1, 1) = " << val1 << " (expected: -3)" << std::endl;
    
    // Тест 2: точка (x=2, y=3)
    // 3*4*3 + 2*2*9 - 5*2 + 4*3 - 7 = 36 + 36 - 10 + 12 - 7 = 67
    std::vector<double> point2 = {2.0, 3.0};
    double val2 = P.evaluate(point2);
    std::cout << "P(2, 3) = " << val2 << " (expected: 67)" << std::endl;
    
    // Тест 3: точка (x=0, y=0)
    // -7
    std::vector<double> point3 = {0.0, 0.0};
    double val3 = P.evaluate(point3);
    std::cout << "P(0, 0) = " << val3 << " (expected: -7)" << std::endl;
    
    // Тест 4: точка (x=0, y=5)
    // 4*5 - 7 = 20 - 7 = 13
    std::vector<double> point4 = {0.0, 5.0};
    double val4 = P.evaluate(point4);
    std::cout << "P(0, 5) = " << val4 << " (expected: 13)" << std::endl;
    
    // Тест 5: точка (x=3, y=0)
    // -5*3 - 7 = -15 - 7 = -22
    std::vector<double> point5 = {3.0, 0.0};
    double val5 = P.evaluate(point5);
    std::cout << "P(3, 0) = " << val5 << " (expected: -22)" << std::endl;
    
    std::cout << std::endl;

    PolyTrie<double> f3_R3({"x", "y", "z"});
    f1_R3.add_term({2, 3, 7}, -3.5);
    f1_R3.add_term({2, 5, 7}, -3.5);
    f1_R3.add_term({2, 6, 2}, -3.5);
    f1_R3.add_term({1, 3, 6}, 10.5);
    PolyTrie<double> H = f1_R3.homogeneous_part(10);
    H.print();

    return 0;
}