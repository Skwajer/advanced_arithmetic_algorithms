#include "trie_poly.hpp"
#include <complex>
#include <vector>

using complex = std::complex<double>;

#include <complex>
#include <vector>

using complex = std::complex<double>;

void print_supp(const std::map<std::vector<int>, double>& supp, const std::vector<std::string>& var_names) {
    if (supp.empty()) {
        std::cout << "  (empty)";
        return;
    }
    
    bool first = true;
    for (const auto& [degrees, coeff] : supp) {
        if (!first) std::cout << " + ";
        first = false;
        
        std::cout << coeff;
        for (size_t i = 0; i < degrees.size(); ++i) {
            if (degrees[i] != 0) {
                std::cout << var_names[i];
                if (degrees[i] != 1) {
                    std::cout << "^" << degrees[i];
                }
            }
        }
    }
}

void test_get_supp() 
{
    std::cout << "Testing get_supp..." << std::endl;
    
    // Test 1: Empty polynomial
    PolyTrie<double> p1({"x", "y"});
    auto supp1 = p1.get_supp();
    std::cout << "  Empty polynomial: " << (supp1.empty() ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp1, {"x", "y"});
    std::cout << std::endl;
    
    // Test 2: Single term
    PolyTrie<double> p2({"x", "y"});
    p2.add_term({2, 1}, 3.5);
    auto supp2 = p2.get_supp();
    bool test2 = (supp2.size() == 1 && supp2.begin()->first == std::vector<int>{2, 1} && supp2.begin()->second == 3.5);
    std::cout << "  Single term: " << (test2 ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp2, {"x", "y"});
    std::cout << std::endl;
    
    // Test 3: Multiple terms
    PolyTrie<double> p3({"x", "y"});
    p3.add_term({2, 0}, 3.0);
    p3.add_term({0, 3}, 4.0);
    p3.add_term({1, 1}, 2.0);
    auto supp3 = p3.get_supp();
    bool test3 = (supp3.size() == 3 && 
                  supp3[{2,0}] == 3.0 && 
                  supp3[{0,3}] == 4.0 && 
                  supp3[{1,1}] == 2.0);
    std::cout << "  Multiple terms: " << (test3 ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp3, {"x", "y"});
    std::cout << std::endl;
    
    // Test 4: Zero coefficients are ignored
    PolyTrie<double> p4({"x", "y"});
    p4.add_term({2, 0}, 3.0);
    p4.add_term({1, 1}, 0.0);
    p4.add_term({0, 2}, 5.0);
    auto supp4 = p4.get_supp();
    bool test4 = (supp4.size() == 2 && 
                  supp4.find({1,1}) == supp4.end());
    std::cout << "  Zero coefficients ignored: " << (test4 ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp4, {"x", "y"});
    std::cout << std::endl;
    
    // Test 5: Negative coefficients
    PolyTrie<double> p5({"x", "y"});
    p5.add_term({2, 0}, -3.0);
    p5.add_term({0, 1}, -5.0);
    auto supp5 = p5.get_supp();
    bool test5 = (supp5.size() == 2 && 
                  supp5[{2,0}] == -3.0 && 
                  supp5[{0,1}] == -5.0);
    std::cout << "  Negative coefficients: " << (test5 ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp5, {"x", "y"});
    std::cout << std::endl;
    
    // Test 6: Large degrees
    PolyTrie<double> p6({"x", "y"});
    p6.add_term({10, 5}, 2.5);
    p6.add_term({0, 20}, 1.5);
    auto supp6 = p6.get_supp();
    bool test6 = (supp6.size() == 2 && 
                  supp6[{10,5}] == 2.5 && 
                  supp6[{0,20}] == 1.5);
    std::cout << "  Large degrees: " << (test6 ? "✓" : "✗");
    std::cout << " -> ";
    print_supp(supp6, {"x", "y"});
    std::cout << std::endl;
    
    std::cout << "Test complete!" << std::endl;
}

int main() 
{
    test_get_supp();
    return 0;
}