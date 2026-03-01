#include <iostream>
#include <vector>
#include <cmath>
#include "cartesian_sum.hpp"

int main() 
{
    
    // Test 1
    {
        std::vector<int> a = {1, 2, 3};
        std::vector<int> b = {2, 3, 4};
        
        auto result = cartesian_sum(a, b, 3);
        
        std::cout << "Test 1:\n";
        std::cout << "Result:\n";
        for (size_t i = 0; i < result.size(); i++) 
        {
            if (std::round(result[i].real()) > 0)
            {
                std::cout << "  sum " << i << ": " << std::round(result[i].real()) << "\n";
            }
        }
        std::cout << "Expected:\n";
        std::cout << "  sum 3: 1\n";
        std::cout << "  sum 4: 2\n";
        std::cout << "  sum 5: 3\n";
        std::cout << "  sum 6: 2\n";
        std::cout << "  sum 7: 1\n\n";
    }
    
    // Test 2
    {
        std::vector<int> a = {1, 1, 2};
        std::vector<int> b = {2, 3};
        auto result = cartesian_sum(a, b, 3);
        
        std::cout << "Test 2:\n";
        std::cout << "Result:\n";
        for (size_t i = 0; i < result.size(); i++) {
            if (std::round(result[i].real()) > 0) {
                std::cout << "  sum " << i << ": " << std::round(result[i].real()) << "\n";
            }
        }
        std::cout << "Expected:\n";
        std::cout << "  sum 3: 2\n";
        std::cout << "  sum 4: 3\n";
        std::cout << "  sum 5: 1\n\n";
    }
    
    // Test 3
    {
        std::vector<int> a = {0, 1, 2};
        std::vector<int> b = {0, 2, 4};
        auto result = cartesian_sum(a, b, 3);
        
        std::cout << "Test 3:\n";
        std::cout << "Result:\n";
        for (size_t i = 0; i < result.size(); i++) {
            if (std::round(result[i].real()) > 0) {
                std::cout << "  sum " << i << ": " << std::round(result[i].real()) << "\n";
            }
        }
        std::cout << "Expected:\n";
        std::cout << "  sum 0: 1\n";
        std::cout << "  sum 1: 1\n";
        std::cout << "  sum 2: 2\n";
        std::cout << "  sum 3: 1\n";
        std::cout << "  sum 4: 2\n";
        std::cout << "  sum 5: 1\n";
        std::cout << "  sum 6: 1\n\n";
    }
    
    // Test 4
    {
        std::vector<int> a = {5, 5, 5};
        std::vector<int> b = {5, 5};
        auto result = cartesian_sum(a, b, 5);
        
        std::cout << "Test 4:\n";
        std::cout << "Result:\n";
        for (size_t i = 0; i < result.size(); i++) {
            if (std::round(result[i].real()) > 0) {
                std::cout << "  sum " << i << ": " << std::round(result[i].real()) << "\n";
            }
        }
        std::cout << "Expected:\n";
        std::cout << "  sum 10: 6\n\n";
    }
    
    // Test 5
    {
        std::vector<int> a = {0, 20};
        std::vector<int> b = {0, 20};
        auto result = cartesian_sum(a, b, 2);
        
        std::cout << "Test 5:\n";
        std::cout << "Result:\n";
        for (size_t i = 0; i < result.size(); i++) {
            if (std::round(result[i].real()) > 0) {
                std::cout << "  sum " << i << ": " << std::round(result[i].real()) << "\n";
            }
        }
        std::cout << "Expected:\n";
        std::cout << "  sum 0: 1\n";
        std::cout << "  sum 20: 2\n";
        std::cout << "  sum 40: 1\n\n";
    }
    
    return 0;
}