#include <iostream>
#include <stdexcept>
#include <utility>



std::pair<double, double> divide(double a0, double a1, double b0, double b1)
{
    if (!b0 && !b1)
    {
        throw std::invalid_argument("denuminator cannot be zero");
    }

    auto m1 = a0*b0;
    auto m2 = a1*b1;
    auto m3 = (a0 + a1) * (b0 - b1);
    double m4 = b0*b0;
    double m5 = b1*b1;
    double den_sum1 = m4 + m5;
    double Im = (m3 - m1 + m2) / den_sum1;
    double Re = (m1 + m2) / den_sum1;
    return {Re, Im};
}
bool approx_equal(double a, double b, double epsilon = 1e-10)
{
    return std::abs(a - b) < epsilon;
}

void test_case(double a0, double a1, double b0, double b1, double expected_re, double expected_im)
{
    std::cout << "Testing (" << a0 << " + " << a1 << "i) / (" << b0 << " + " << b1 << "i)" << std::endl;
    
    try {
        auto result = divide(a0, a1, b0, b1);
        
        std::cout << "  Result:   (" << result.first << " + " << result.second << "i)" << std::endl;
        std::cout << "  Expected: (" << expected_re << " + " << expected_im << "i)" << std::endl;
        
        if (approx_equal(result.first, expected_re) && approx_equal(result.second, expected_im)) {
            std::cout << "  ✅ PASS" << std::endl;
        } else {
            std::cout << "  ❌ FAIL" << std::endl;
        }
        std::cout << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "  ❌ EXCEPTION: " << e.what() << std::endl << std::endl;
    }
}

int main()
{
    std::cout << "COMPLEX DIVISION TESTS\n";
    std::cout << "======================\n\n";
    
    // Test case 1: Simple real division
    test_case(6.0, 0.0, 2.0, 0.0, 3.0, 0.0);
    
    // Test case 2: Simple imaginary division
    test_case(0.0, 8.0, 0.0, 2.0, 4.0, 0.0);
    
    // Test case 3: (1 + i) / (1 + i) = 1
    test_case(1.0, 1.0, 1.0, 1.0, 1.0, 0.0);
    
    // Test case 4: (1 + 2i) / (3 + 4i) = 0.44 + 0.08i
    test_case(1.0, 2.0, 3.0, 4.0, 0.44, 0.08);
    
    // Test case 5: (3 - 2i) / (1 + 2i) = -0.2 - 1.6i
    test_case(3.0, -2.0, 1.0, 2.0, -0.2, -1.6);
    
    // Test case 6: (2 + 3i) / (2 - 3i) = (-5/13) + (12/13)i
    test_case(2.0, 3.0, 2.0, -3.0, -5.0/13.0, 12.0/13.0);
    
    // Test case 7: Zero denominator test
    std::cout << "Testing zero denominator case:\n";
    try {
        divide(1.0, 1.0, 0.0, 0.0);
        std::cout << "  ❌ FAIL - Should have thrown exception" << std::endl << std::endl;
    }
    catch (const std::invalid_argument& e) {
        std::cout << "  ✅ PASS - Correctly threw: " << e.what() << std::endl << std::endl;
    }
    
    // Test case 8: (3 + 4i) / (0 + 2i) = 2 - 1.5i
    test_case(3.0, 4.0, 0.0, 2.0, 2.0, -1.5);
    
    // Test case 9: (0 + 5i) / (3 + 4i) = 0.8 + 0.6i
    test_case(0.0, 5.0, 3.0, 4.0, 0.8, 0.6);
    
    // Test case 10: (7 + 0i) / (0 + 3i) = 0 - (7/3)i ≈ 0 - 2.33333i
    test_case(7.0, 0.0, 0.0, 3.0, 0.0, -7.0/3.0);
    
    // Test case 11: (-4 + 2i) / (1 - 3i) = (-10 - 10i)/10 = -1 - i
    test_case(-4.0, 2.0, 1.0, -3.0, -1.0, -1.0);
    
    // Test case 12: (1 + 0i) / (0 + 1i) = 0 - i
    test_case(1.0, 0.0, 0.0, 1.0, 0.0, -1.0);
    
    // Test case 13: (0 + 1i) / (1 + 0i) = 0 + i
    test_case(0.0, 1.0, 1.0, 0.0, 0.0, 1.0);
    
    // Test case 14: (5 + 5i) / (2 + 2i) = 2.5 + 0i
    test_case(5.0, 5.0, 2.0, 2.0, 2.5, 0.0);
    
    // Test case 15: (2 + 0i) / (1 + i) = 1 - i
    test_case(2.0, 0.0, 1.0, 1.0, 1.0, -1.0);
    
    // Test case 16: (0 + 2i) / (1 + i) = 1 + i
    test_case(0.0, 2.0, 1.0, 1.0, 1.0, 1.0);
    
    // Test case 17: (-3 - 3i) / (1 - 2i) = (3 - 9i)/5 = 0.6 - 1.8i
    test_case(-3.0, -3.0, 1.0, -2.0, 0.6, -1.8);
    
    // Test case 18: (4 - 2i) / (-1 + 3i) = (-10 - 10i)/10 = -1 - i
    test_case(4.0, -2.0, -1.0, 3.0, -1.0, -1.0);
    
    // Test case 19: (100 + 0i) / (0 + 10i) = 0 - 10i
    test_case(100.0, 0.0, 0.0, 10.0, 0.0, -10.0);
    
    // Test case 20: (1e6 + 1e6i) / (1e3 + 1e3i) = 1000 + 0i
    test_case(1000000.0, 1000000.0, 1000.0, 1000.0, 1000.0, 0.0);
    
    // Test case 21: (1e-6 + 1e-6i) / (1e-3 + 1e-3i) = 0.001 + 0i
    test_case(1e-6, 1e-6, 1e-3, 1e-3, 0.001, 0.0);
    
    // Test case 22: (1 + i) / (2 + 0i) = 0.5 + 0.5i
    test_case(1.0, 1.0, 2.0, 0.0, 0.5, 0.5);
    
    // Test case 23: (1 + 2i) / (0 + 4i) = 0.5 - 0.25i
    test_case(1.0, 2.0, 0.0, 4.0, 0.5, -0.25);
    
    // Test case 24: (3 + 4i) / (3 - 4i) = (-7 + 24i)/25 = -0.28 + 0.96i
    test_case(3.0, 4.0, 3.0, -4.0, -7.0/25.0, 24.0/25.0);
    
    // Test case 25: (1 + 0i) / (-1 + 0i) = -1 + 0i
    test_case(1.0, 0.0, -1.0, 0.0, -1.0, 0.0);
    
    // Test case 26: (0 + 1i) / (0 - 1i) = -1 + 0i
    test_case(0.0, 1.0, 0.0, -1.0, -1.0, 0.0);

    // Test case 28: (7 + 8i) / (9 + 10i) = (143 + 2i)/181 ≈ 0.79006 + 0.01105i
    test_case(7.0, 8.0, 9.0, 10.0, 143.0/181.0, 2.0/181.0);
    
    // Test case 29: (0.5 + 0.5i) / (0.5 + 0.5i) = 1 + 0i
    test_case(0.5, 0.5, 0.5, 0.5, 1.0, 0.0);
    
    // Test case 30: (π + ei) / (1 + i) where π≈3.14159, e≈2.71828
    test_case(3.14159, 2.71828, 1.0, 1.0, (3.14159+2.71828)/2.0, (2.71828-3.14159)/2.0);
    
    return 0;
}