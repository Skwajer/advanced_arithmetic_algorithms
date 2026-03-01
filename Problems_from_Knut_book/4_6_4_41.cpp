#include <utility>
#include <complex>
#include <iostream>
#include <iomanip>

std::pair<double, double> mul_complex_3mul(double a, double b, double c, double d) {
    double p1 = a * c;
    double p2 = b * d;
    double p3 = (a + b) * (c + d);
    
    double real_part = p1 - p2;
    double imag_part = p3 - p1 - p2;
    
    return {real_part, imag_part};
}

void test_case(double a, double b, double c, double d) {
    std::pair<double, double> result = mul_complex_3mul(a, b, c, d);
    std::complex<double> expected = std::complex<double>(a, b) * std::complex<double>(c, d);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "(" << a << " + " << b << "i) * (" << c << " + " << d << "i) = ";
    std::cout << "(" << result.first << " + " << result.second << "i)";
    std::cout << "  [expected: (" << expected.real() << " + " << expected.imag() << "i)]";
    
    if (std::abs(result.first - expected.real()) < 1e-10 && 
        std::abs(result.second - expected.imag()) < 1e-10) {
        std::cout << " ✓";
    } else {
        std::cout << " ✗";
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "=== TESTS ===\n\n";
    
    test_case(3, 4, 2, 5);
    test_case(1, 0, 5, 7);
    test_case(0, 1, 3, 2);
    test_case(2, 3, -1, 4);
    test_case(-2, -3, 4, -5);
    test_case(0, 0, 10, 20);
    test_case(1.5, 2.5, 3.5, 4.5);
    test_case(100, 0, 0, 100);
    test_case(1, 1, 1, 1);
    test_case(1, -1, 1, -1);
    test_case(2, 0, 3, 0);
    test_case(0, 2, 0, 3);
    test_case(7, 8, 3, 2);
    test_case(10, 20, 30, 40);
    test_case(-5, 5, -5, 5);
    
    return 0;
}