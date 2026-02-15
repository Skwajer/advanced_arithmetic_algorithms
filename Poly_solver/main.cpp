#include "include/BelongingToLinearShell.h"

int main()
{
     Polynomial g1 = {1.0, 1.0, -2.0}; // 1 + x
     Polynomial g2 = {1.0, -1.0, 1.0}; // 1 - x
     Polynomial g3 = {1.0, -2.0, 2.0};

     Polynomial f = {2.0, 0.0, 2};

     auto result = isInLinearSpan(f, {g1, g2, g3});

     if (result.has_value())
     {
     std::cout << "f belongs to the linear shell\n";
     std::cout << "decomposition coefficients:\n";
     for (size_t i = 0; i < result->size(); ++i)
     {
     std::cout << "A" << i+1 << " = " << (*result)[i] << "\n";
     }
     }

     std::vector<double> coeffs{4, 12, 7, 0, 1};  // x^4 + 0x^3 + 7x^2 + 12x^1 + 4
    
     Polynomial poly(coeffs);
     poly.print_decomposition_to_degrees();
     poly.get_decomposition_to_degrees(-3);
     Polynomial shifted_poly = poly.get_shifted_representation(5);
     shifted_poly.print_decomposition_to_degrees();

     {
        std::cout << "=== Test 1: Simple limit ===\n";
        Polynomial f({1, 2, 1}, 0.0);     // x² + 2x + 1
        Polynomial g({-1, 0, 1}, 0.0);    // x² - 1
        double x = 2.0;
        
        double result = Polynomial::limit_at_endpoint(f, g, x);
        std::cout << "f(x) = "; f.print_decomposition_to_degrees();
        std::cout << "g(x) = "; g.print_decomposition_to_degrees();
        std::cout << "x → " << x << "\n";
        std::cout << "Result: " << result << "\n";
        std::cout << "Expected: 3.0\n\n";
    }
    
    // Тест 2: Знаменатель 0, числитель не 0 → бесконечность
     {
        std::cout << "=== Test 2: Denominator zero, numerator non-zero ===\n";
        Polynomial f({1, 0, 0}, 0.0);     // 1
        Polynomial g({0, 1, 0}, 0.0);    // x
        double x = 0.0;
    
        double result = Polynomial::limit_at_endpoint(f, g, x);
        std::cout << "f(x) = "; f.print_decomposition_to_degrees();
        std::cout << "g(x) = "; g.print_decomposition_to_degrees();
        std::cout << "x → " << x << "\n";
        std::cout << "Result: ";
        if (std::isinf(result)) {
            std::cout << (result > 0 ? "+∞" : "-∞");
        } else {
            std::cout << result;
        }
        std::cout << "\nExpected: +∞\n\n";
    }

    // Тест 6: Кратный корень, разная кратность
{
    std::cout << "=== Test 6: Multiple root, different multiplicity ===\n";
    Polynomial f({1, -3, 3, -1}, 0);    // (x-1)³
    Polynomial g({1, -2, 1}, 0);        // (x-1)²
    double x = 1.0;
    
    double result = Polynomial::limit_at_endpoint(f, g, x);
    std::cout << "Result: " << result << "\n";
    std::cout << "Expected: 0.0\n\n";
}

// Тест 7: Знаменатель обнуляется после сокращения
{
    std::cout << "=== Test 7: Denominator still zero after cancellation ===\n";
    Polynomial f({1, -1});           // x - 1
    Polynomial g({1, -2, 1});        // (x-1)²
    double x = 1.0;
    
    double result = Polynomial::limit_at_endpoint(f, g, x);
    std::cout << "Result: ";
    if (std::isinf(result)) {
        std::cout << (result > 0 ? "+∞" : "-∞");
    } else {
        std::cout << result;
    }
    std::cout << "\nExpected: +∞\n\n";
}

   std::cout << "=== T(x) = f1^k(s1(x)) / f2^l(s2(x)) limits ===\n\n";
    
    // Тест 1: f1(t)=t, s1(x)=x, k=2 → x
    //        f2(t)=1, s2(x)=x, l=1 → 1
    // T(x) = x, lim_{x→2} = 2
    {
        std::cout << "Test 1: Identity\n";
        Polynomial f1({0, 1});     // t
        Polynomial s1({0, 1});     // x
        int k = 2;
        
        Polynomial f2({1});        // 1
        Polynomial s2({0, 1});     // x
        int l = 1;
        
        double result = Polynomial::limit_T_at_point(f1, k, s1, f2, l, s2, 2.0);
        std::cout << "Result: " << result << "\n";
        std::cout << "Expected: 2.0\n\n";
    }

    // Тест 2: f1(t)=t+1, s1(x)=x, k=2 → x+2
    //        f2(t)=t, s2(x)=x, l=1 → x
    // T(x) = (x+2)/x, lim_{x→∞} = 1
    {
        std::cout << "Test 2: Infinity limit\n";
        Polynomial f1({1, 1});     // 1 + t
        Polynomial s1({0, 1});     // x
        int k = 2;
        
        Polynomial f2({0, 1});     // t
        Polynomial s2({0, 1});     // x
        int l = 1;
        
        double result = Polynomial::limit_T_at_infinity(f1, k, s1, f2, l, s2, true);
        std::cout << "Result: " << result << "\n";
        std::cout << "Expected: 1.0\n\n";
    }
    
    // Тест 3: f1(t)=t, s1(x)=x², k=2 → x²
    //        f2(t)=t, s2(x)=x, l=1 → x
    // T(x) = x²/x = x, lim_{x→0} = 0
    {
        std::cout << "Test 3: Quadratic inner\n";
        Polynomial f1({0, 1});     // t
        Polynomial s1({0, 0, 1});  // x²
        int k = 2;
        
        Polynomial f2({0, 1});     // t
        Polynomial s2({0, 1});     // x
        int l = 1;
        
        double result = Polynomial::limit_T_at_point(f1, k, s1, f2, l, s2, 0.0);
        std::cout << "Result: " << result << "\n";
        std::cout << "Expected: 0.0\n\n";
    }
    return 0;
}