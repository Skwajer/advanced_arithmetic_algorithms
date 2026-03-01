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
        
        Polynomial f2({1});            // 1
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

    // Тест 1: Простой многочлен в обычном базисе
    // P(x) = 1 + 2x + 3x^2
    // В факториальном базисе: x^2 = x(x-1) + x, поэтому
    // P(x) = 1 + 2x + 3(x(x-1) + x) = 1 + 5x + 3x(x-1)
    // Коэффициенты: [1, 5, 3] в порядке [a0, a1, a2]
    {
        Polynomial p({1.0, 5.0, 3.0});  // already in factorial basis
        
        double x = 2.0;
        double expected = 1.0 + 5.0*2.0 + 3.0*2.0*1.0;  // = 1 + 10 + 6 = 17
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 1: P(x) = 1 + 5x + 3x(x-1)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }
    
        // Тест 4: Многочлен с отрицательными коэффициентами
    // P(x) = 2 - 3x + 4x(x-1) - 5x(x-1)(x-2)
    // Коэффициенты: [2, -3, 4, -5]
    {
        Polynomial p({2.0, -3.0, 4.0, -5.0});
        
        double x = 3.0;
        // Вычисляем:
        // x = 3
        // x(x-1) = 3*2 = 6
        // x(x-1)(x-2) = 6*1 = 6
        // P(3) = 2 + (-3)*3 + 4*6 + (-5)*6 = 2 - 9 + 24 - 30 = -13
        double expected = -13.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 4: P(x) = 2 - 3x + 4x(x-1) - 5x(x-1)(x-2)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 5: Дробные коэффициенты и дробный x
    // P(x) = 0.5 + 1.5x + 2.5x(x-1)
    // Коэффициенты: [0.5, 1.5, 2.5]
    {
        Polynomial p({0.5, 1.5, 2.5});
        
        double x = 1.5;
        // Вычисляем:
        // x = 1.5
        // x(x-1) = 1.5 * 0.5 = 0.75
        // P(1.5) = 0.5 + 1.5*1.5 + 2.5*0.75 = 0.5 + 2.25 + 1.875 = 4.625
        double expected = 4.625;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 5: Дробные коэффициенты\n";
        std::cout << "P(x) = 0.5 + 1.5x + 2.5x(x-1)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 6: Многочлен высокой степени
    // P(x) = 1 + x + x(x-1) + x(x-1)(x-2) + x(x-1)(x-2)(x-3)
    // Коэффициенты: [1, 1, 1, 1, 1]
    {
        Polynomial p({1.0, 1.0, 1.0, 1.0, 1.0});
        
        double x = 4.0;
        // Вычисляем:
        // x = 4
        // x(x-1) = 4*3 = 12
        // x(x-1)(x-2) = 12*2 = 24
        // x(x-1)(x-2)(x-3) = 24*1 = 24
        // P(4) = 1 + 4 + 12 + 24 + 24 = 65
        double expected = 65.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 6: Многочлен степени 4\n";
        std::cout << "P(x) = 1 + x + x(x-1) + x(x-1)(x-2) + x(x-1)(x-2)(x-3)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 7: Нулевой x
    // P(x) = 3 + 2x + 4x(x-1)
    // Коэффициенты: [3, 2, 4]
    {
        Polynomial p({3.0, 2.0, 4.0});
        
        double x = 0.0;
        // Вычисляем:
        // x = 0
        // x(x-1) = 0 * (-1) = 0
        // P(0) = 3 + 2*0 + 4*0 = 3
        double expected = 3.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 7: x = 0\n";
        std::cout << "P(x) = 3 + 2x + 4x(x-1)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 8: Отрицательный x
    // P(x) = 1 + 2x + 3x(x-1)
    // Коэффициенты: [1, 2, 3]
    {
        Polynomial p({1.0, 2.0, 3.0});
        
        double x = -2.0;
        // Вычисляем:
        // x = -2
        // x(x-1) = (-2) * (-3) = 6
        // P(-2) = 1 + 2*(-2) + 3*6 = 1 - 4 + 18 = 15
        double expected = 15.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 8: Отрицательный x\n";
        std::cout << "P(x) = 1 + 2x + 3x(x-1)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 9: Случай, когда x меньше степени
    // P(x) = 1 + 2x + 3x(x-1) + 4x(x-1)(x-2)
    // Коэффициенты: [1, 2, 3, 4]
    {
        Polynomial p({1.0, 2.0, 3.0, 4.0});
        
        double x = 1.0;
        // Вычисляем:
        // x = 1
        // x(x-1) = 1*0 = 0
        // x(x-1)(x-2) = 0*(-1) = 0
        // P(1) = 1 + 2*1 + 3*0 + 4*0 = 3
        double expected = 3.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 9: x меньше степени (x = 1, степень 3)\n";
        std::cout << "P(x) = 1 + 2x + 3x(x-1) + 4x(x-1)(x-2)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 10: Константный многочлен
    // P(x) = 42
    // Коэффициенты: [42]
    {
        Polynomial p({42.0});
        
        double x = 5.0;
        double expected = 42.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 10: Константный многочлен\n";
        std::cout << "P(x) = 42\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n\n";
    }

    // Тест 11: Ещё один тест с x равным индексу
    // P(x) = 1 + x + x(x-1) + x(x-1)(x-2)
    // Коэффициенты: [1, 1, 1, 1]
    {
        Polynomial p({1.0, 1.0, 1.0, 1.0});
        
        double x = 2.0;
        // Вычисляем:
        // x = 2
        // x(x-1) = 2*1 = 2
        // x(x-1)(x-2) = 2*1*0 = 0
        // P(2) = 1 + 2 + 2 + 0 = 5
        double expected = 5.0;
        double result = p.evaluateUsingFactorialPowers(x);
        
        std::cout << "Тест 11: x = 2 (равно индексу)\n";
        std::cout << "P(x) = 1 + x + x(x-1) + x(x-1)(x-2)\n";
        std::cout << "x = " << x << "\n";
        std::cout << "Ожидаемое: " << expected << "\n";
        std::cout << "Результат: " << result << "\n";
        std::cout << (std::abs(result - expected) < 1e-10 ? "ПРОЙДЕН" : "НЕ ПРОЙДЕН") << "\n";
    }
    return 0;
}