#include "../Numerical_solvers/Newtons_method.hpp"
#include "../Numerical_solvers/TheDichotomy_method.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <string>

struct RootInfo {
    std::string name;
    double a, b;
    double x0;
    std::function<double(double)> f;
};

void printSeparator() {
    std::cout << "\n" << std::string(60, '=') << "\n";
}

void compareMethods(const RootInfo& root, double eps) {
    printSeparator();
    std::cout << "РЕШЕНИЕ ДЛЯ: " << root.name << "\n";
    std::cout << "Точность: " << eps << "\n";
    printSeparator();
    
    std::cout << "\n--- МЕТОД ДИХОТОМИИ ---\n";
    try {
        double root_bisection = Dichotomy_method(root.a, root.b, 
                                const_cast<std::function<double(double)>&>(root.f), eps);
        std::cout << "Найденный корень (дихотомия): " << std::setprecision(10) << root_bisection << "\n";
    } catch (const std::runtime_error& e) {
        std::cout << "Дихотомия не применима: " << e.what() << "\n";
        std::cout << "Для этого уравнения используем только метод Ньютона.\n";
    }
    
    std::cout << "\n--- МЕТОД НЬЮТОНА ---\n";
    auto deriv = [&root](double x) {
        return computeDeravative(x, const_cast<std::function<double(double)>&>(root.f));
    };
    
    double root_newton = NewtonsMethod(root.x0, deriv, 
                                        const_cast<std::function<double(double)>&>(root.f), eps);
    std::cout << "Найденный корень (Ньютон): " << std::setprecision(10) << root_newton << "\n";
    
    std::cout << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    
    int n = 8;
    double eps = std::pow(10.0, -n);
    
    std::cout << "========================================================\n";
    std::cout << "РЕШЕНИЕ УРАВНЕНИЙ МЕТОДАМИ ДИХОТОМИИ И НЬЮТОНА\n";
    std::cout << "Точность: 10^-" << n << " = " << eps << "\n";
    std::cout << "========================================================\n";
    
    std::vector<RootInfo> roots;
    

    int n_pow = 3;
    double a_val = 8.0;
    double root_exact = std::pow(a_val, 1.0/n_pow);
    
    roots.push_back({
        "b) x^3 = 8 (корень x ≈ 2)",
        root_exact - 0.5, root_exact + 0.5,
        root_exact,
        [n_pow, a_val](double x) { return std::pow(x, n_pow) - a_val; }
    });
    

    auto func_c = [](double x) { return std::sqrt(1 - x*x) - std::exp(x) + 0.1; };
    
    roots.push_back({
        "c) √(1-x²) - e^x + 0.1 = 0 (корень 1, отриц.)",
        -0.99, -0.95,
        -0.97,
        func_c
    });
    
    roots.push_back({
        "c) √(1-x²) - e^x + 0.1 = 0 (корень 2, полож.)",
        0.05, 0.1,
        0.075,
        func_c
    });
    

    auto func_d = [](double x) { 
        double x3 = x*x*x;
        return x3*x3 - 5*x3 - 2;  // x^6 - 5x^3 - 2
    };
    
    roots.push_back({
        "d) x^6 = 5x^3 + 2 (корень 1, полож.)",
        1.7, 1.8,
        1.75,
        func_d
    });
    
    roots.push_back({
        "d) x^6 = 5x^3 + 2 (корень 2, отриц.)",
        -0.8, -0.7,
        -0.75,
        func_d
    });
    

    auto func_e = [](double x) { 
        return std::log2(x) - 1.0/(1.0 + x*x); 
    };
    
    roots.push_back({
        "e) log₂(x) = 1/(1+x²)",
        1.29, 1.3,
        1.295,
        func_e
    });
    

    auto func_f = [](double x) { return std::sin(x/2.0) - 1.0; };
    
    roots.push_back({
        "f) sin(x/2) = 1 (корень 1, π)",
        3.1, 3.2,
        3.14,
        func_f
    });
    
    roots.push_back({
        "f) sin(x/2) = 1 (корень 2, 5π)",
        15.6, 15.8,
        15.7,
        func_f
    });
    
    roots.push_back({
        "f) sin(x/2) = 1 (корень 3, 9π)",
        28.2, 28.4,
        28.27,
        func_f
    });
    

    auto func_g = [](double x) { return std::log(x) - 1.0; };
    
    roots.push_back({
        "g) ln(x) = 1",
        2.7, 2.8,
        2.71828,
        func_g
    });
    
    for (const auto& root : roots) {
        compareMethods(root, eps);
    }
    
    printSeparator();
    std::cout << "ВСЕ КОРНИ НАЙДЕНЫ!\n";
    printSeparator();
    
    return 0;
}