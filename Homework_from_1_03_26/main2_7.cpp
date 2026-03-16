#include "../Numerical_solvers/Newtons_method.hpp"
#include <cmath>

int main() {
    std::cout << std::fixed;
    std::cout << "=============================================================\n";
    std::cout << "ЗАДАЧА 7: МЕТОД НЬЮТОНА ДЛЯ КРАТНЫХ КОРНЕЙ\n";
    std::cout << "Уравнение: (x-1)^3 * sin(πx) * (cos(2πx)-1) = 0\n";
    std::cout << "Демонстрация линейной сходимости со знаменателем 1 - 1/p\n";
    std::cout << "=============================================================\n";
    
    double eps = 1e-8;
    int maxIter = 50;
    
    std::vector<double> roots = {
        0.2,
        1.2,
        2.2
    };
    
    auto f = [](double x) {return (x-1)*(x-1)*(x-1) * sin(M_PI*x) * (cos(2*M_PI*x) - 1);};
    std::function<double(double)> func = f;

    for (auto root : roots) 
    {
        auto deriv = [&func](double x) 
        {
            return computeDeravative(x, func);
        };
        NewtonsMethod(root, deriv, func, eps);
        std::cout << "------------------------\n";
    }
    
    return 0;
}