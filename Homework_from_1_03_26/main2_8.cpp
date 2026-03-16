#include <cmath>
#include <functional>
#include <iostream>
#include "../Numerical_solvers/Newtons_method.hpp"

int main() 
{
    std::cout << std::fixed;
    std::cout << "=============================================================\n";
    std::cout << "ЗАДАЧА 8: МЕТОД НЬЮТОНА ДЛЯ КРАТНЫХ КОРНЕЙ\n";
    std::cout << "Уравнение: (x-1)^3 * sin(πx) * (cos(2πx)-1) = 0\n";
    std::cout << "=============================================================\n";
    
    double eps = 1e-8;
    int maxIter = 50;
    
    std::vector<std::pair<double, size_t>> roots = {
        {0.2, 3},
        {1.2, 6},
        {2.2, 3}
    };
    
    auto f = [](double x) {return (x-1)*(x-1)*(x-1) * sin(M_PI*x) * (cos(2*M_PI*x) - 1);};
    std::function<double(double)> func = f;

    for (auto [root, multi] : roots) 
    {
        auto deriv = [&func](double x) 
        {
            return computeDeravative(x, func);
        };
        NewtonsMethod(root, deriv, func, eps, multi);
        std::cout << "------------------------\n";
    }
    
    return 0;
}