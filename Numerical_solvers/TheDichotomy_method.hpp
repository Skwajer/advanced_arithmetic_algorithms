#include <cmath>
#include <stdexcept>
#include "helpers.hpp"
#include <iostream>


static double Dichotomy_method(double a, double b, std::function<double(double)> &f , double eps, size_t max_iters = 1000)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb >= 0)
    {
        throw std::runtime_error("func must be diff signs at the endpoints");
    }

    double fc;
    double c = 0;
    size_t iter = 0;
    do 
    {
        std::cout << "iter = " << iter << std::endl;
        c = (a + b) / 2;
        fc = f(c);
        if ((std::abs(fc) < eps) || ((a - b) / 2 < eps))
        {
            return c;
        }
        std::cout << "current approximation = " << c << std::endl;

        if (fa * fc < 0)
        {
            fb = fc;
            b = c;
        }else 
        {
            fa = fc;
            a = c;
        }

        iter++;
    } while (iter < max_iters);
    throw std::runtime_error("the dichotomy method does not converge within max_iters");
}