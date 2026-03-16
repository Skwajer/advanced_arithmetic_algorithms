#include "helpers.hpp"
#include <functional>
#include <cmath>
#include <iostream>
#include <stdexcept>

static double NewtonsMethod(
        double x0, 
        std::function<double(double x0)> compute_deravitive,
        std::function<double(double)> &f,
        double eps, 
        size_t multiplicity = 1,
        size_t max_iters = 1000)
{
    size_t iter = 0;

    do 
    {
        std::cout << "iter = " << iter << std::endl;
        auto x_next = x0 - multiplicity * (f(x0) / compute_deravitive(x0));
        if (std::abs(x_next - x0) < eps)
        {
            return x_next;
        }
        std::cout << "current approximation = " << x_next << std::endl;
        x0 = x_next;

        iter++;
    } while (iter < max_iters);

    throw std::runtime_error("Newtons method does not converge within max_iters");
}