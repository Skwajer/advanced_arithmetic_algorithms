#include <cstdio>
#include <stdexcept>
#include "helpers.hpp"

static double computeLambda(double x0, std::function<double(double)> &f)
{
    return (-1 / computeDeravative(x0, f));

}

static double simpleIteration(double x0, std::function<double(double x)> f, double eps = 1e-6, size_t max_iters = 100000)
{
    double x1;
    size_t iter = 0;
    auto lambda = computeLambda(x0, f);
    printf("derivative = %f\n", computeDeravative(x0, f));

    printf("lambda = %f\n", lambda);
    do
    {
        x1 = x0 + lambda * f(x0);
        if (iter % 2 ==0) printf("%llu iter, x1 = %f\n", iter, x1);
        if (std::abs(x1 - x0) < eps)
        {
            return x1;
        }
        x0 = x1;
        iter++;

    } while(iter < max_iters);

    throw std::runtime_error("the simple iteration method does not converge within max_iters");
}