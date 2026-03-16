#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <functional>

static double computeDeravative(double x0, std::function<double(double)> &f, double h = 1e-8)
{
    return  (f(x0 + h) - f(x0 - h)) / (2*h);
}

#endif
