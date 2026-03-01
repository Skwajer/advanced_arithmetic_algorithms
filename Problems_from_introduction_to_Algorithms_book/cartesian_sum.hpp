#include "fft.hpp"
#include <complex>
#include <cstddef>
#include <vector>

inline std::vector<std::complex<double>> cartesian_sum(
    std::vector<int> a,
    std::vector<int> b,
    size_t n)
{
    size_t max_sum = n*20;
    auto even_degree  = 1;
    while (even_degree <= max_sum)
    {
        even_degree <<= 1;
    }

    std::vector<std::complex<double>> extended_a(even_degree);
    std::vector<std::complex<double>> extended_b(even_degree);
    for (size_t i = 0; i < a.size(); i++)
    {
        extended_a[a[i]] += 1.0;    
    }

    for (size_t i = 0; i < b.size(); i++)
    {
        extended_b[b[i]] += 1;
    }
    std::vector<std::complex<double>> Fa(even_degree), Fb(even_degree);
    Fa = fft(extended_a);
    Fb = fft(extended_b);

    std::vector<std::complex<double>> Fc(even_degree);
    for (size_t i = 0; i < even_degree; i++)
    {
        Fc[i] =  Fa[i] * Fb[i];
    }
    std::vector<std::complex<double>> c(even_degree);
    c = fft(Fc, true);

    return c;
}