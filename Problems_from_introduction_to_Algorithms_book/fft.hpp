#ifndef FFT_HPP
#define FFT_HPP
#include <cmath>
#include <cstddef>
#include <vector>
#include <complex>

const double PI = acos(-1.0);

inline std::vector<std::complex<double>> fft(std::vector<std::complex<double>> a, bool inverse = false)
{
    size_t n = a.size();
    if (n == 1)
    {
        return a;
    }
    
    std::vector<std::complex<double>> a_0(n/2), a_1(n/2);
    for (size_t i = 0; i < n/2; i++)
    {
        a_0[i] = a[i*2];
        a_1[i] = a[i*2 + 1];
    }
    
    std::vector<std::complex<double>> y_0 = fft(a_0, inverse);
    std::vector<std::complex<double>> y_1 = fft(a_1, inverse);
    
    double ang = 2 * PI / n * (inverse ? -1 : 1);
    std::complex<double> w(1), wn(cos(ang), sin(ang));
    
    std::vector<std::complex<double>> y(n);
    for (size_t k = 0; k < n/2; k++)
    {
        y[k] = y_0[k] + w * y_1[k];
        y[k + n/2] = y_0[k] - w * y_1[k];
        
        if (inverse) {
            y[k] /= 2;
            y[k + n/2] /= 2;
        }
        
        w *= wn;
    }
    return y;
}
#endif //FFT_HPP