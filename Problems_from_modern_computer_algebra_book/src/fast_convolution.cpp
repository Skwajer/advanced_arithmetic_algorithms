#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

using coeffs = std::vector<std::complex<double>>;

std::pair<std::vector<std::vector<std::complex<double>>>,
          std::vector<std::vector<std::complex<double>>>>
precompute_omega_powers(size_t k_max)
{
    std::vector<std::vector<std::complex<double>>> omega_pow(k_max + 1);
    std::vector<std::vector<std::complex<double>>> omega_inv(k_max + 1);
    
    for (size_t k = 1; k <= k_max; ++k)
    {
        size_t n = 1 << k;
        size_t m = n / 2;
        
        double angle = 2.0 * M_PI / n;
        std::complex<double> omega_n(cos(angle), sin(angle));
        std::complex<double> omega_n_inv = conj(omega_n);
        
        omega_pow[k].resize(m);
        omega_inv[k].resize(m);
        
        std::complex<double> pow = 1.0;
        std::complex<double> inv_pow = 1.0;
        
        for (size_t j = 0; j < m; ++j)
        {
            omega_pow[k][j] = pow;
            omega_inv[k][j] = inv_pow;
            pow *= omega_n;
            inv_pow *= omega_n_inv;
        }
    }
    
    return {omega_pow, omega_inv};
}

coeffs convolution(
    coeffs const &f,
    coeffs const &g,
    size_t k,
    size_t d,
    std::vector<std::vector<std::complex<double>>> const& omega_pow,
    std::vector<std::vector<std::complex<double>>> const& omega_inv)
{
    size_t n = 1 << k;
    
    if (k <= d)
    {
        coeffs prod(f.size() + g.size() - 1, 0.0);
        for (size_t i = 0; i < f.size(); ++i)
            for (size_t j = 0; j < g.size(); ++j)
                prod[i + j] += f[i] * g[j];
        
        coeffs result(n, 0.0);
        for (size_t i = 0; i < prod.size(); ++i)
            result[i % n] += prod[i];
        
        return result;
    }
    
    size_t m = n / 2;
    
    coeffs f0(m, 0.0), f1(m, 0.0);
    coeffs g0(m, 0.0), g1(m, 0.0);
    
    for (size_t i = 0; i < m; ++i)
    {
        f0[i] = f[i] + f[i + m];
        f1[i] = f[i] - f[i + m];
        g0[i] = g[i] + g[i + m];
        g1[i] = g[i] - g[i + m];
    }
        
    auto h0 = convolution(f0, g0, k - 1, d, omega_pow, omega_inv);
    
    const auto &pow_k = omega_pow[k];
    coeffs f1_twisted(m, 0.0), g1_twisted(m, 0.0);
    for (size_t j = 0; j < m; ++j)
    {
        f1_twisted[j] = f1[j] * pow_k[j];
        g1_twisted[j] = g1[j] * pow_k[j];
    }
    
    auto H = convolution(f1_twisted, g1_twisted, k - 1, d, omega_pow, omega_inv);
    
    const auto &inv_k = omega_inv[k];
    coeffs h1(m, 0.0);
    for (size_t j = 0; j < m; ++j)
    {
        h1[j] = H[j] * inv_k[j];
    }
    
    coeffs result(n, 0.0);
    for (size_t i = 0; i < m; ++i)
    {
        result[i] = (h0[i] + h1[i]) * 0.5;
        result[i + m] = (h0[i] - h1[i]) * 0.5;
    }
    
    return result;
}

void print_coeffs(const coeffs& c)
{
    for (auto& x : c)
        std::cout << x << " ";
    std::cout << "\n";
}

void test_case(const coeffs& result, const coeffs& expected, const std::string& test_name)
{
    std::cout << test_name << ": ";
    bool ok = true;
    for (size_t i = 0; i < result.size(); ++i)
    {
        double diff_real = std::abs(result[i].real() - expected[i].real());
        double diff_imag = std::abs(result[i].imag() - expected[i].imag());
        if (diff_real > 1e-10 || diff_imag > 1e-10)
        {
            ok = false;
            break;
        }
    }
    if (ok)
        std::cout << "✓ ПРОЙДЕН\n";
    else
        std::cout << "✗ НЕ ПРОЙДЕН\n";
}

int main()
{
    size_t k = 3;
    size_t d = 1;
    
    std::complex<double> omega = std::exp(2.0 * M_PI * std::complex<double>(0, 1) / 8.0);
    
    auto [omega_pow, omega_inv] = precompute_omega_powers(k);
    
    // ------------------------------------------------------------
    // ТЕСТ 1: f = 1 + x, g = 1 + x
    // (1+x)^2 = 1 + 2x + x^2
    // ------------------------------------------------------------
    {
        coeffs f = {1, 1, 0, 0, 0, 0, 0, 0};
        coeffs g = {1, 1, 0, 0, 0, 0, 0, 0};
        coeffs expected = {1, 2, 1, 0, 0, 0, 0, 0};
        print_coeffs(expected);
        auto result = convolution(f, g, k, d, omega_pow, omega_inv);
        print_coeffs(result);
        test_case(result, expected, "Тест 1 (1+x) * (1+x)");
    }
    
    // ------------------------------------------------------------
    // ТЕСТ 2: f = 1 + x + x^2, g = 1 + x
    // (1+x+x^2)(1+x) = 1 + 2x + 2x^2 + x^3
    // ------------------------------------------------------------
    {
        coeffs f = {1, 1, 1, 0, 0, 0, 0, 0};
        coeffs g = {1, 1, 0, 0, 0, 0, 0, 0};
        coeffs expected = {1, 2, 2, 1, 0, 0, 0, 0};
        print_coeffs(expected);
        auto result = convolution(f, g, k, d, omega_pow, omega_inv);
        print_coeffs(result);
        test_case(result, expected, "Тест 2 (1+x+x^2) * (1+x)");
    }
    
    // ------------------------------------------------------------
    // ТЕСТ 3: Единичные многочлены: f = 1, g = 1
    // Результат: 1
    // ------------------------------------------------------------
    {
        coeffs f = {1, 0, 0, 0, 0, 0, 0, 0};
        coeffs g = {1, 0, 0, 0, 0, 0, 0, 0};
        coeffs expected = {1, 0, 0, 0, 0, 0, 0, 0};
        print_coeffs(expected);
        auto result = convolution(f, g, k, d, omega_pow, omega_inv);
        print_coeffs(result);
        test_case(result, expected, "Тест 3 (1) * (1)");
    }
    
    // ------------------------------------------------------------
    // ТЕСТ 4: f = x, g = x^7 (проверка циклического переноса)
    // x * x^7 = x^8 ≡ 1 (mod x^8 - 1)
    // ------------------------------------------------------------
    {
        coeffs f = {0, 1, 0, 0, 0, 0, 0, 0};        // x
        coeffs g = {0, 0, 0, 0, 0, 0, 0, 1};        // x^7
        coeffs expected = {1, 0, 0, 0, 0, 0, 0, 0}; // 1 (после циклического переноса)
        print_coeffs(expected);
        auto result = convolution(f, g, k, d, omega_pow, omega_inv);
        print_coeffs(result);
        test_case(result, expected, "Тест 4 (x) * (x^7)");
    }
    
    // ------------------------------------------------------------
    // ТЕСТ 5: f = 1 + x^4, g = 1 + x^4
    // (1+x^4)^2 = 1 + 2x^4 + x^8 ≡ 2 + 2x^4 (mod x^8-1)
    // ------------------------------------------------------------
    {
        coeffs f = {1, 0, 0, 0, 1, 0, 0, 0};
        coeffs g = {1, 0, 0, 0, 1, 0, 0, 0};
        coeffs expected = {2, 0, 0, 0, 2, 0, 0, 0}; // x^8 ≡ 1, поэтому 1 + 1 = 2 в x^0
        print_coeffs(expected);
        
        auto result = convolution(f, g, k, d, omega_pow, omega_inv);
        print_coeffs(result);
        test_case(result, expected, "Тест 5 (1+x^4) * (1+x^4)");
    }

    // ------------------------------------------------------------
// ТЕСТ 6: f = 1 + x + x^2 + x^3, g = 1 + x 
// (1 + x + x^2 + x^3)(1 + x) = 1 + 2x + 2x^2 + 2x^3 + x^4
// ------------------------------------------------------------
{
    coeffs f = {1, 1, 1, 1, 0, 0, 0, 0};
    coeffs g = {1, 1, 0, 0, 0, 0, 0, 0};
    coeffs expected = {1, 2, 2, 2, 1, 0, 0, 0};
    print_coeffs(expected);
    auto result = convolution(f, g, k, d, omega_pow, omega_inv);
    print_coeffs(result);
    test_case(result, expected, "Тест 6 (1+x+x^2+x^3) * (1+x)");
}

// ------------------------------------------------------------
// ТЕСТ 7: f = 1 - x, g = 1 + x 
// (1-x)(1+x) = 1 - x^2
// ------------------------------------------------------------
{
    coeffs f = {1, -1, 0, 0, 0, 0, 0, 0};
    coeffs g = {1, 1, 0, 0, 0, 0, 0, 0};
    coeffs expected = {1, 0, -1, 0, 0, 0, 0, 0};
    print_coeffs(expected);
    auto result = convolution(f, g, k, d, omega_pow, omega_inv);
    print_coeffs(result);
    test_case(result, expected, "Тест 7 (1-x) * (1+x)");
}

// ------------------------------------------------------------
// ТЕСТ 8: f = x^2, g = x^6 
// x^2 * x^6 = x^8 ≡ 1 (mod x^8 - 1)
// ------------------------------------------------------------
{
    coeffs f = {0, 0, 1, 0, 0, 0, 0, 0};        // x^2
    coeffs g = {0, 0, 0, 0, 0, 0, 1, 0};        // x^6
    coeffs expected = {1, 0, 0, 0, 0, 0, 0, 0}; // 1 (после циклического переноса)
    print_coeffs(expected);
    auto result = convolution(f, g, k, d, omega_pow, omega_inv);
    print_coeffs(result);
    test_case(result, expected, "Тест 8 (x^2) * (x^6)");
}

// ------------------------------------------------------------
// ТЕСТ 9: f = 1 + x^2 + x^4 + x^6, g = 1 + x^2
// (1+x^2+x^4+x^6)(1+x^2) = 1 + 2x^2 + 2x^4 + 2x^6 + x^8
// x^8 ≡ 1, поэтому: (1+1) + 2x^2 + 2x^4 + 2x^6 = 2 + 2x^2 + 2x^4 + 2x^6
// ------------------------------------------------------------
{
    coeffs f = {1, 0, 1, 0, 1, 0, 1, 0};
    coeffs g = {1, 0, 1, 0, 0, 0, 0, 0};
    coeffs expected = {2, 0, 2, 0, 2, 0, 2, 0};
    print_coeffs(expected);
    auto result = convolution(f, g, k, d, omega_pow, omega_inv);
    print_coeffs(result);
    test_case(result, expected, "Тест 9 (1+x^2+x^4+x^6) * (1+x^2)");
}

// ------------------------------------------------------------
// ТЕСТ 10: f = 1 + x^3, g = 1 + x^5
// (1+x^3)(1+x^5) = 1 + x^3 + x^5 + x^8 ≡ 1 + x^3 + x^5 + 1 = 2 + x^3 + x^5
// ------------------------------------------------------------
{
    coeffs f = {1, 0, 0, 1, 0, 0, 0, 0};
    coeffs g = {1, 0, 0, 0, 0, 1, 0, 0};
    coeffs expected = {2, 0, 0, 1, 0, 1, 0, 0};
    print_coeffs(expected);
    auto result = convolution(f, g, k, d, omega_pow, omega_inv);
    print_coeffs(result);
    test_case(result, expected, "Тест 10 (1+x^3) * (1+x^5)");
}
    
    return 0;
}