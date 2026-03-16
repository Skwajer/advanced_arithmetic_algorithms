#include "../Numerical_solvers/Series_inverse.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

void print_coeffs(const std::vector<double>& coeffs, const std::string& name) {
    std::cout << name << " = [";
    for (size_t i = 0; i < coeffs.size(); ++i) {
        std::cout << std::fixed << std::setprecision(6) << coeffs[i];
        if (i != coeffs.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

int main() {
    try {
        // (truncated to M = 10)
        int M = 10;
        int n = 10;
        std::vector<double> a_exp(M + 1);
        double fact = 1.0;
        for (int k = 0; k <= M; ++k) {
            if (k > 0) fact *= k;
            a_exp[k] = 1.0 / fact;
        }
        
        std::vector<double> b_exp = invert_series(a_exp, n);
        print_coeffs(b_exp, "Обратный к экспоненте");

        std::vector<double> a_poly = {-1.0, -1.0, 1.0};
        a_poly.resize(n, 0.0);
        
        std::vector<double> b_poly = invert_series(a_poly, n);

        std::cout << "\nПроверка для экспоненты (f * g до x^M):" << std::endl;
        for (int k = 0; k < M; ++k) {
            double sum = 0.0;
            for (int i = 0; i <= k; ++i) {
                if (i < a_exp.size() && (k - i) < b_exp.size()) 
                {
                    sum += a_exp[i] * b_exp[k - i];
                }
            }
            std::cout << "Коэфф при x^" << k << ": " << sum << std::endl;
        }

        std::cout << "\nПроверка для x^2 - x - 1 (f * g до x^3):" << std::endl;
        for (int k = 0; k < 3; ++k) {
            double sum = 0.0;
            for (int i = 0; i <= M; ++i) {
                if (i < a_poly.size() && (k - i) < b_poly.size()) {
                    sum += a_poly[i] * b_poly[k - i];
                }
            }
            std::cout << "Коэфф при x^" << k << ": " << sum << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}