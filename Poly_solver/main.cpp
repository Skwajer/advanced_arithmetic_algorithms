#include "include/BelongingToLinearShell.h"

int main()
{
    // Polynomial g1 = {1.0, 1.0, -2.0}; // 1 + x
    // Polynomial g2 = {1.0, -1.0, 1.0}; // 1 - x
    // Polynomial g3 = {1.0, -2.0, 2.0};

    // Polynomial f = {2.0, 0.0, 2};

    // auto result = isInLinearSpan(f, {g1, g2, g3});

    // if (result.has_value())
    // {
    // std::cout << "f belongs to the linear shell\n";
    // std::cout << "decomposition coefficients:\n";
    // for (size_t i = 0; i < result->size(); ++i)
    // {
    // std::cout << "A" << i+1 << " = " << (*result)[i] << "\n";
    // }
    // }

    std::vector<double> coeffs{4, 12, 7, 0, 1};  // x^4 + 0x^3 + 7x^2 + 12x^1 + 4
    
    Polynomial poly(coeffs);
    poly.print_decomposition_to_degrees();
    poly.get_decomposition_to_degrees(-3);
    Polynomial shifted_poly = poly.get_shifted_representation(5);
    shifted_poly.print_decomposition_to_degrees();
    return 0;
}