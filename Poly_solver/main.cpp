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

    std::vector<double> coeffs = {4, 3, 2, 1};  // x³ + 2x² + 3x + 4
    Polynomial poly(coeffs);
    
    poly.get_decomposition_by_degrees(2);
    return 0;
}