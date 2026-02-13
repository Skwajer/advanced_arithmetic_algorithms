#include "../core/polynoms.h"
#include <optional>
std::optional<std::vector<double>> gaussianElimination(
    std::vector<std::vector<double>> A, 
    std::vector<double> b);

std::optional<std::vector<double>> isInLinearSpan(
    const Polynomial& f, 
    const std::vector<Polynomial>& basis);