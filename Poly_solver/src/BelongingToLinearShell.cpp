#include "../include/BelongingToLinearShell.h"

std::optional<std::vector<double>> gaussianElimination(
    std::vector<std::vector<double>> A, 
    std::vector<double> b)
{
    auto equations_count = A.size();
    auto variables_count = A[0].size();

    if (variables_count > equations_count)
    {
        return std::nullopt;
    }

    for (auto col = 0, row = 0; col < variables_count && row < equations_count; col++)
    {
        auto pivot = row;
        while (pivot < equations_count && std::abs(pivot) < Polynomial::EPSILON)
        {
            pivot++;
        }

        if (pivot == equations_count)
        {
            continue;
        }

        std::swap(A[row], A[pivot]);
        std::swap(b[row], b[pivot]);

        auto pivot_element = A[row][col];
        for(auto i = col; i < variables_count; i++)
        {
            A[row][i] /= pivot_element;
        }
        b[row] /= pivot_element;

        for (auto i = 0; i < equations_count; i++)
        {
            if (i != row && std::abs(A[i][col]) > Polynomial::EPSILON)
            {
                auto multiplier = A[i][col];
                for (auto j = col; j < variables_count; j ++)
                {
                    A[i][j] -= multiplier * A[row][j];
                }
                b[i] -= multiplier * b[row];
            }
        }
        row++;
    }

    for (size_t i = 0; i < equations_count; ++i)
    {
        bool all_zero = true;
        for (size_t j = 0; j < variables_count; ++j)
        {
            if (std::abs(A[i][j]) > Polynomial::EPSILON)
            {
                all_zero = false;
                break;
            }
        }
        
        if (all_zero && std::abs(b[i]) > Polynomial::EPSILON)
        {
            return std::nullopt;
        }
    }
    
    std::vector<double> x(variables_count, 0.0);
    for (size_t i = 0; i < variables_count; ++i)
    {
        for (size_t j = 0; j < equations_count; ++j)
        {
            if (std::abs(A[j][i] - 1.0) < Polynomial::EPSILON)
            {
                x[i] = b[j];
                break;
            }
        }
    }
    
    return x;
}

std::optional<std::vector<double>> isInLinearSpan(
    const Polynomial& f, 
    const std::vector<Polynomial>& basis)
{
    if (basis.empty()) return std::nullopt;
    
    size_t max_degree = f.getCoeffs().size();
    for (const auto& g : basis)
    {
        max_degree = std::max(max_degree, g.getCoeffs().size());
    }
    
    size_t n = basis.size();
    size_t m = max_degree;
    
    std::vector<std::vector<double>> A(m, std::vector<double>(n, 0.0));
    std::vector<double> b(m, 0.0);
    
    for (size_t j = 0; j < n; ++j)
    {
        const auto& coeffs = basis[j].getCoeffs();
        for (size_t i = 0; i < coeffs.size(); ++i)
        {
            A[i][j] = coeffs[i];
        }
    }
    
    const auto& f_coeffs = f.getCoeffs();
    for (size_t i = 0; i < f_coeffs.size(); ++i)
    {
        b[i] = f_coeffs[i];
    }
    
    auto solution = gaussianElimination(A, b);
    
    if (!solution.has_value())
    {
        return std::nullopt;
    }
    
    return solution;
}