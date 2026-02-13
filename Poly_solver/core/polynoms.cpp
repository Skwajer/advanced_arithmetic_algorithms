#include "polynoms.h"
#include <cstddef>
#include <vector>
#include <iostream>

void Polynomial::normalize()
{
    while( (coeffs.size() > 1) && (std::abs(coeffs.back()) < EPSILON) )
    {
        coeffs.pop_back();
    }
}

Polynomial::Polynomial(const std::vector<double>& coefficients)
    : coeffs(coefficients)
    {
        normalize();
    }

Polynomial::Polynomial(std::initializer_list<double> init) 
    : coeffs(init)
    {
        normalize();
    }

Polynomial& Polynomial::operator+=(const Polynomial& other) 
{
    size_t new_size = std::max(coeffs.size(), other.coeffs.size());
    coeffs.resize(new_size, 0.0);
    
    for (size_t i = 0; i < other.coeffs.size(); ++i) 
    {
        coeffs[i] += other.coeffs[i];
    }
    normalize();
    return *this;
}

Polynomial Polynomial::operator+(const Polynomial& other) const
{
    Polynomial res = *this;
    res += other;
    return res;
}

Polynomial& Polynomial::operator-=(const Polynomial& other) 
{
    size_t new_size = std::max(coeffs.size(), other.coeffs.size());
    coeffs.resize(new_size, 0.0);
    
    for (size_t i = 0; i < other.coeffs.size(); ++i) 
    {
        coeffs[i] -= other.coeffs[i];
    }
    normalize();
    return *this;
}

Polynomial Polynomial::operator-(const Polynomial& other) const
{
    Polynomial res = *this;
    res -= other;
    return res;
}

Polynomial& Polynomial::operator*=(double scalar)
{
    if (std::abs(scalar) < EPSILON)
    {
        coeffs = {0.0};
        return *this;
    }
    
    for (double& coeff : coeffs)
    {
        coeff *= scalar;
    }
    
    normalize();
    return *this;
}

Polynomial Polynomial::operator*(double scalar) const
{
    Polynomial result = *this;
    result *= scalar;
    return result;
}

Polynomial operator*(double scalar, const Polynomial& p)
{
    return p * scalar;
}

Polynomial Polynomial::operator-() const
{
    Polynomial result = *this;
    
    for (double& coeff : result.coeffs)
    {
        coeff = -coeff;
    }
    
    return result;
}

std::pair<std::vector<double>, double> Polynomial::HornerDivide(double base_coeff) const
{
    size_t n = coeffs.size();
    if (n == 1)
    {
        return {{}, coeffs[0]};
    }
    std::vector<double> quotient(n - 1);
    quotient[n - 2] = coeffs[n - 1];

    for (auto i = n - 2; i > 0 ;i--)
    {
        quotient[i - 1] = quotient[i] * base_coeff + coeffs[i];
    }
    auto remainder = quotient[0] * base_coeff + coeffs[0];

    return {quotient, remainder};
}

std::vector<double> Polynomial::TaylorExpansion(double base_coeff) const
{
    std::vector<double> result;
    std::vector<double> current_divisible = coeffs;

    while(current_divisible.size() > 0)
    {
        if (current_divisible.size() == 1)
        {
            result.push_back(current_divisible[0]);
        }

        Polynomial poly(current_divisible);
        auto [quotient, remainder] = poly.HornerDivide(base_coeff);
        result.push_back(remainder);
        current_divisible = quotient;
    }
    return result;
}


void Polynomial::get_decomposition_by_degrees(double a)
{
    std::vector<double> expansion; 
    expansion = TaylorExpansion(a);
    
    std::cout << "Разложение по степеням (x - " << a << "):\n";
    
    for (size_t k = 0; k < expansion.size(); k++) {

        if (std::abs(expansion[k]) < Polynomial::EPSILON) continue;
        
        std::cout << expansion[k];
        if (k > 0) 
        {
            std::cout << "(x - " << a << ")";
            if (k > 1) 
            {
                std::cout << "^" << k;
            }
        }
        
        if (k < expansion.size() - 1) std::cout << " + ";
    }
    std::cout << std::endl;
}

