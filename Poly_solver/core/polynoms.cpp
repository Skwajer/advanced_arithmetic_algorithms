#include "polynoms.h"
#include <cmath>
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
    : coeffs(coefficients), expansion_center(0.0)
    {
        normalize();
    }

Polynomial::Polynomial(const std::vector<double>& coefficients, double init_center)
    : coeffs(coefficients), expansion_center(init_center)
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

    for (auto i = n - 2; i >= 1; i--)
    {
        quotient[i - 1] = quotient[i] * base_coeff + coeffs[i];
    }
    auto remainder = quotient[0] * base_coeff + coeffs[0];

    return {quotient, remainder};
}

std::vector<double> Polynomial::TaylorExpansion(double base_coeff) const
{
    std::vector<double> result{};
    std::vector<double> current_divisible = coeffs;

    while(current_divisible.size() > 0)
    {
        if (current_divisible.size() == 1)
        {
            result.push_back(current_divisible[0]);
            break;
        }

        Polynomial poly(current_divisible);
        auto [quotient, remainder] = poly.HornerDivide(base_coeff);
        result.push_back(remainder);
        current_divisible = quotient;
    }
    return result;
}


void Polynomial::get_decomposition_to_degrees(double a)
{
    if (coeffs.size() == 0)
    {
        std::cout << "polynomial is empty" << std::endl;
        return;
    }
    std::vector<double> expansion; 
    expansion = TaylorExpansion(a);
    coeffs = expansion;
    expansion_center = a;
    print_decomposition_to_degrees();
}

void Polynomial::print_decomposition_to_degrees() const
{
    expansion_center >= 0 ?
    std::cout << "decomposition to degrees (x - " << expansion_center << "):\n"
    :
    std::cout << "decomposition to degrees (x + " << -expansion_center << "):\n";
    
    for (size_t k = 0; k < coeffs.size(); k++) 
    {

        if (std::abs(coeffs[k]) < Polynomial::EPSILON) continue;
        
        std::cout << coeffs[k];
        if (k > 0) 
        {
            if (expansion_center > 0)
            {
                std::cout << "(x - " << expansion_center << ")";
            }
            else if (expansion_center < 0)
            {
                std::cout << "(x + " << -expansion_center << ")";
            }
            else 
            {
                std::cout << "x";
            }
            if (k > 1) 
            {
                std::cout << "^" << k;
            }
        }
        
        if (k < coeffs.size() - 1) std::cout << " + ";
    }
    std::cout << std::endl;
}

Polynomial Polynomial::get_shifted_representation(double new_center) const
{
    double delta = new_center - expansion_center;
     if (std::abs(delta) < Polynomial::EPSILON) 
        {
            return *this;
        }
    int n = coeffs.size() - 1;
    std::vector<double> new_coeffs(n + 1, 0.0); 

    for (int k = 0; k <= n; k++) 
    {         
        double binom = 1.0;
        double delta_pow = 1.0;
        
        delta_pow = pow(delta, k);
        
        for (int j = 0; j <= k; j++) 
        {
            double contribution = coeffs[k] * binom * delta_pow;
            new_coeffs[j] += contribution;
            
            
            if (j < k) 
            {
                binom = binom * (k - j) / (j + 1);
                delta_pow /= delta;
            }
        }
        std::cout << "\n";
    }

std::cout << "DEBUG: new_coeffs = ";
for (double x : new_coeffs) std::cout << x << " ";
std::cout << " | delta = " << delta << "\n";
    return Polynomial(new_coeffs, new_center);
}


