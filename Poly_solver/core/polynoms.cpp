#include "polynoms.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <tuple>
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

Polynomial Polynomial::operator*(const Polynomial& other) const 
{
    std::vector<double> res(coeffs.size() + other.coeffs.size() - 1, 0.0);
    
    for (size_t i = 0; i < coeffs.size(); i++) 
    {
        for (size_t j = 0; j < other.coeffs.size(); j++) 
        {
            res[i + j] += coeffs[i] * other.coeffs[j];
        }
    }
    
    return Polynomial(res);
}


std::pair<std::vector<double>, double> Polynomial::HornerDivide(double base_coeff) const
{
    size_t n = coeffs.size();
    if (n == 0) 
    {
        return {{}, 0.0};
    }
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
    int n = coeffs.size();
    std::vector<double> new_coeffs = coeffs; 

    for (int i = 0; i < n; i++) 
    {
    for (int j = n - 1; j > i; j--) 
    {
        new_coeffs[j - 1] += delta * new_coeffs[j];
    }
}
    return Polynomial(new_coeffs, new_center);
}

double Polynomial::evaluateAt(double x) const
{
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; i--)
    {
        result = result * x + coeffs[i];
    }
    return result;
}

double Polynomial::limit_at_infinity(Polynomial const &f, Polynomial const &g, bool is_negate_infinity)
{
    if ((f.coeffs.empty() && g.coeffs.empty()))
    {
        throw std::invalid_argument("empty polynomials");
    }
    double sign = (f.coeffs.back() / g.coeffs.back()) > 0? 1 : -1;
    if (f.coeffs.size() - 1 > g.coeffs.size() - 1) 
    {
        if (is_negate_infinity)
        {
            return (f.coeffs.size() - g.coeffs.size()) % 2 == 0? std::numeric_limits<double>::infinity() * sign :
            std::numeric_limits<double>::infinity() * sign * -1;
        }
        else 
        {
            return std::numeric_limits<double>::infinity() * sign;
        }
    }

    else if (f.coeffs.size() - 1 < g.coeffs.size() - 1) {return 0.0;}
    else
    {
        return (f.coeffs.back() / g.coeffs.back());
    }
}

 double Polynomial::limit_at_endpoint(Polynomial const &f, Polynomial const &g, double x)
{
    Polynomial f_cpy = f;
    Polynomial g_cpy = g;
    Polynomial f_cp = f;
    Polynomial g_cp = g;

    if ((f.coeffs.empty() && g.coeffs.empty()))
    {
        throw std::invalid_argument("empty polynomials");
    }

    double g_abs_evaluated = std::abs(g.evaluateAt(x));
    double f_abs_evaluated = std::abs(f.evaluateAt(x));

    if (((g_abs_evaluated > Polynomial::EPSILON) && (f_abs_evaluated > Polynomial::EPSILON))
     || ((g_abs_evaluated > Polynomial::EPSILON) && (f_abs_evaluated < Polynomial::EPSILON)))
    {
        return (f.evaluateAt(x) / g.evaluateAt(x));
    }
    else if ((g_abs_evaluated < Polynomial::EPSILON) && (f_abs_evaluated > Polynomial::EPSILON)) 
    {
        return std::numeric_limits<double>::infinity() * ((f_abs_evaluated > Polynomial::EPSILON)? 1 : -1);
    }
    else //if (!f.evaluateAt(x) && !g.evaluateAt(x))
    {

        double f_remainder = 0.0, g_remainder = 0.0;
        std::vector<double> f_quotient, g_quotient;
        while (true)
        {
            std::tie(f_quotient, f_remainder) = f_cpy.HornerDivide(x);
            std::tie(g_quotient, g_remainder) = g_cpy.HornerDivide(x);
            if(std::abs(f_remainder) > Polynomial::EPSILON || std::abs(g_remainder) > Polynomial::EPSILON)
            {
                break;
            }

            f_cpy = Polynomial(f_quotient, 0.0);
            g_cpy = Polynomial(g_quotient, 0.0);
        }
        if (((std::abs(f_cpy.evaluateAt(x)) > Polynomial::EPSILON) && (std::abs(g_cpy.evaluateAt(x)) > Polynomial::EPSILON)) || 
        ((std::abs(f_cpy.evaluateAt(x)) < Polynomial::EPSILON) && std::abs(g_cpy.evaluateAt(x)) > Polynomial::EPSILON))
        {
            return (f_cpy.evaluateAt(x) / g_cpy.evaluateAt(x));
        }
        else
        {
            return std::numeric_limits<double>::infinity() * (f_cpy.evaluateAt(x) > 0? 1 : -1);
        }
    }
}

Polynomial Polynomial::compose(const Polynomial& g) const
{
    
    Polynomial result({0});
    for (int i = coeffs.size() - 1; i >= 0; i--) 
    {

        Polynomial term = result * g;
        std::vector<double> new_coeffs = term.coeffs;
        
        if (new_coeffs.size() < 1) new_coeffs.resize(1, 0.0);
        new_coeffs[0] += coeffs[i];
        
        result = Polynomial(new_coeffs);
    }
    return result;
}

Polynomial Polynomial::iterated_compose(const Polynomial& f, int k, const Polynomial& inner) 
{
    Polynomial result = inner;
    for (int i = 0; i < k; i++) {
        result = f.compose(result);
    }
    return result;
}

double Polynomial::limit_T_at_point(
    const Polynomial& f1, int k, const Polynomial& s1,
    const Polynomial& f2, int l, const Polynomial& s2,
    double A)
{
    Polynomial num = iterated_compose(f1, k, s1);
    Polynomial den = iterated_compose(f2, l, s2);
    return limit_at_endpoint(num, den, A);
}

 double Polynomial::limit_T_at_infinity(
    const Polynomial& f1, int k, const Polynomial& s1,
    const Polynomial& f2, int l, const Polynomial& s2,
    bool positive_infinity)
{
    Polynomial num = Polynomial::iterated_compose(f1, k, s1);
    Polynomial den = Polynomial::iterated_compose(f2, l, s2);
    return Polynomial::limit_at_infinity(num, den, positive_infinity);
}


