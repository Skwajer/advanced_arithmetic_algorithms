#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <boost/multiprecision/cpp_int.hpp>

using bigint = boost::multiprecision::cpp_int;

class fraction 
{
private:
    bigint _numerator;
    bigint _denominator;

    static bigint gcd(bigint enumerator, bigint const &denominator) 
    {
        bigint abs_denominator = bigint_abs(denominator);
        while (abs_denominator != 0) 
        {
            bigint temp = abs_denominator;
            abs_denominator = enumerator % abs_denominator;
            enumerator = std::move(temp);
        }
        return enumerator;
    }

    static bigint bigint_abs(bigint const &val) 
    {
        return (val < 0 ? -val : val);
    }

    void mutual_simplicity() 
    {
        if (_numerator == 0) 
        {
            _denominator = 1;
            return;
        }

        if (_numerator < 0) 
        {
            _numerator = -_numerator;
            _denominator = -_denominator;
        }

        bigint divider = gcd(_numerator, bigint_abs(_denominator));
        if (divider == 1) 
        {
            return;
        }

        _numerator /= divider;
        _denominator /= divider;
    }

public:
    fraction() : _numerator(0), _denominator(1) {}

    fraction(bigint const &numerator, bigint const &denominator) 
    {
        if (denominator == 0) {
            throw std::invalid_argument("Denominator can not be 0");
        }

        _numerator = bigint_abs(numerator);
        _denominator = bigint_abs(denominator);

        if (_numerator == 0) {
            _denominator = 1;
            return;
        }

        _denominator *= (((numerator > 0) ^ (denominator > 0)) ? -1 : 1);

        mutual_simplicity();
    }

    fraction(bigint const &other) : _numerator(bigint_abs(other)), _denominator(other < 0 ? -1 : 1) {}

    fraction(int other) : fraction(bigint(other)) {}

    fraction(fraction const &other) 
    {
        _numerator = bigint_abs(other._numerator);
        _denominator = other._denominator;
    }

    fraction &operator=(fraction const &other) 
    {
        if (this == &other) 
        {
            return *this;
        }
        _numerator = other._numerator;
        _denominator = other._denominator;
        return *this;
    }

    ~fraction() {}

    fraction abs() const 
    {
        return sign() == 1 ? (*this) : -(*this);
    }

    int sign() const 
    {
        if (_numerator == 0) 
        {
            return 0;
        }
        return _denominator < 0 ? -1 : 1;
    }

    fraction &operator+=(fraction const &other) 
    {
        _numerator *= other._denominator;
        bigint new_other_numerator = other._numerator * _denominator;
        _numerator += new_other_numerator;
        _denominator *= other._denominator;
        mutual_simplicity();
        return *this;
    }

    fraction operator+(fraction const &other) const 
    {
        fraction temp = *this;
        temp += other;
        return temp;
    }

    fraction &operator-=(fraction const &other) 
    {
        fraction temp = -other;
        *this += temp;
        return *this;
    }

    fraction operator-(fraction const &other) const 
    {
        fraction temp = *this;
        temp -= other;
        return temp;
    }

    fraction &operator*=(fraction const &other) 
    {
        _numerator *= other._numerator;
        _denominator *= other._denominator;
        mutual_simplicity();
        return *this;
    }

    fraction operator*(fraction const &other) const 
    {
        fraction temp = *this;
        temp *= other;
        return temp;
    }

    fraction &operator/=(fraction const &other) 
    {
        _numerator *= other._denominator;
        _denominator *= other._numerator;
        mutual_simplicity();
        return *this;
    }

    fraction operator/(fraction const &other) const 
    {
        fraction temp = *this;
        temp /= other;
        return temp;
    }

    fraction operator-() const 
    {
        fraction temp = *this;
        temp._denominator = -temp._denominator;
        return temp;
    }

    bool operator==(fraction const &other) const 
    {
        return (_numerator == other._numerator && _denominator == other._denominator);
    }

    bool operator!=(fraction const &other) const 
    {
        return !(*this == other);
    }

    bool operator>=(fraction const &other) const 
    {
        return !(*this < other);
    }

    bool operator>(fraction const &other) const 
    {
        bigint this_new_numerator = _numerator * other._denominator;
        bigint other_new_numerator = other._numerator * _denominator;
        return (this_new_numerator > other_new_numerator);
    }

    bool operator<=(fraction const &other) const 
    {
        return !(*this > other);
    }

    bool operator<(fraction const &other) const 
    {
        bigint this_new_numerator = _numerator * other._denominator;
        bigint other_new_numerator = other._numerator * _denominator;
        return (this_new_numerator < other_new_numerator);
    }

    fraction pow(size_t degree) const 
    {
        if (degree == 0) 
        {
            return 1;
        }
        fraction result = 1;
        fraction x = *this;
        size_t deg = degree;

        while (deg > 0) 
        {
            if (deg & 1) 
            {
                result *= x;
            }
            x *= x;
            deg >>= 1;
        }
        return result;
    }

    fraction arctg(fraction const &epsilon, long long &ops_count, int &terms_used) const 
    {
        fraction result = *this;
        ops_count = 0;
        terms_used = 0;
        
        if (result == 1) 
        {
            fraction frac_of1_5(1, 5);
            fraction frac_of1_239(1, 239);
            long long ops1 = 0, ops2 = 0;
            int terms1 = 0, terms2 = 0;
            fraction arctg1_5 = frac_of1_5.arctg(epsilon, ops1, terms1);
            fraction arctg1_239 = frac_of1_239.arctg(epsilon, ops2, terms2);
            result = fraction(4) * (arctg1_5 - arctg1_239);
            ops_count = ops1 + ops2 + 3;
            terms_used = terms1 + terms2;
            return result;
        }

        fraction term = result;
        fraction x = term;
        bigint n = 1;
        ops_count = 1;
        terms_used = 1;

        while (term.abs() > epsilon) 
        {
            term = fraction(-term * x * x * fraction(2 * n - 1) / fraction(2 * n + 1));
            result += term;
            ops_count += 5;
            terms_used++;
            ++n;
        }
        return result;
    }

    friend std::ostream &operator<<(std::ostream &stream, fraction const &num);
};

std::ostream &operator<<(std::ostream &stream, fraction const &num) 
{
    if (num == 0) {
        stream << 0;
        return stream;
    }
    if (num.sign() < 0) 
    {
        stream << "-";
    }
    
    fraction abs_num = num.abs();
    
    if (abs_num._denominator != 1) 
    {
        stream << abs_num._numerator << "/" << abs_num._denominator;
    } else {
        stream << abs_num._numerator;
    }
    return stream;
}

fraction calculate_pi(const fraction &epsilon, long long &total_ops, int &terms_used, double &time_ms) 
{
    fraction frac_1(1, 1);
    fraction one_fifth(1, 5);
    fraction one_239th(1, 239);
    
    long long ops1 = 0, ops2 = 0;
    int terms1 = 0, terms2 = 0;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    fraction arctan1_5 = one_fifth.arctg(epsilon, ops1, terms1);
    fraction arctan1_239 = one_239th.arctg(epsilon, ops2, terms2);
    
    fraction pi = frac_1 * 4;
    pi = pi * (frac_1 * 4 * arctan1_5 - arctan1_239);
    
    auto end = std::chrono::high_resolution_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    total_ops = ops1 + ops2 + 4;
    terms_used = terms1 + terms2;
    
    return pi;
}

fraction calculate_e_power_r(int r, const fraction &epsilon, long long &total_ops, int &terms_used, double &time_ms) 
{
    auto start = std::chrono::high_resolution_clock::now();
    
    fraction e(1);
    fraction term(1);
    fraction r_fraction(r);
    
    terms_used = 1;
    total_ops = 1;
    
    int k = 1;
    
    while (term.abs() > epsilon) 
    {
        term = term * r_fraction / k;
        e = e + term;
        total_ops += 3;
        terms_used++;
        k++;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    return e;
}

fraction precision_to_epsilon(int precision_digits) 
{
    bigint denominator = 1;
    for (int i = 0; i < precision_digits; ++i) 
    {
        denominator *= 10;
    }
    return fraction(1, denominator);
}


int main() 
{    
    int r;    
    std::cout << "Write natural number r (for e^r): ";
    std::cin >> r;
    
    if (r <= 0) {
        std::cout << "Error: r must be natural (positive) number\n";
        return 1;
    }
    
    std::vector<int> precisions = {5, 10, 15, 20, 30, 40, 50};
    
    std::cout << "\n============================================================\n";
    std::cout << "           Complexity parameters for e^" << r << " and π calculation\n";
    std::cout << "============================================================\n\n";
    
    for (int precision : precisions) 
    {
        fraction eps = precision_to_epsilon(precision);
        
        long long ops_e = 0, ops_pi = 0;
        int terms_e = 0, terms_pi = 0;
        double time_e_ms = 0, time_pi_ms = 0;
        
        fraction e_val = calculate_e_power_r(r, eps, ops_e, terms_e, time_e_ms);
        fraction pi_val = calculate_pi(eps, ops_pi, terms_pi, time_pi_ms);
        
        std::cout << "Accuracy: " << precision << " decimal Places\n";
        std::cout << "  e^" << r << " = " << e_val << "\n";
        std::cout << "  π = " << pi_val << "\n";
        std::cout << "  e^" << r << ": Time=" << std::fixed << std::setprecision(3) << time_e_ms 
                  << " ms, Terms=" << terms_e << ", Operations=" << ops_e << "\n";
        std::cout << "  π: Time=" << time_pi_ms 
                  << " ms, Terms=" << terms_pi << ", Operations=" << ops_pi << "\n";
        std::cout << "------------------------------------------------------------\n";
    }
    
    return 0;
}