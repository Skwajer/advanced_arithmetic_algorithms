#include <iostream>
#include <string>
#include <cmath>
#include <numeric>

class Rational {
private:
    long long num;
    long long den;

    void normalize() {
        if (den == 0) {
            throw std::invalid_argument("Denominator cannot be zero");
        }
        
        if (den < 0) {
            num = -num;
            den = -den;
        }
        
        long long g = std::gcd(llabs(num), den);
        num /= g;
        den /= g;
    }

public:
    Rational() : num(0), den(1) {}
    
    Rational(long long n) : num(n), den(1) {}
    
    Rational(long long n, long long d) : num(n), den(d) {
        normalize();
    }
    
    Rational(const Rational& other) = default;
    Rational& operator=(const Rational& other) = default;
    
    Rational operator+(const Rational& other) const {
        return Rational(
            num * other.den + other.num * den,
            den * other.den
        );
    }
    
    Rational operator-(const Rational& other) const {
        return Rational(
            num * other.den - other.num * den,
            den * other.den
        );
    }
    
    Rational operator*(const Rational& other) const {
        return Rational(
            num * other.num,
            den * other.den
        );
    }
    
    Rational operator/(const Rational& other) const {
        if (other.num == 0) {
            throw std::invalid_argument("Division by zero");
        }
        return Rational(
            num * other.den,
            den * other.num
        );
    }
    
    // Унарный минус
    Rational operator-() const {
        return Rational(-num, den);
    }
    
    // Составные присваивания
    Rational& operator+=(const Rational& other) {
        *this = *this + other;
        return *this;
    }
    
    Rational& operator-=(const Rational& other) {
        *this = *this - other;
        return *this;
    }
    
    Rational& operator*=(const Rational& other) {
        *this = *this * other;
        return *this;
    }
    
    Rational& operator/=(const Rational& other) {
        *this = *this / other;
        return *this;
    }
    
    // Операторы сравнения
    bool operator==(const Rational& other) const {
        return num == other.num && den == other.den;
    }
    
    bool operator!=(const Rational& other) const {
        return !(*this == other);
    }
    
    bool operator<(const Rational& other) const {
        return num * other.den < other.num * den;
    }
    
    bool operator>(const Rational& other) const {
        return other < *this;
    }
    
    bool operator<=(const Rational& other) const {
        return !(*this > other);
    }
    
    bool operator>=(const Rational& other) const {
        return !(*this < other);
    }
    
    // Сравнение с целыми числами
    bool operator==(long long n) const {
        return *this == Rational(n);
    }
    
    // Преобразование в double
    double to_double() const {
        return static_cast<double>(num) / den;
    }
    
    // Получение числителя и знаменателя
    long long numerator() const { return num; }
    long long denominator() const { return den; }
    
    // Абсолютное значение
    Rational abs() const {
        return Rational(llabs(num), den);
    }
    
    // Ввод/вывод
    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        if (r.den == 1) {
            os << r.num;
        } else {
            os << r.num << "/" << r.den;
        }
        return os;
    }
    
    friend std::istream& operator>>(std::istream& is, Rational& r) {
        long long n, d = 1;
        char slash;
        
        is >> n;
        if (is.peek() == '/') {
            is >> slash >> d;
        }
        r = Rational(n, d);
        return is;
    }
};

// Дополнительные математические функции для Rational
Rational abs(const Rational& r) {
    return r.abs();
}

Rational pow(const Rational& base, int exponent) {
    if (exponent == 0) return Rational(1);
    if (exponent < 0) return Rational(1) / pow(base, -exponent);
    
    Rational result = base;
    for (int i = 1; i < exponent; ++i) {
        result = result * base;
    }
    return result;
}

// Специализация std::abs для Rational
namespace std {
    Rational abs(const Rational& r) {
        return r.abs();
    }
}