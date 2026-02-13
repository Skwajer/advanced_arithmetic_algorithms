#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

class Polynomial
{
    private:
    std::vector<double> coeffs;
    double expansion_center;

    public:
    static constexpr double EPSILON = 1e-10;

    public:
    Polynomial() : coeffs({0.0}) 
    {}
    
    Polynomial(const std::vector<double>& coefficients);
    Polynomial(const std::vector<double>& coefficients, double init_center);
    Polynomial(std::initializer_list<double> init);
    
    Polynomial(const Polynomial& other) = default;

    public:
    void normalize();

    public:
    Polynomial& operator=(const Polynomial& other) = default;
    
    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);
    Polynomial& operator*=(double scalar);
    
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;
    Polynomial operator*(double scalar) const;
    
    Polynomial operator+() const;
    Polynomial operator-() const;

    public:
    double getCoeff(size_t i) const;
    const std::vector<double>& getCoeffs() const { return coeffs; }
    
    public:
    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;

    private:
    std::pair<std::vector<double>, double> HornerDivide(double base_coeff) const;
    std::vector<double> TaylorExpansion(double base_coeff) const;

    public:
    void print_decomposition_to_degrees() const;

    public:
    void get_decomposition_to_degrees(double a);
    [[nodiscard]] Polynomial get_shifted_representation(double new_center) const;

    
    public:
    friend std::ostream& operator<<(std::ostream & os, const Polynomial & p);
};

Polynomial operator*(double scalar, const Polynomial& p);
