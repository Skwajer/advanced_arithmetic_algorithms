#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>
#include <cmath>

class Polynomial
{
    private:
    std::vector<double> coeffs; // from junior to senior degrees
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
    double evaluateAt(double x) const;
    static double limit_at_infinity(Polynomial const &f, Polynomial const &g, bool is_negate_infinity = false);
    static double limit_at_endpoint(Polynomial const &f, Polynomial const &g, double x);
    Polynomial compose(const Polynomial& g) const;
    static Polynomial iterated_compose(const Polynomial& f, int k, const Polynomial& inner);
    static double limit_T_at_point(
    const Polynomial& f1, int k, const Polynomial& s1,
    const Polynomial& f2, int l, const Polynomial& s2,
    double A);
    static double limit_T_at_infinity(
    const Polynomial& f1, int k, const Polynomial& s1,
    const Polynomial& f2, int l, const Polynomial& s2,
    bool positive_infinity = true);
    




    public:
    friend std::ostream& operator<<(std::ostream &os, const Polynomial &p);
};

Polynomial operator*(double scalar, const Polynomial& p);
