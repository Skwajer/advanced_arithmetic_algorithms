#include "../Numerical_solvers/Newton_system.hpp"

std::vector<double> solve_part_a(double A, double alpha2, double beta2, 
                                 const std::vector<double>& x0, double eps)
{
    auto F = [A, alpha2, beta2](const std::vector<double>& x) -> std::vector<double> {
        double xv = x[0], yv = x[1];
        return {
            tan(xv * yv + A) - xv * xv,
            xv * xv / alpha2 + yv * yv / beta2 - 1.0
        };
    };
    
    auto J = [A, alpha2, beta2](const std::vector<double>& x) -> std::vector<std::vector<double>> {
        double xv = x[0], yv = x[1];
        double sec2 = 1.0 / (cos(xv * yv + A) * cos(xv * yv + A));
        
        return {
            {yv * sec2 - 2 * xv, xv * sec2},
            {2 * xv / alpha2, 2 * yv / beta2}
        };
    };
    
    return NewtonSystem(F, J, x0, eps);
}

std::vector<double> solve_part_b(const std::vector<double>& x0, double eps)
{
    auto F = [](const std::vector<double>& x) -> std::vector<double> {
        double xv = x[0], yv = x[1];
        return {
            xv * xv + yv * yv - 2.0,
            exp(xv - 1.0) + yv * yv * yv - 2.0
        };
    };
    
    auto J = [](const std::vector<double>& x) -> std::vector<std::vector<double>> {
        double xv = x[0], yv = x[1];
        return {
            {2 * xv, 2 * yv},
            {exp(xv - 1.0), 3 * yv * yv}
        };
    };
    
    return NewtonSystem(F, J, x0, eps);
}

int main()
{
    double eps = 1e-8;
    
    std::cout << "============================================================\n";
    std::cout << "МЕТОД НЬЮТОНА ДЛЯ СИСТЕМ НЕЛИНЕЙНЫХ УРАВНЕНИЙ\n";
    std::cout << "Критерий остановки: ||F|| < eps И ||Δ|| < eps\n";
    std::cout << "eps = " << eps << "\n";
    std::cout << "============================================================\n";

    
    std::cout << "\n----------- ЧАСТЬ (a) -----------\n";
    
    {
        std::vector<double> x0 = {0.85, 0.5};
        auto sol = solve_part_a(0.2, 1.0/0.6, 1.0/2.0, x0, eps);
        std::cout << "\na-i) Решение: x = " << sol[0] << ", y = " << sol[1] << "\n";
    }
    
    {
        std::vector<double> x0 = {0.88, 0.43};
        auto sol = solve_part_a(0.4, 1.0/0.8, 1.0/2.0, x0, eps);
        std::cout << "\na-ii) Решение: x = " << sol[0] << ", y = " << sol[1] << "\n";
    }
    
    {
        std::vector<double> x0 = {1.0, 0.52};
        auto sol = solve_part_a(0.3, 1.0/0.2, 1.0/3.0, x0, eps);
        std::cout << "\na-iii) Решение: x = " << sol[0] << ", y = " << sol[1] << "\n";
    }
    
    {
        std::vector<double> x0 = {0.65, 0.6};
        auto sol = solve_part_a(0.0, 1.0/0.6, 1.0/2.0, x0, eps);
        std::cout << "\na-iv) Решение: x = " << sol[0] << ", y = " << sol[1] << "\n";
    }

    std::cout << "\n----------- ЧАСТЬ (b) -----------\n";
    
    {
        std::vector<double> x0 = {1.1, 0.9};
        auto sol = solve_part_b(x0, eps);
        std::cout << "\nb) Решение: x = " << sol[0] << ", y = " << sol[1] << "\n";
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "ВСЕ РАСЧЁТЫ ВЫПОЛНЕНЫ\n";
    std::cout << "============================================================\n";
    
    return 0;
}
