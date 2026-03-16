#include <functional>
#include <optional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

static std::optional<std::vector<double>> gaussianElimination(
    std::vector<std::vector<double>> A, 
    std::vector<double> b,
    double eps = 1e-10)
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
        while (pivot < equations_count && std::abs(pivot) < eps)
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
            if (i != row && std::abs(A[i][col]) > eps)
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
            if (std::abs(A[i][j]) > eps)
            {
                all_zero = false;
                break;
            }
        }
        
        if (all_zero && std::abs(b[i]) > eps)
        {
            return std::nullopt;
        }
    }
    
    std::vector<double> x(variables_count, 0.0);
    for (size_t i = 0; i < variables_count; ++i)
    {
        for (size_t j = 0; j < equations_count; ++j)
        {
            if (std::abs(A[j][i] - 1.0) < eps)
            {
                x[i] = b[j];
                break;
            }
        }
    }
    
    return x;
}

static std::vector<std::vector<double>> inverseMatrix(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    
    if (n == 0 || A[0].size() != n) {
        throw std::invalid_argument("Матрица должна быть квадратной");
    }
    
    std::vector<std::vector<double>> augmented(n, std::vector<double>(2 * n));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(augmented[k][i]) > fabs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }
        
        if (fabs(augmented[maxRow][i]) < 1e-10) {
            throw std::invalid_argument("Матрица вырожденная, обратной не существует");
        }
        
        swap(augmented[i], augmented[maxRow]);
        
        double pivot = augmented[i][i];
        for (int j = 0; j < 2 * n; j++) {
            augmented[i][j] /= pivot;
        }
        
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }
    
    std::vector<std::vector<double>> inverse(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmented[i][j + n];
        }
    }
    
    return inverse;
}


struct SystemResult {
    std::vector<double> solution;
    size_t iterations;
    bool converged;
};

static std::vector<double> NewtonSystem(
    const std::function<std::vector<double>(const std::vector<double>&)>& F,
    const std::function<std::vector<std::vector<double>>(const std::vector<double>&)>& J, 
    const std::vector<double>& x0,
    double eps = 1e-10,
    size_t max_iters = 100000,
    bool verbose = true 
)
{
    std::vector<double> x = x0;
    std::vector<double> x_prev;
    size_t iter = 0;
    
    if (verbose) {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << "\nНачальное приближение: ";
        for (size_t i = 0; i < x.size(); ++i) {
            std::cout << "x" << i << " = " << x[i] << " ";
        }
        std::cout << "\n";
    }
    
    do
    {
        auto fx = F(x);
        auto Jx = J(x);
        
        for (auto& val : fx) {
            val = -val;
        }
        
        auto delta_opt = gaussianElimination(Jx, fx);
        
        if (!delta_opt.has_value()) {
            throw std::runtime_error("Jacobian is singular at iteration " + std::to_string(iter));
        }
        
        auto delta = delta_opt.value();
        
        x_prev = x;
        for (size_t i = 0; i < x.size(); ++i) 
        {
            x[i] += delta[i];
        }
        
        double norm_delta = 0.0;
        for (size_t i = 0; i < delta.size(); ++i) 
        {
            norm_delta += delta[i] * delta[i];
        }
        norm_delta = std::sqrt(norm_delta);
        
        double norm_F = 0.0;
        for (size_t i = 0; i < fx.size(); ++i) 
        {
            norm_F += fx[i] * fx[i];
        }
        norm_F = std::sqrt(norm_F);
        
        if (verbose) {
            std::cout << "iter " << std::setw(3) << iter << ": ";
            for (size_t i = 0; i < x.size(); ++i) {
                std::cout << "x" << i << " = " << std::setw(12) << x[i] << " ";
            }
            std::cout << "|F| = " << std::setw(12) << norm_F
                      << " |Δ| = " << std::setw(12) << norm_delta << "\n";
        }
        
        if (norm_F < eps && norm_delta < eps) {
            if (verbose) {
                std::cout << "СОШЛОСЬ за " << iter + 1 << " итераций\n";
            }
            return x;
        }
        
        iter++;
    } while (iter < max_iters);
    
    throw std::runtime_error("Newton method did not converge within " + std::to_string(max_iters) + " iterations");
}

static double normInf(const std::vector<double>& x, const std::vector<double>& y) {
    double max_diff = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double diff = fabs(x[i] - y[i]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    return max_diff;
}

// Метод Ньютона с обратной матрицей Якоби
static std::vector<double> NewtonSystem(
    const std::function<std::vector<double>(const std::vector<double>&)>& F,
    const std::function<std::vector<std::vector<double>>(const std::vector<double>&)>& J,
    const std::vector<double>& x0,
    int n_accuracy = 10,  // n для критерия остановки: 10^(-n)
    size_t max_iters = 1000,
    bool verbose = true
)
{
    std::vector<double> x = x0;
    std::vector<double> x_prev = x0;
    size_t dim = x0.size();
    
    double tolerance = pow(10.0, -n_accuracy);
    
    if (verbose) {
        std::cout << "Метод Ньютона для системы нелинейных уравнений\n";
        std::cout << "Критерий остановки: ||X^k - X^{k-1}||_∞ < 10^(-" << n_accuracy << ") = " << tolerance << "\n";
        std::cout << "Начальное приближение: ";
        for (size_t i = 0; i < dim; i++) {
            std::cout << "x" << i+1 << " = " << x[i] << " ";
        }
        std::cout << "\n\n";
        std::cout << std::setw(5) << "iter" << std::setw(15) << "||Δx||_∞" << std::setw(15) << "||F(x)||";
        for (size_t i = 0; i < dim; i++) {
            std::cout << std::setw(15) << "x" << i+1;
        }
        std::cout << "\n";
    }
    
    for (size_t iter = 0; iter < max_iters; iter++) {
        std::vector<double> Fx = F(x);
        
        double F_norm = 0.0;
        for (size_t i = 0; i < dim; i++) {
            F_norm += Fx[i] * Fx[i];
        }
        F_norm = sqrt(F_norm);
        
        std::vector<std::vector<double>> Jx = J(x);
        std::vector<std::vector<double>> J_inv = inverseMatrix(Jx);
        
        std::vector<double> dx(dim, 0.0);
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = 0; j < dim; j++) {
                dx[i] -= J_inv[i][j] * Fx[j];
            }
        }
        
        x_prev = x;
        for (size_t i = 0; i < dim; i++) {
            x[i] += dx[i];
        }
        
        double diff_norm = normInf(x, x_prev);
        
        if (verbose) {
            std::cout << std::setw(5) << iter 
                 << std::setw(15) << std::scientific << std::setprecision(5) << diff_norm
                 << std::setw(15) << std::scientific << std::setprecision(5) << F_norm;
            for (size_t i = 0; i < dim; i++) {
                std::cout << std::setw(15) << std::fixed << std::setprecision(8) << x[i];
            }
            std::cout << "\n";
        }
        
        if (diff_norm < tolerance) 
        {
            if (verbose) {
                std::cout << "\nРешение найдено за " << iter + 1 << " итераций\n";
            }
            return x;
        }
    }
    
    throw std::runtime_error("Newtons method does not converge within max_iters");
}

