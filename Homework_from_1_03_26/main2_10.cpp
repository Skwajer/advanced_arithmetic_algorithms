#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <optional>
#include <stdexcept>

using namespace std;

// ============================================================
// ОБРАЩЕНИЕ МАТРИЦЫ МЕТОДОМ ГАУССА (с присоединением единичной)
// ============================================================
std::optional<std::vector<std::vector<double>>> invertMatrix(
    std::vector<std::vector<double>> A)
{
    size_t n = A.size();
    
    // Создаём расширенную матрицу [A | I]
    vector<vector<double>> augmented(n, vector<double>(2 * n));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = A[i][j];              // левая часть: A
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0; // правая часть: I
        }
    }
    
    // Прямой ход метода Гаусса
    for (size_t col = 0; col < n; ++col) {
        // Поиск ведущего элемента
        size_t pivot = col;
        while (pivot < n && abs(augmented[pivot][col]) < 1e-15) {
            ++pivot;
        }
        
        if (pivot == n) {
            return std::nullopt;  // матрица вырождена
        }
        
        // Меняем строки местами
        if (pivot != col) {
            swap(augmented[col], augmented[pivot]);
        }
        
        // Нормализуем строку (делим на ведущий элемент)
        double pivot_val = augmented[col][col];
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented[col][j] /= pivot_val;
        }
        
        // Зануляем остальные строки в этом столбце
        for (size_t row = 0; row < n; ++row) {
            if (row != col && abs(augmented[row][col]) > 1e-15) {
                double factor = augmented[row][col];
                for (size_t j = 0; j < 2 * n; ++j) {
                    augmented[row][j] -= factor * augmented[col][j];
                }
            }
        }
    }
    
    // Извлекаем обратную матрицу из правой части
    vector<vector<double>> inv(n, vector<double>(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inv[i][j] = augmented[i][j + n];
        }
    }
    
    return inv;
}

// ============================================================
// МЕТОД НЬЮТОНА ДЛЯ СИСТЕМЫ 3x3
// ============================================================

static vector<double> NewtonSystem3D(
    const function<vector<double>(const vector<double>&)>& F,
    const function<vector<vector<double>>(const vector<double>&)>& J,
    const vector<double>& x0,
    double eps,
    size_t max_iters = 1000,
    bool verbose = true)
{
    vector<double> x = x0;
    vector<double> x_prev;
    size_t iter = 0;
    
    if (verbose) {
        cout << fixed << setprecision(10);
        cout << "\nНачальное приближение: ";
        for (size_t i = 0; i < x.size(); ++i) {
            cout << "x" << i+1 << " = " << x[i] << " ";
        }
        cout << "\n";
    }
    
    do
    {
        auto fx = F(x);
        auto Jx = J(x);
        
        auto Jinv_opt = invertMatrix(Jx);
        
        if (!Jinv_opt.has_value()) {
            throw runtime_error("Jacobian is singular at iteration " + to_string(iter));
        }
        
        auto Jinv = Jinv_opt.value();
        
        vector<double> delta(x.size(), 0.0);
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t j = 0; j < x.size(); ++j) {
                delta[i] -= Jinv[i][j] * fx[j];
            }
        }
        
        x_prev = x;
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += delta[i];
        }
        
        double norm_inf = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            norm_inf = max(norm_inf, abs(x[i] - x_prev[i]));
        }
        
        if (verbose) {
            cout << "iter " << setw(2) << iter << ": ";
            for (size_t i = 0; i < x.size(); ++i) {
                cout << "x" << i+1 << " = " << setw(12) << x[i] << " ";
            }
            cout << "||Δ||∞ = " << setw(12) << norm_inf << "\n";
        }
        
        if (norm_inf < eps) {
            if (verbose) {
                cout << "СОШЛОСЬ за " << iter + 1 << " итераций\n";
            }
            return x;
        }
        
        iter++;
    } while (iter < max_iters);
    
    throw runtime_error("Newton method did not converge within " + to_string(max_iters) + " iterations");
}



static void solve_part_a(double n, bool verbose = true)
{
    double eps = pow(10.0, -n);
    
    auto F = [](const vector<double>& x) -> vector<double> {
        double x1 = x[0], x2 = x[1], x3 = x[2];
        return {
            x1*x1*x1 + x1*x1*x2 - x1*x3 + 6.0,
            exp(x1) + exp(x2) - x3,
            x2*x2 - 2.0*x1*x3 - 4.0
        };
    };
    
    auto J = [](const vector<double>& x) -> vector<vector<double>> {
        double x1 = x[0], x2 = x[1], x3 = x[2];
        
        return {
            {3*x1*x1 + 2*x1*x2 - x3, x1*x1, -x1},
            {exp(x1), exp(x2), -1.0},
            {-2.0*x3, 2.0*x2, -2.0*x1}
        };
    };
    
    vector<double> x0 = {0.5, 0.5, 2.5};
    
    cout << "\n============================================================\n";
    cout << "ЧАСТЬ (a)\n";
    cout << "Точность: 10^-" << n << " = " << eps << "\n";
    cout << "============================================================\n";
    
    auto sol = NewtonSystem3D(F, J, x0, eps, 1000, verbose);
    
    cout << "\nРешение (a):\n";
    for (size_t i = 0; i < sol.size(); ++i) {
        cout << "x" << i+1 << " = " << sol[i] << "\n";
    }
}


static void solve_part_b(double n, bool verbose = true)
{
    double eps = pow(10.0, -n);
    
    auto F = [](const vector<double>& x) -> vector<double> {
        double x1 = x[0], x2 = x[1], x3 = x[2];
        return {
            6.0*x1 - 2.0*cos(x2 * x3) - 1.0,
            9.0*x2 + sqrt(x1*x1 + sin(x3)) + 1.06 + 0.9 - 0.0,
            60.0*x3 + 3.0*exp(-x1 * x2) + 10.0*M_PI - 3.0
        };
    };
    
    auto J = [](const vector<double>& x) -> vector<vector<double>> {
        double x1 = x[0], x2 = x[1], x3 = x[2];
        
        double denom = sqrt(x1*x1 + sin(x3));
        double dF2dx1 = (denom > 1e-12) ? x1 / denom : 0.0;
        double dF2dx3 = (denom > 1e-12) ? cos(x3) / (2.0 * denom) : 0.0;
        
        return {
            {6.0, 2.0 * sin(x2 * x3) * x3, 2.0 * sin(x2 * x3) * x2},
            {dF2dx1, 9.0, dF2dx3},
            {-3.0 * exp(-x1 * x2) * x2, -3.0 * exp(-x1 * x2) * x1, 60.0}
        };
    };
    
    vector<double> x0 = {0.2, -0.2, -0.5};
    
    cout << "\n============================================================\n";
    cout << "ЧАСТЬ (b)\n";
    cout << "Точность: 10^-" << n << " = " << eps << "\n";
    cout << "============================================================\n";
    
    auto sol = NewtonSystem3D(F, J, x0, eps, 1000, verbose);
    
    cout << "\nРешение (b):\n";
    for (size_t i = 0; i < sol.size(); ++i) {
        cout << "x" << i+1 << " = " << sol[i] << "\n";
    }
}


int main(int argc, char* argv[])
{
    double n = 8;
    
    if (argc > 1) {
        n = atof(argv[1]);
    }
    
    cout << "============================================================\n";
    cout << "МЕТОД НЬЮТОНА ДЛЯ СИСТЕМ 3x3\n";
    cout << "Критерий остановки: ||ΔX||∞ < 10^-n\n";
    cout << "n = " << n << "\n";
    cout << "============================================================\n";
    
    try {
        solve_part_a(n, true);
        solve_part_b(n, true);
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
        return 1;
    }
    
    cout << "\n============================================================\n";
    cout << "ВСЕ РАСЧЁТЫ ВЫПОЛНЕНЫ\n";
    cout << "============================================================\n";
    
    return 0;
}