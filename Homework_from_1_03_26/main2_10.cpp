#include <iostream>
#include <iomanip>
#include "../Numerical_solvers/Newton_system.hpp"

using namespace std;
void solveSystemA() 
{
    cout << "========================================\n";
    cout << "Система a:\n";
    cout << "x1^3 + x1^2*x2 - x1*x3 + 6 = 0\n";
    cout << "e^x1 + e^x2 - x3 = 0\n";
    cout << "x2^2 - 2*x1*x3 = 4\n\n";
    
    auto F = [](const vector<double>& x) {
        vector<double> f(3);
        f[0] = x[0]*x[0]*x[0] + x[0]*x[0]*x[1] - x[0]*x[2] + 6.0;
        f[1] = exp(x[0]) + exp(x[1]) - x[2];
        f[2] = x[1]*x[1] - 2.0*x[0]*x[2] - 4.0;
        return f;
    };
    
    auto J = [](const vector<double>& x) {
        vector<vector<double>> Jx(3, vector<double>(3));
        
        Jx[0][0] = 3.0*x[0]*x[0] + 2.0*x[0]*x[1] - x[2];
        Jx[0][1] = x[0]*x[0];
        Jx[0][2] = -x[0];
        
        Jx[1][0] = exp(x[0]);
        Jx[1][1] = exp(x[1]);
        Jx[1][2] = -1.0;
        
        Jx[2][0] = -2.0*x[2];
        Jx[2][1] = 2.0*x[1];
        Jx[2][2] = -2.0*x[0];
        
        return Jx;
    };
    
    vector<vector<double>> initial_guesses = {
        {-1.0, 1.0, 2.0},
        {0.5, 0.5, 3.0},
        {-2.0, 2.0, 1.0}
    };
    
    for (size_t i = 0; i < initial_guesses.size(); i++) {
        try {
            cout << "\n--- Попытка " << i+1 << " ---\n";
            vector<double> solution = NewtonSystem(F, J, initial_guesses[i], 10, 100, true);
            
            cout << "\nНайденное решение:\n";
            for (size_t j = 0; j < solution.size(); j++) {
                cout << "x" << j+1 << " = " << setprecision(10) << solution[j] << "\n";
            }
            
            vector<double> check = F(solution);
            cout << "\nНевязки:\n";
            for (size_t j = 0; j < check.size(); j++) {
                cout << "f" << j+1 << " = " << scientific << check[j] << "\n";
            }
            cout << "----------------------------------------\n";
            
        } catch (const exception& e) {
            cerr << "Ошибка: " << e.what() << "\n";
        }
    }
}

void solveSystemB() {
    cout << "========================================\n";
    cout << "Система b:\n";
    cout << "6x1 - 2cos(x2*x3) - 1 = 0\n";
    cout << "9x2 + sqrt(x1^2 + sin(x3)) + 1.06 + 0.9 = 0\n";
    cout << "60x3 + 3e^(-x1*x2) + 10π - 3 = 0\n\n";
    
    auto F = [](const vector<double>& x) {
        vector<double> f(3);
        f[0] = 6.0*x[0] - 2.0*cos(x[1]*x[2]) - 1.0;
        
        double under_sqrt = x[0]*x[0] + sin(x[2]);
        if (under_sqrt < 0) under_sqrt = 0;
        f[1] = 9.0*x[1] + sqrt(under_sqrt) + 1.06 + 0.9;
        
        f[2] = 60.0*x[2] + 3.0*exp(-x[0]*x[1]) + 10.0*M_PI - 3.0;
        return f;
    };
    
    auto J = [](const vector<double>& x) {
        vector<vector<double>> Jx(3, vector<double>(3));
        
        Jx[0][0] = 6.0;
        Jx[0][1] = 2.0 * x[2] * sin(x[1]*x[2]);
        Jx[0][2] = 2.0 * x[1] * sin(x[1]*x[2]);
        
        double under_sqrt = x[0]*x[0] + sin(x[2]);
        if (under_sqrt < 0) under_sqrt = -under_sqrt;
        
        double denom = sqrt(under_sqrt);
        if (denom < 1e-15) denom = 1e-15;
        
        Jx[1][0] = x[0] / denom;
        Jx[1][1] = 9.0;
        Jx[1][2] = cos(x[2]) / (2.0 * denom);
        
        Jx[2][0] = -3.0 * x[1] * exp(-x[0]*x[1]);
        Jx[2][1] = -3.0 * x[0] * exp(-x[0]*x[1]);
        Jx[2][2] = 60.0;
        
        return Jx;
    };
    
    vector<vector<double>> initial_guesses = {
        {0.1, -0.15, -0.1},
        {0.2, -0.2, -0.15},
        {0.15, -0.18, -0.12}
    };
    
    for (size_t i = 0; i < initial_guesses.size(); i++) {
        try {
            cout << "\n--- Попытка " << i+1 << " ---\n";
            vector<double> solution = NewtonSystem(F, J, initial_guesses[i], 10, 100, true);
            
            cout << "\nНайденное решение:\n";
            for (size_t j = 0; j < solution.size(); j++) {
                cout << "x" << j+1 << " = " << setprecision(10) << solution[j] << "\n";
            }
            
            // Проверка
            vector<double> check = F(solution);
            cout << "\nНевязки:\n";
            for (size_t j = 0; j < check.size(); j++) {
                cout << "f" << j+1 << " = " << scientific << check[j] << "\n";
            }
            cout << "----------------------------------------\n";
            
        } catch (const exception& e) {
            cerr << "Ошибка: " << e.what() << "\n";
        }
    }
}

void solveSystemC() {
    cout << "========================================\n";
    cout << "Система c (задача на собственные значения):\n";
    cout << "4x1 - x2 + x3 = x1*x4\n";
    cout << "-x1 + 3x2 - 2x3 = x2*x4\n";
    cout << "x1 - 2x2 + 3x3 = x3*x4\n";
    cout << "x1^2 + x2^2 + x3^2 = 1\n\n";
    
    auto F = [](const vector<double>& x) {
        vector<double> f(4);
        f[0] = 4.0*x[0] - x[1] + x[2] - x[0]*x[3];
        f[1] = -x[0] + 3.0*x[1] - 2.0*x[2] - x[1]*x[3];
        f[2] = x[0] - 2.0*x[1] + 3.0*x[2] - x[2]*x[3];
        f[3] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0;
        return f;
    };
    
    auto J = [](const vector<double>& x) {
        vector<vector<double>> Jx(4, vector<double>(4));
        
        Jx[0][0] = 4.0 - x[3];
        Jx[0][1] = -1.0;
        Jx[0][2] = 1.0;
        Jx[0][3] = -x[0];
        
        Jx[1][0] = -1.0;
        Jx[1][1] = 3.0 - x[3];
        Jx[1][2] = -2.0;
        Jx[1][3] = -x[1];
        
        Jx[2][0] = 1.0;
        Jx[2][1] = -2.0;
        Jx[2][2] = 3.0 - x[3];
        Jx[2][3] = -x[2];
        
        Jx[3][0] = 2.0*x[0];
        Jx[3][1] = 2.0*x[1];
        Jx[3][2] = 2.0*x[2];
        Jx[3][3] = 0.0;
        
        return Jx;
    };
    
    vector<vector<double>> initial_guesses = {
        {0.5, 0.5, 0.5, 2.0},
        {0.7, 0.3, -0.5, 4.0},
        {0.3, -0.7, 0.5, 6.0},
        {-0.5, 0.5, -0.5, 1.0}
    };
    
    for (size_t i = 0; i < initial_guesses.size(); i++) {
        try {
            cout << "\n--- Попытка " << i+1 << " ---\n";
            vector<double> solution = NewtonSystem(F, J, initial_guesses[i], 10, 100, true);
            
            cout << "\nНайденное решение (собственный вектор и собственное значение):\n";
            cout << "x1 = " << setprecision(10) << solution[0] << "\n";
            cout << "x2 = " << setprecision(10) << solution[1] << "\n";
            cout << "x3 = " << setprecision(10) << solution[2] << "\n";
            cout << "λ = x4 = " << setprecision(10) << solution[3] << "\n";
            
            double norm = sqrt(solution[0]*solution[0] + solution[1]*solution[1] + solution[2]*solution[2]);
            cout << "Норма вектора (должна быть 1): " << norm << "\n";
            
            vector<double> check = F(solution);
            cout << "\nНевязки:\n";
            for (size_t j = 0; j < check.size(); j++) {
                cout << "f" << j+1 << " = " << scientific << check[j] << "\n";
            }
            cout << "----------------------------------------\n";
            
        } catch (const exception& e) {
            cerr << "Ошибка: " << e.what() << "\n";
        }
    }
}

int main() {
    try {
        solveSystemA();
        cout << "\n";
        solveSystemB();
        cout << "\n";
        solveSystemC();
        
    } catch (const exception& e) {
        cerr << "Ошибка: " << e.what() << endl;
    }
    
    return 0;
}