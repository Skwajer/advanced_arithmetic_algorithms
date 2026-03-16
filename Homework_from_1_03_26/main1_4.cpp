#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#include <vector>
#include <string>

// ============================================================
// ДЕМОНСТРАЦИЯ МЕТОДА ПРОСТОЙ ИТЕРАЦИИ
// ============================================================

/**
 * Метод простой итерации с подробным выводом
 * @param phi итерационная функция
 * @param x0 начальное приближение
 * @param eps точность
 * @param max_iters максимальное число итераций
 * @param param_desc описание параметров для вывода
 */
void simpleIterationDemo(const std::function<double(double)>& phi, 
                         double x0, 
                         double eps,
                         size_t max_iters,
                         const std::string& param_desc) {
    
    std::cout << "\n" << std::string(50, '-') << "\n";
    std::cout << param_desc << "\n";
    std::cout << "Начальное приближение: x0 = " << x0 << "\n";
    std::cout << std::string(50, '-') << "\n";
    
    double x = x0;
    double x_next;
    size_t iter = 0;
    
    std::cout << std::left << std::setw(5) << "iter" 
              << std::setw(15) << "x" 
              << std::setw(15) << "φ(x)" 
              << std::setw(15) << "|x - x_prev|" << "\n";
    std::cout << std::string(50, '-') << "\n";
    
    do {
        x_next = phi(x);
        double diff = std::abs(x_next - x);
        
        std::cout << std::fixed << std::setprecision(8);
        std::cout << std::left << std::setw(5) << iter 
                  << std::setw(15) << x 
                  << std::setw(15) << x_next 
                  << std::setw(15) << diff << "\n";
        
        if (diff < eps) {
            std::cout << std::string(50, '-') << "\n";
            std::cout << "СОШЛОСЬ! Корень ≈ " << x_next 
                      << " за " << iter+1 << " итераций\n";
            return;
        }
        
        x = x_next;
        iter++;
        
        if (iter > max_iters) {
            std::cout << std::string(50, '-') << "\n";
            std::cout << "НЕ СОШЛОСЬ за " << max_iters << " итераций\n";
            return;
        }
        
    } while (true);
}

double phiA(double x, double a, double b, double c) {
    return c + a * sin(x) * sin(x) + b * cos(x) * cos(x);
}

double phiADerivative(double x, double a, double b) {
    return (a - b) * sin(2 * x);
}

void demonstrateFunctionA() {
    std::cout << "\n\n";
    std::cout << "============================================================\n";
    std::cout << "ФУНКЦИЯ (a): φ(x) = c + a·sin²x + b·cos²x\n";
    std::cout << "УСЛОВИЕ СХОДИМОСТИ: |a - b| < 1\n";
    std::cout << "============================================================\n";
    
    
    {
        double a = 0.3, b = 0.8, c = 1.0;
        double max_deriv = std::abs(a - b);
        
        std::string desc = "ТЕСТ 1: |a-b| = " + std::to_string(max_deriv) + 
                          " < 1 (СХОДИТСЯ)\n"
                          "a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", c=" + std::to_string(c);
        
        auto phi = [a, b, c](double x) { return phiA(x, a, b, c); };
        
        simpleIterationDemo(phi, 0.0, 1e-8, 50, desc);
        simpleIterationDemo(phi, 2.0, 1e-8, 50, desc + " (x0=2)");
        simpleIterationDemo(phi, -1.5, 1e-8, 50, desc + " (x0=-1.5)");
    }
    
    {
        double a = 2.0, b = 1.0, c = 0.0;
        double max_deriv = std::abs(a - b);
        
        std::string desc = "ТЕСТ 2: |a-b| = " + std::to_string(max_deriv) + 
                          " = 1 (ГРАНИЧНЫЙ)\n"
                          "a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", c=" + std::to_string(c);
        
        auto phi = [a, b, c](double x) { return phiA(x, a, b, c); };
        
        simpleIterationDemo(phi, 1.0, 1e-8, 50, desc);
    }
    
    {
        double a = 3.0, b = 1.0, c = 0.0;
        double max_deriv = std::abs(a - b);
        
        std::string desc = "ТЕСТ 3: |a-b| = " + std::to_string(max_deriv) + 
                          " > 1 (РАСХОДИТСЯ)\n"
                          "a=" + std::to_string(a) + ", b=" + std::to_string(b) + ", c=" + std::to_string(c);
        
        auto phi = [a, b, c](double x) { return phiA(x, a, b, c); };
        
        simpleIterationDemo(phi, 1.0, 1e-8, 20, desc);
    }
}


double phiB(double x, double a, double b, double c) {
    return c + a * exp(-b * x * x);
}

double phiBDerivative(double x, double a, double b) {
    return -2 * a * b * x * exp(-b * x * x);
}

double maxPhiBDerivative(double a, double b) {
    if (b <= 0) return 1e10;
    return std::abs(a) * sqrt(2 * b / M_E);
}

void demonstrateFunctionB() {
    std::cout << "\n\n";
    std::cout << "============================================================\n";
    std::cout << "ФУНКЦИЯ (b): φ(x) = c + a·e^{-b·x²}\n";
    std::cout << "УСЛОВИЕ СХОДИМОСТИ: b > 0 И |a| < √(e/(2b))\n";
    std::cout << "============================================================\n";
    
    {
        double a = 1.0, b = 2.0, c = 0.0;
        double limit = sqrt(M_E / (2 * b));
        double max_deriv = maxPhiBDerivative(a, b);
        
        std::string desc = "ТЕСТ 1: b=" + std::to_string(b) + " > 0, |a|=" + 
                          std::to_string(std::abs(a)) + " > √(e/(2b))=" + 
                          std::to_string(limit) + "\nmax|φ'|=" + 
                          std::to_string(max_deriv) + " > 1 (РАСХОДИТСЯ)";
        
        auto phi = [a, b, c](double x) { return phiB(x, a, b, c); };
        
        simpleIterationDemo(phi, 0.5, 1e-8, 50, desc);
        simpleIterationDemo(phi, 2.0, 1e-8, 50, desc + " (x0=2)");
    }
    
    {
        double a = 0.5, b = 1.0, c = 1.0;
        double limit = sqrt(M_E / (2 * b));
        double max_deriv = maxPhiBDerivative(a, b);
        
        std::string desc = "ТЕСТ 2: b=" + std::to_string(b) + " > 0, |a|=" + 
                          std::to_string(std::abs(a)) + " < √(e/(2b))=" + 
                          std::to_string(limit) + "\nmax|φ'|=" + 
                          std::to_string(max_deriv) + " < 1 (СХОДИТСЯ)";
        
        auto phi = [a, b, c](double x) { return phiB(x, a, b, c); };
        
        simpleIterationDemo(phi, -1.0, 1e-8, 50, desc);
    }
    
    {
        double b = 1.0;
        double a = sqrt(M_E / (2 * b));
        double c = 0.0;
        double max_deriv = maxPhiBDerivative(a, b);
        
        std::string desc = "ТЕСТ 3: b=" + std::to_string(b) + " > 0, |a| = √(e/(2b))=" + 
                          std::to_string(a) + "\nmax|φ'|=" + 
                          std::to_string(max_deriv) + " = 1 (ГРАНИЧНЫЙ)";
        
        auto phi = [a, b, c](double x) { return phiB(x, a, b, c); };
        
        simpleIterationDemo(phi, 1.0, 1e-8, 30, desc);
    }
    
    {
        double b = 1.0;
        double a = 2.0; // > 1.165
        double c = 0.0;
        double limit = sqrt(M_E / (2 * b));
        double max_deriv = maxPhiBDerivative(a, b);
        
        std::string desc = "ТЕСТ 4: b=" + std::to_string(b) + " > 0, |a|=" + 
                          std::to_string(std::abs(a)) + " > √(e/(2b))=" + 
                          std::to_string(limit) + "\nmax|φ'|=" + 
                          std::to_string(max_deriv) + " > 1 (РАСХОДИТСЯ)";
        
        auto phi = [a, b, c](double x) { return phiB(x, a, b, c); };
        
        simpleIterationDemo(phi, 0.5, 1e-8, 15, desc);
    }
    
    {
        double a = 1.0, b = -0.5, c = 0.0;
        
        std::string desc = "ТЕСТ 5: b=" + std::to_string(b) + " < 0 (ЗАВЕДОМО РАСХОДИТСЯ)";
        
        auto phi = [a, b, c](double x) { return phiB(x, a, b, c); };
        
        simpleIterationDemo(phi, 0.1, 1e-8, 50, desc);
    }
}


int main() {
    std::cout << std::setprecision(8);
    
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "ДЕМОНСТРАЦИЯ МЕТОДА ПРОСТОЙ ИТЕРАЦИИ\n";
    std::cout << "Исследование сходимости для различных значений параметров\n";
    std::cout << "================================================================\n";
    

    demonstrateFunctionA();
    
    demonstrateFunctionB();
    
    std::cout << "\n\n";
    std::cout << "================================================================\n";
    std::cout << "ИТОГОВЫЕ УСЛОВИЯ СХОДИМОСТИ:\n";
    std::cout << "================================================================\n";
    std::cout << "(a) φ(x) = c + a·sin²x + b·cos²x\n";
    std::cout << "    |a - b| < 1\n\n";
    std::cout << "(b) φ(x) = c + a·e^{-b·x²}\n";
    std::cout << "    b > 0  И  |a| < √(e/(2b))\n";
    std::cout << "================================================================\n";
    
    return 0;
}