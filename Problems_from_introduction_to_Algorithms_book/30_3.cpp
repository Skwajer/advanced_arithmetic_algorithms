#include "fft.hpp"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

std::vector<std::vector<std::complex<double>>> fft2d(
    std::vector<std::vector<std::complex<double>>>& a, 
    bool inverse = false)
{
    size_t n1 = a.size();
    size_t n2 = a[0].size();
    
    for (size_t i = 0; i < n1; i++) {
        a[i] = fft(a[i], inverse);
    }
    
    for (size_t j = 0; j < n2; j++) {
        std::vector<std::complex<double>> column(n1);
        for (size_t i = 0; i < n1; i++) {
            column[i] = a[i][j];
        }
        
        column = fft(column, inverse);
        
        for (size_t i = 0; i < n1; i++) {
            a[i][j] = column[i];
        }
    }
    
    return a;
}

std::vector<std::vector<std::vector<std::complex<double>>>> fft3d(
    std::vector<std::vector<std::vector<std::complex<double>>>>& a,
    bool inverse = false)
{
    size_t n1 = a.size();
    size_t n2 = a[0].size();
    size_t n3 = a[0][0].size();
    
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
            a[i][j] = fft(a[i][j], inverse);
        }
    }
    
    for (size_t i = 0; i < n1; i++) {
        for (size_t k = 0; k < n3; k++) {
            std::vector<std::complex<double>> slice(n2);
            for (size_t j = 0; j < n2; j++) {
                slice[j] = a[i][j][k];
            }
            
            slice = fft(slice, inverse);
            
            for (size_t j = 0; j < n2; j++) {
                a[i][j][k] = slice[j];
            }
        }
    }
    
    for (size_t j = 0; j < n2; j++) {
        for (size_t k = 0; k < n3; k++) {
            std::vector<std::complex<double>> slice(n1);
            for (size_t i = 0; i < n1; i++) {
                slice[i] = a[i][j][k];
            }
            
            slice = fft(slice, inverse);
            
            for (size_t i = 0; i < n1; i++) {
                a[i][j][k] = slice[i];
            }
        }
    }
    
    return a;
}

void print2d(const std::vector<std::vector<std::complex<double>>>& arr, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : arr) {
        for (const auto& val : row) {
            printf("%6.1f ", val.real());
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void print3d(const std::vector<std::vector<std::vector<std::complex<double>>>>& arr, const std::string& name) {
    std::cout << name << ":\n";
    for (size_t i = 0; i < arr.size(); i++) {
        std::cout << "Слой " << i << ":\n";
        for (size_t j = 0; j < arr[i].size(); j++) {
            for (size_t k = 0; k < arr[i][j].size(); k++) {
                printf("%6.1f ", arr[i][j][k].real());
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

int main() {
    {
        std::cout << "=== ТЕСТ 1: 2D БПФ (4x4) ===\n";
        
        std::vector<std::vector<std::complex<double>>> image = {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16}
        };
        
        print2d(image, "Исходное изображение");
        
        auto freq = fft2d(image, false);
        print2d(freq, "После 2D БПФ (частотная область)");
        
        auto restored = fft2d(freq, true);
        print2d(restored, "После обратного БПФ");
    }
    
    {
        std::cout << "=== ТЕСТ 2: 2D БПФ (3x4) с другим порядком ===\n";
        
        std::vector<std::vector<std::complex<double>>> arr = {
            {1, 1, 1, 1},
            {1, -1, 1, -1},
            {1, 1, -1, -1}
        };
        
        print2d(arr, "Исходный массив");
        
        auto arr1 = arr;
        for (size_t i = 0; i < arr1.size(); i++) {
            arr1[i] = fft(arr1[i], false);
        }
        for (size_t j = 0; j < arr1[0].size(); j++) {
            std::vector<std::complex<double>> col(arr1.size());
            for (size_t i = 0; i < arr1.size(); i++) {
                col[i] = arr1[i][j];
            }
            col = fft(col, false);
            for (size_t i = 0; i < arr1.size(); i++) {
                arr1[i][j] = col[i];
            }
        }
        
        auto arr2 = arr;
        for (size_t j = 0; j < arr2[0].size(); j++) {
            std::vector<std::complex<double>> col(arr2.size());
            for (size_t i = 0; i < arr2.size(); i++) {
                col[i] = arr2[i][j];
            }
            col = fft(col, false);
            for (size_t i = 0; i < arr2.size(); i++) {
                arr2[i][j] = col[i];
            }
        }
        for (size_t i = 0; i < arr2.size(); i++) {
            arr2[i] = fft(arr2[i], false);
        }
        
        std::cout << "Результат (сначала строки, потом столбцы):\n";
        for (const auto& row : arr1) {
            for (const auto& val : row) {
                printf("%6.1f ", val.real());
            }
            std::cout << "\n";
        }
        
        std::cout << "\nРезультат (сначала столбцы, потом строки):\n";
        for (const auto& row : arr2) {
            for (const auto& val : row) {
                printf("%6.1f ", val.real());
            }
            std::cout << "\n";
        }
        
        bool equal = true;
        for (size_t i = 0; i < arr1.size() && equal; i++) {
            for (size_t j = 0; j < arr1[i].size() && equal; j++) {
                if (std::abs(arr1[i][j].real() - arr2[i][j].real()) > 1e-10) {
                    equal = false;
                }
            }
        }
        std::cout << "\nРезультаты " << (equal ? "СОВПАДАЮТ" : "НЕ СОВПАДАЮТ") 
                  << " (порядок измерений не важен)\n\n";
    }
    
    {
        std::cout << "=== TEST 3: 3D БПФ (2x2x2) ===\n";
        
        std::vector<std::vector<std::vector<std::complex<double>>>> cube(2, 
            std::vector<std::vector<std::complex<double>>>(2, 
                std::vector<std::complex<double>>(2)));
        
        for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < 2; j++) {
                for (size_t k = 0; k < 2; k++) {
                    cube[i][j][k] = i + j + k + 1;
                }
            }
        }
        
        print3d(cube, "Исходный 3D массив");
        
        auto freq3d = fft3d(cube, false);
        print3d(freq3d, "После 3D БПФ");
        
        auto restored3d = fft3d(freq3d, true);
        print3d(restored3d, "После обратного 3D БПФ");
    }
    
    return 0;
}