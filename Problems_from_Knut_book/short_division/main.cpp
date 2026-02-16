#include <stdio.h>
#include <stdlib.h>

#define MAX_DIGITS 100

void short_division(int u[], int n, int v, int b, int w[], int *r) {
    printf("\n=== SHORT DIVISION ALGORITHM DEMONSTRATION (Base b = %d) ===\n\n", b);
    printf("Dividend: ");
    for (int i = n-1; i >= 0; i--) {
        printf("%d ", u[i]);
    }
    printf(" (digits from least significant to most significant)\n");
    printf("Divisor (single precision): %d\n", v);
    printf("Number system base: %d\n\n", b);
    
    printf("--- Step-by-step execution ---\n\n");
    
    *r = 0;
    
    for (int j = n-1; j >= 0; j--) {
        int dividend = (*r) * b + u[j];
        
        w[j] = dividend / v;
        
        *r = dividend % v;
        
        printf("Step %d (digit %d):\n", n-j, j);
        printf("  Current dividend = (remainder r=%d) * %d + digit u[%d]=%d = %d\n", 
               (dividend - u[j]) / b, b, j, u[j], dividend);
        printf("  Quotient digit w[%d] = %d / %d = %d\n", j, dividend, v, w[j]);
        printf("  New remainder r = %d %% %d = %d\n\n", dividend, v, *r);
    }
    
    printf("--- Result ---\n");
    printf("Quotient (digits from least to most significant): ");
    for (int i = n-1; i >= 0; i--) {
        printf("%d ", w[i]);
    }
    printf("\nRemainder: %d\n", *r);
}

void print_number(int digits[], int n, char* name) {
    printf("%s: ", name);
    for (int i = n-1; i >= 0; i--) {
        printf("%d", digits[i]);
    }
    printf("\n");
}

int main() {
    int b;
    int n;
    int v;
    int u[MAX_DIGITS];
    int w[MAX_DIGITS];
    int r;
    
    printf("============================================\n");
    printf("SHORT DIVISION ALGORITHM (from Knuth, Vol.2, 4.3.1 Ex.16)\n");
    printf("============================================\n\n");
    
    printf("Enter the number system base b (e.g., 10): ");
    scanf("%d", &b);
    
    printf("Enter the number of digits n: ");
    scanf("%d", &n);
    
    printf("Enter the digits of the number from MOST significant to LEAST significant:\n");
    printf("(for example, for number 1234 enter: 1 2 3 4)\n");
    printf("Digits: ");
    for (int i = n-1; i >= 0; i--) {
        int digit;
        scanf("%d", &digit);
        u[i] = digit;
    }
    
    printf("Enter the divisor v (0 < v < %d): ", b);
    scanf("%d", &v);
    
    if (v <= 0 || v >= b) {
        printf("Error: divisor must satisfy 0 < v < %d\n", b);
        return 1;
    }
    
    for (int i = 0; i < n; i++) {
        if (u[i] < 0 || u[i] >= b) {
            printf("Error: digit u[%d] = %d is invalid for base %d\n", i, u[i], b);
            return 1;
        }
    }
    
    printf("\n--- Initial data ---\n");
    print_number(u, n, "Dividend");
    printf("Divisor: %d\n", v);
    printf("Base: %d\n", b);
    
    short_division(u, n, v, b, w, &r);
    
    printf("\n--- Verification ---\n");
    
    int check = 0;
    for (int i = n-1; i >= 0; i--) {
        check = check * b + w[i];
    }
    check = check * v + r;
    
    int original = 0;
    for (int i = n-1; i >= 0; i--) {
        original = original * b + u[i];
    }
    
    printf("Original dividend: %d\n", original);
    printf("Verification: quotient * divisor + remainder = %d\n", check);
    
    if (original == check) {
        printf("✓ RESULT IS CORRECT!\n");
    } else {
        printf("✗ ERROR IN CALCULATIONS!\n");
    }
    
    return 0;
}