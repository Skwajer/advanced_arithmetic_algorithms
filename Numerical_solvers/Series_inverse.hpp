#include <stdexcept>
#include <vector>

static std::vector<double> invert_series(const std::vector<double>& a, int n) 
{
    std::vector<double> b(n, 0.0);
    double a0_inv = 1.0 / a[0];
    b[0] = a0_inv;

    for (int k = 1; k < n; ++k) 
    {
        double sum = 0.0;
        for (int i = 1; i <= k; ++i) 
        {
            if (i < a.size()) 
            {
                sum += a[i] * b[k - i];
            }
        }
        b[k] = -a0_inv * sum;
    }

    return b;
}