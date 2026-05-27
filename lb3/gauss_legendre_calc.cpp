#include <iostream>
#include <iomanip>
#include <cmath>
#include "gauss_table.h"

using namespace std;

// Вычисление подынтегральной функции в точке
double integral_func(double theta, double k)
{
    double sin_theta = sin(theta);
    double k_sin = k * sin_theta;
    return 1.0 / sqrt(1.0 - k_sin * k_sin);
}

// Вычисление интеграла
double ellip_GaussLegendre(double k)
{
    const double a = 0.0;
    const double b = M_PI_2;
    const double factor = (b - a) / 2.0;
    const double mid = (a + b) / 2.0;

    vector<GaussPoint> table = getGaussLegendreTable();

    double result = 0.0;
    for (size_t i = 0; i < table.size(); ++i)
    {
        const GaussPoint &point = table[i];
        double theta = factor * point.x + mid;
        double f = integral_func(theta, k);
        result += point.w * f;
    }

    return result * factor;
}

int main()
{
    double k;

    cout << "Gauss-Legendre quadrature (n=100)" << endl;
    cout << "Enter k (0 <= k < 1): ";
    cin >> k;

    if (k < 0 || k >= 1)
    {
        cout << "Error: k must be in [0, 1)" << endl;
        return 1;
    }

    double result = ellip_GaussLegendre(k);

    cout << fixed << setprecision(16);
    cout << "K(" << k << ") = " << result << endl;

    return 0;
}