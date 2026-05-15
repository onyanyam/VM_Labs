#ifndef FIXED_INTEGRATOR_H
#define FIXED_INTEGRATOR_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "rossler.h"

template<typename Butcher>
void integrate_fixed(double t0, double t_end, double y0[3], double h, const std::string& filename) {
    std::ofstream out(filename);
    out << std::fixed << std::setprecision(8);
    out << "# t x y z\n";
    double t = t0;
    double y[3] = {y0[0], y0[1], y0[2]};
    int steps = 0;
    constexpr int s = Butcher::stages;

    while (t <= t_end + 1e-12) {
        // out << t << " " << y[0] << " " << y[1] << " " << y[2] << "\n";
        double k[s][3];
        double ytmp[3];
        for (int i = 0; i < s; ++i) {
            if (i == 0) {
                rossler(t, y, k[0]);
            } else {
                for (int dim = 0; dim < 3; ++dim) {
                    ytmp[dim] = y[dim];
                    for (int j = 0; j < i; ++j)
                        ytmp[dim] += Butcher::a[i][j] * k[j][dim];
                }
                rossler(t + Butcher::c[i] * h, ytmp, k[i]);
            }
            for (int dim = 0; dim < 3; ++dim) k[i][dim] *= h;
        }
        for (int dim = 0; dim < 3; ++dim) {
            y[dim] = y[dim];
            for (int i = 0; i < s; ++i)
                y[dim] += Butcher::b[i] * k[i][dim];
        }
        t += h;
        ++steps;
    }
    out.close();
    // std::cout << "Steps: " << steps << "\n";
}
#endif