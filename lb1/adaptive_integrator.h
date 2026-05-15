#ifndef ADAPTIVE_INTEGRATOR_H
#define ADAPTIVE_INTEGRATOR_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include "rossler.h"

template<typename Butcher>
void integrate_adaptive(double t0, double t_end, double y0[3], double h0, double tol, const std::string& filename) {
    std::ofstream out(filename);
    out << std::fixed << std::setprecision(8);
    out << "# t x y z\n";
    double t = t0;
    double y[3] = {y0[0], y0[1], y0[2]};
    double h = h0;
    int steps = 0, rejected = 0;
    constexpr int s = Butcher::stages;

    while (t < t_end) {
        out << t << " " << y[0] << " " << y[1] << " " << y[2] << "\n";
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
        double y_high[3], y_low[3];
        for (int dim = 0; dim < 3; ++dim) {
            y_high[dim] = y[dim];
            y_low[dim] = y[dim];
            for (int i = 0; i < s; ++i) {
                y_high[dim] += Butcher::b_high[i] * k[i][dim];
                y_low[dim]  += Butcher::b_low[i]  * k[i][dim];
            }
        }
        double max_err = 0.0;
        for (int dim = 0; dim < 3; ++dim) {
            double err = std::abs(y_high[dim] - y_low[dim]);
            if (err > max_err) max_err = err;
        }
        double new_h = 0.9 * h * std::pow(tol / max_err, 1.0/5.0);
        if (new_h > 4.0 * h) new_h = 4.0 * h;
        if (new_h < 0.25 * h) new_h = 0.25 * h;
        if (max_err <= tol) {
            for (int dim = 0; dim < 3; ++dim) y[dim] = y_high[dim];
            t += h;
            h = new_h;
            ++steps;
        } else {
            h = new_h;
            ++rejected;
        }
        if (h < 1e-15) break;
    }
    out.close();
    std::cout << "Steps accepted: " << steps << ", rejected: " << rejected << "\n";
}

#endif