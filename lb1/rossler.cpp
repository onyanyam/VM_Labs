#include "rossler.h"

double a, b, c;

void rossler(double t, const double variables[], double derivatives[]) {
    derivatives[0] = -variables[1] - variables[2];
    derivatives[1] = variables[0] + a * variables[1];
    derivatives[2] = b + variables[2] * (variables[0] - c);
}