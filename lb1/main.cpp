#include <iostream>
#include <string>
#include "rossler.h"

#include "rk4.h"
#include "rk38.h"
#include "ralston3.h"
#include "rkf45.h"
#include "dopri45.h"

#include "fixed_integrator.h"
#include "adaptive_integrator.h"

int main() {
    std::cout << "Rossler system parameters a, b, c: ";
    std::cin >> a >> b >> c;
    double x0, y0, z0;
    std::cout << "Initial conditions x0, y0, z0: ";
    std::cin >> x0 >> y0 >> z0;
    double t_end;
    std::cout << "Integration time t_end: ";
    std::cin >> t_end;
    double y0_arr[3] = {x0, y0, z0};

    std::cout << "Choose method:\n"
              << "1 - RK4\n2 - RK3/8\n3 - Ralston3\n4 - RKF45\n5 - DOPRI45\n";
    int choice;
    std::cin >> choice;
    std::string outfile;
    std::cout << "Output filename: ";
    std::cin >> outfile;

    if (choice == 1) {
        double h;
        std::cout << "Fixed step size h: ";
        std::cin >> h;
        integrate_fixed<RK4>(0.0, t_end, y0_arr, h, outfile);
    } else if (choice == 2) {
        double h;
        std::cout << "Fixed step size h: ";
        std::cin >> h;
        integrate_fixed<RK38>(0.0, t_end, y0_arr, h, outfile);
    } else if (choice == 3) {
        double h;
        std::cout << "Fixed step size h: ";
        std::cin >> h;
        integrate_fixed<Ralston3>(0.0, t_end, y0_arr, h, outfile);
    } else if (choice == 4) {
        double tol, h0;
        std::cout << "Tolerance: ";
        std::cin >> tol;
        std::cout << "Initial step size h0: ";
        std::cin >> h0;
        integrate_adaptive<RKF45>(0.0, t_end, y0_arr, h0, tol, outfile);
    } else if (choice == 5) {
        double tol, h0;
        std::cout << "Tolerance: ";
        std::cin >> tol;
        std::cout << "Initial step size h0: ";
        std::cin >> h0;
        integrate_adaptive<DOPRI45>(0.0, t_end, y0_arr, h0, tol, outfile);
    } else {
        std::cerr << "Invalid choice\n";
        return 1;
    }

    std::cout << "Results saved to " << outfile << std::endl;
    return 0;
}