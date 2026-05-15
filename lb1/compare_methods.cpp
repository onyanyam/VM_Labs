#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

#include "rossler.h"
#include "rk4.h"
#include "rk38.h"
#include "ralston3.h"
#include "rkf45.h"
#include "dopri45.h"
#include "fixed_integrator.h"
#include "adaptive_integrator.h"

// Глобальные параметры системы 
extern double a, b, c;

using namespace std;

// Структура для хранения одной точки траектории
struct Point3D {
    double t;
    double x, y, z;
};

// Чтение данных из файла, созданного интегратором
vector<Point3D> load_points(const string& filename) {
    vector<Point3D> points;
    ifstream in(filename);

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        double t, x, y, z;
        istringstream iss(line);
        iss >> t >> x >> y >> z;
        points.push_back({t, x, y, z});
    }
    return points;
}

// Линейная интерполяция для нахождения значения в момент t
Point3D interpolate(const vector<Point3D>& ref, double t) {
    if (t <= ref.front().t) return ref.front();
    if (t >= ref.back().t) return ref.back();
    
    size_t i = 0;
    while (i+1 < ref.size() && ref[i+1].t < t) ++i;
    double t1 = ref[i].t, t2 = ref[i+1].t;
    double alpha = (t - t1) / (t2 - t1);
    double x = ref[i].x + alpha * (ref[i+1].x - ref[i].x);
    double y = ref[i].y + alpha * (ref[i+1].y - ref[i].y);
    double z = ref[i].z + alpha * (ref[i+1].z - ref[i].z);
    return {t, x, y, z};
}


// Вычисление ошибок между тестовым решением и эталоном
void compute_errors(const vector<Point3D>& test, const vector<Point3D>& ref,
                    double& abs_err, double& rel_err) {
    abs_err = 0.0;
    rel_err = 0.0;
    size_t n = test.size();
    if (n == 0) return;
    
    for (const auto& p : test) {
        Point3D ref_p = interpolate(ref, p.t);
        double dx = p.x - ref_p.x;
        double dy = p.y - ref_p.y;
        double dz = p.z - ref_p.z;
        double abs = sqrt(dx*dx + dy*dy + dz*dz);
        double rel = abs / (sqrt(ref_p.x*ref_p.x + ref_p.y*ref_p.y + ref_p.z*ref_p.z) + 1e-12);
        abs_err += abs;
        rel_err += rel;
    }
}

// Запуск фиксированного метода (RK4, RK38, Ralston3) и запись в файл
void run_fixed_method(const string& method_name, double h, const string& outfile,
                      double t0, double t_end, double y0[3]) {
    if (method_name == "RK4") {
        integrate_fixed<RK4>(t0, t_end, y0, h, outfile);
    } else if (method_name == "RK38") {
        integrate_fixed<RK38>(t0, t_end, y0, h, outfile);
    } else if (method_name == "Ralston3") {
        integrate_fixed<Ralston3>(t0, t_end, y0, h, outfile);
    }
}

// Запуск адаптивного метода (RKF45, DOPRI45) и запись в файл
void run_adaptive_method(const string& method_name, double tol, double h0,
                         const string& outfile, double t0, double t_end, double y0[3]) {
    if (method_name == "RKF45") {
        integrate_adaptive<RKF45>(t0, t_end, y0, h0, tol, outfile);
    } else if (method_name == "DOPRI45") {
        integrate_adaptive<DOPRI45>(t0, t_end, y0, h0, tol, outfile);
    }
}

int main() {
    cout << "COMPARISON OF NUMERICAL METHODS FOR ROSSLER SYSTEM" << endl;
    
    // Ввод параметров системы
    cout << "Enter Rossler parameters a, b, c:\n";
    cin >> a >> b >> c;
    
    // Начальные условия
    double x0, y0, z0;
    cout << "Enter initial conditions x0, y0, z0:\n";
    cin >> x0 >> y0 >> z0;
    double y0_arr[3] = {x0, y0, z0};
    
    double t_end;
    cout << "Enter integration time t_end:\n";
    cin >> t_end;
    
    // Генерация эталонного решения (DOPRI45 )
    cout << "\nGenerating reference solution (DOPRI45, tolerance = 1e-12)" << endl;
    double ref_tol = 1e-12;
    double ref_h0 = 1e-5;   // начальный шаг
    string ref_file = "DOPRI45_reference.txt";
    integrate_adaptive<DOPRI45>(0.0, t_end, y0_arr, ref_h0, ref_tol, ref_file);
    vector<Point3D> ref_points = load_points(ref_file);
    cout << "Reference points: " << ref_points.size() << endl;
    
    // Определение тестируемых методов и их параметров
    struct TestCase {
        string name;
        bool is_adaptive;
        double param;   // для фиксированных – h, для адаптивных – tolerance
        double h0;      // для адаптивных – начальный шаг
    };
    
    vector<TestCase> tests;
    
    // Фиксированные методы
    cout << "\nEnter fixed step size for RK4, RK38, Ralston3 (e.g., 1e-3):\n";
    double h_fixed;
    cin >> h_fixed;
    tests.push_back({"RK4", false, h_fixed, 0.0});
    tests.push_back({"RK38", false, h_fixed, 0.0});
    tests.push_back({"Ralston3", false, h_fixed, 0.0});
    
    // Адаптивные методы
    cout << "Enter tolerance for RKF45 (e.g., 1e-8):\n";
    double tol_rkf;
    cin >> tol_rkf;
    cout << "Enter initial step for RKF45 (e.g., 1e-5):\n";
    double h0_rkf;
    cin >> h0_rkf;
    tests.push_back({"RKF45", true, tol_rkf, h0_rkf});
    
    cout << "Enter tolerance for DOPRI45 (e.g., 1e-8):\n";
    double tol_dop;
    cin >> tol_dop;
    cout << "Enter initial step for DOPRI45 (e.g., 1e-5):\n";
    double h0_dop;
    cin >> h0_dop;
    tests.push_back({"DOPRI45", true, tol_dop, h0_dop});
    
    // Запуск каждого метода и сравнение
    cout << "\nRunning comparisons...\n" << endl;
    cout << fixed << setprecision(10);

    cout << left << setw(15) << "Method"
         << right << setw(20) << "Absolute error"
         << setw(20) << "Relative error" << endl;
    cout << string(55, '-') << endl;
    
    for (const auto& test : tests) {
        string outfile = test.name + "_temp.txt";
        
        // Запуск метода
        if (!test.is_adaptive) {
            run_fixed_method(test.name, test.param, outfile, 0.0, t_end, y0_arr);
        } else {
            run_adaptive_method(test.name, test.param, test.h0, outfile, 0.0, t_end, y0_arr);
        }
        
        // Загрузка полученных точек
        vector<Point3D> test_points = load_points(outfile);
        
        // Вычисление ошибок
        double abs_err = 0.0;
        double rel_err = 0.0;
        compute_errors(test_points, ref_points, abs_err, rel_err);
        
        // Вывод результата
        cout << left << setw(15) << test.name
             << right << setw(20) << abs_err
             << setw(20) << rel_err << endl;
    }
    
    cout << "\nComparison completed.\n";
    return 0;
}