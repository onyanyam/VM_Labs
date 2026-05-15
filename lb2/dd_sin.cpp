#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <iostream>
#include <iomanip>
#include <string>

class DDouble {
public:
    double hi;  // старшая часть (основная точность)
    double lo;  // младшая часть (коррекция), |l0| <= 0.5 ulp(hi)

    DDouble() : hi(0.0), lo(0.0) {}
    DDouble(double x) : hi(x), lo(0.0) {} 
    DDouble(double h, double l) {
    // Нормализация через Fast2Sum (стр. 126)
        hi = h + l;
        double v = hi - h;
        lo = (h - (hi - v)) + (l - v);
    }
};

// стр. 127-130: 2Sum
// Сложение двух double
DDouble two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    double e = (a - (s - v)) + (b - v);
    return DDouble(s, e);
}

// стр. 152: 2MultFMA (Fused Multiply-Add - уножение-сложение с однократным округлением)
// Умножение двух double -> DDouble
DDouble two_mul(double a, double b) {
    double p = a * b;
    double e = std::fma(a, b, -p);
    return DDouble(p, e);
}

// стр. 497: Сложение двух DD чисел
DDouble add_dd(const DDouble &x, const DDouble &y) {
    DDouble s = two_sum(x.hi, y.hi);   // (sh, sl)
    DDouble t = two_sum(x.lo, y.lo);   // (th, tl)
    double c = s.lo + t.hi;
    DDouble v(s.hi, c);
    double w = t.lo + v.lo;
    DDouble z(v.hi, w);
    return z;
}

// стр. 498: Умножение двух DD чисел
DDouble mul_dd(const DDouble &x, const DDouble &y) {
    DDouble p = two_mul(x.hi, y.hi);  // (ph, pl) - точное умножение hi*yhi
    double t1 = x.hi * y.lo;
    double t2 = x.lo * y.hi;
    double c = p.lo + t1 + t2;
    return DDouble(p.hi, c);           // Нормализация
}

// Вспомогательные функции
DDouble neg_dd(const DDouble &x) {
    return DDouble(-x.hi, -x.lo);
}

DDouble abs_dd(const DDouble &x) {
    return (x.hi < 0.0) ? neg_dd(x) : x;
}

DDouble div_dd_double(const DDouble &x, double a) {
    return DDouble(x.hi / a, x.lo / a);
}

DDouble quick_two_sum(double a, double b) {
    double s = a + b;
    double e = b - (s - a);
    return DDouble(s, e);
}

// Деление double-double на double-double с сохранением lo
DDouble div_dd(const DDouble &x, const DDouble &y) {
    double q1 = x.hi / y.hi;
    DDouble p = mul_dd(DDouble(q1), y);
    DDouble r = add_dd(x, neg_dd(p));
    double q2 = r.hi / y.hi;
    return quick_two_sum(q1, q2);  //  lo сохраняется
}

bool less_dd(const DDouble &x, const DDouble &y) {
    return x.hi < y.hi || (x.hi == y.hi && x.lo < y.lo);
}

bool eq_dd(const DDouble &x, const DDouble &y) {
    return x.hi == y.hi && x.lo == y.lo;
}

// Используем метод разделения константы (Cody and Waite, 379-380)
// 2/pi, разбитая на две части
double TWO_OVER_PI_HI = 0.63661977236758138;
double TWO_OVER_PI_LO = -3.3794354165343723e-17;
    
// pi/2, разбитое на две части
double PIO2_H0 = 1.5707963267948966e+00;
double PIO2_H1 = 6.1232339957367660e-17;

// Точное умножение double на DD-константу 2/pi
DDouble mul_power2(double a) {
    DDouble t = two_mul(a, TWO_OVER_PI_HI);
    t.lo += a * TWO_OVER_PI_LO;
    return DDouble(t.hi, t.lo); 
}

// Редукция аргумента: x -> r в [-pi/4, pi/4], возвращает k (квадрант)
int reduce_argument(const DDouble &x, DDouble &r) {

    // Шаг 1: k = round(x * 2/pi)
    DDouble prod = mul_power2(x.hi);
    double k = std::round(prod.hi + prod.lo);  // целое
    int k_int = static_cast<int>(k);
    
    // Шаг 2: Вычисляем x - k * pi/2 в DD точности
    // Умножаем k на каждую часть pi/2 и вычитаем из x
    DDouble t0 = two_mul(k, PIO2_H0);
    DDouble t1 = two_mul(k, PIO2_H1);
    
    DDouble s0 = two_sum(x.hi, -t0.hi);
    double s1 = x.lo - t0.lo - t1.hi - t1.lo;
    
    r = DDouble(s0.hi, s0.lo + s1); 
    return k_int;
}

// стр. 390-394: Полиномиальная аппроксимация
// Общая функция для вычисления полинома по схеме Горнера
// coeffs - массив коэффициентов (от старшей степени к младшей)
// n - количество коэффициентов
DDouble poly_eval(const DDouble &x, const double coeffs[], int n) {
    DDouble sum(coeffs[n-1]);
    for (int i = n-2; i >= 0; --i) {
        sum = add_dd(DDouble(coeffs[i]), mul_dd(sum, x));
    }
    return sum;
}

// P(t) = 1 - t/3! + t^2/5! - t^3/7! + ... - ряд Тейлора
DDouble sin_poly(const DDouble &r) {
    // Вычисляем r^2 в DD
    DDouble r2 = mul_dd(r, r);
    
    // Коэффициенты полинома для P(r^2)
    double coeffs[] = {
        1.00000000000000000e+00,
        -1.66666666666666657e-01,
        8.33333333333333322e-03,
        -1.98412698412698413e-04,
        2.75573192239858925e-06,
        -2.50521083854417202e-08,
        1.60590438368216133e-10,
        -7.64716373181981641e-13,
        2.81145725434552060e-15,
    };

    int n = sizeof(coeffs) / sizeof(coeffs[0]);
    
    // Вычисляем P(r^2) через схему Горнера
    DDouble P = poly_eval(r2, coeffs, n);
    
    // sin(r) = r * P(r^2)
    return mul_dd(r, P);
}

// Q(t) = 1 - t/2! + t^2/4! - t^3/6! + ... - ряд Тейлора
DDouble cos_poly(const DDouble &r) {
    DDouble r2 = mul_dd(r, r);
    
    // Коэффициенты полинома для cos(r) = Q(r^2)
    double coeffs[] = {
        1.00000000000000000e+00,
        -5.00000000000000000e-01,
        4.16666666666666644e-02,
        -1.38888888888888894e-03,
        2.48015873015873016e-05,
        -2.75573192239858883e-07,
        2.08767569878681002e-09,
        -1.14707455977297245e-11,
        4.77947733238738525e-14
    };

    int n = sizeof(coeffs) / sizeof(coeffs[0]);
    
    // cos(r) = Q(r^2)
    return poly_eval(r2, coeffs, n);
}

// Основная функция sin для double-double
DDouble dd_sin(const DDouble &x) {
    
    // Редукция аргумента: x -> r в [-pi/4, pi/4], k - номер квадранта
    DDouble r;
    int k = reduce_argument(x, r);
    
    // Вычисляем sin(r) и cos(r)
    DDouble sin_r = sin_poly(r);
    DDouble cos_r = cos_poly(r);
    
    // стр. 378: финальная сборка по квадранту
    switch (k & 3) {
        case 0: return sin_r;      // sin(x) = sin(r)
        case 1: return cos_r;      // sin(x) = cos(r)
        case 2: return neg_dd(sin_r); // sin(x) = -sin(r)
        case 3: return neg_dd(cos_r); // sin(x) = -cos(r)
        default: return sin_r;
    }
}

DDouble string_to_dd(const std::string& s) {
    size_t dotPos = s.find('.');
    if (dotPos == std::string::npos) {
        return DDouble(std::stod(s), 0.0);
    }

    std::string integral = s.substr(0, dotPos);
    std::string fractional = s.substr(dotPos + 1);
    
    double num = std::stod(integral + fractional);
    double den = std::pow(10.0, (double)fractional.length());
    
    return div_dd(DDouble(num, 0.0), DDouble(den, 0.0));
}

int main() {
    std::string input_str;
    std::cout << "Введите число: ";
    std::cin >> input_str;
    
    // Парсим строку напрямую в DDouble
    DDouble x = string_to_dd(input_str);
    
    // Вычисляем sin нашим методом
    DDouble res = dd_sin(x);
    
    // Для сравнения берём double-версию
    double x_double = std::stod(input_str);
    double std_sin = std::sin(x_double);

    DDouble exact = DDouble(res.hi, res.lo); 
    double exact_val = exact.hi + exact.lo;
    
    double diff = exact_val - std_sin;
    
    std::cout << std::scientific << std::setprecision(17);
    std::cout << "\n СРАВНЕНИЕ sin(x)" << std::endl;
    std::cout << "Входное x (hi): " << x.hi << std::endl;
    std::cout << "Входное x (lo): " << x.lo << std::endl;
    std::cout << "sin(x) (hi): " << res.hi << std::endl;
    std::cout << "sin(x) (lo): " << res.lo << std::endl;
    std::cout << "sin(x) (реализованный): " << exact_val << std::endl;
    std::cout << "sin(x) (стандартный):   " << std_sin << std::endl;
    std::cout << "Разница: " << diff << std::endl;
    std::cout << "Точность (|lo|/|hi|): " << std::abs(res.lo / res.hi) << std::endl;
    
    return 0;
}