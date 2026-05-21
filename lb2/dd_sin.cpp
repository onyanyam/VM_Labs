#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>

class DDouble {
public:
    double hi, lo;

    DDouble() : hi(0.0), lo(0.0) {}
    DDouble(double x) : hi(x), lo(0.0) {}
    DDouble(double h, double l) {
        // Нормализация через two_sum
        double s = h + l;
        double v = s - h;
        double e = (h - (s - v)) + (l - v);
        hi = s;
        lo = e;
    }
};

// 2/pi  - разбиение на 4 части
static const double two_div_pi[] = {
    0.636619772367581343,
    5.912351290886566e-18,
    -6.046162238423234e-35,
    1.146313195287612e-51
};

// pi/2 - разбиение на 5 частей
static const double pi_div_two[] = {
    1.570796326794896558e+00, 
    6.123233995736766036e-17,
    2.022266248795950632e-33, 
    1.996755123683625902e-49,
    6.037002664545628121e-66
};

// константы pi/2, pi/4
static const DDouble DD_PI_2(1.57079632679489656e+00, 6.12323399573676604e-17);
static const DDouble DD_PI_4(7.85398163397448279e-01, 3.06161699786838302e-17);

// коэффициенты синуса (до 1/41!)
static const DDouble sin_coeffs[] = {
    { 1.00000000000000000e+00, 0.00000000000000000e+00 },
    { -1.66666666666666657e-01, -9.25185853854297066e-18 },
    { 8.33333333333333322e-03, 1.15648231731787138e-19 },
    { -1.98412698412698413e-04, -1.72095582934207053e-22 },
    { 2.75573192239858925e-06, -1.85839327404647208e-22 },
    { -2.50521083854417202e-08, 1.44881407093591197e-24 },
    { 1.60590438368216133e-10, 1.25852945887520981e-26 },
    { -7.64716373181981641e-13, -7.03872877733453001e-30 },
    { 2.81145725434552060e-15, 1.65088427308614326e-31 },
    { -8.22063524662432950e-18, -2.21418941196042654e-34 },
    { 1.95729410633912626e-20, -1.36435038300879085e-36 },
    { -3.86817017063068413e-23, 8.84317765548234385e-40 },
    { 6.44695028438447359e-26, -1.93304042337034648e-42 },
    { -9.18368986379554601e-29, -1.43031503967873220e-45 },
    { 1.13099628864477159e-31, 1.04980154129595060e-47 },
    { -1.21612504155351789e-34, -5.58629056788880577e-51 },
    { 1.15163356207719509e-37, -6.09957445788453978e-54 },
    { -9.67759295863189067e-41, -3.20229554864556196e-57 },
    { 7.26546017915307136e-44, -4.36409714935444569e-61 },
    { -4.90246975651354352e-47, 1.21301910051792795e-63 },
    { 2.98931082714240461e-50, -1.04072477030331555e-66 },
};
    
static const int SIN_COEFFS_COUNT = sizeof(sin_coeffs) / sizeof(sin_coeffs[0]);

// коэффициенты косинуса (до 1/40! )
static const DDouble cos_coeffs[] = {
    { 1.00000000000000000e+00, 0.00000000000000000e+00 },
    { -5.00000000000000000e-01, -0.00000000000000000e+00 },
    { 4.16666666666666644e-02, 2.31296463463574266e-18 },
    { -1.38888888888888894e-03, 5.30054395437357706e-20 },
    { 2.48015873015873016e-05, 2.15119478667758816e-23 },
    { -2.75573192239858883e-07, -2.37677146222502973e-23 },
    { 2.08767569878681002e-09, -1.20734505911325997e-25 },
    { -1.14707455977297245e-11, -2.06555127528307454e-28 },
    { 4.77947733238738525e-14, 4.39920548583408126e-31 },
    { -1.56192069685862253e-16, -1.19106796602737540e-32 },
    { 4.11031762331216484e-19, 1.44129733786595271e-36 },
    { -8.89679139245057408e-22, 7.91140261487237622e-38 },
    { 1.61173757109611839e-24, -3.68465735645097660e-41 },
    { -2.47959626322479759e-27, 1.29537309647652288e-43 },
    { 3.27988923706983776e-30, 1.51175427440298787e-46 },
    { -3.76998762881590539e-33, -2.58703478327503238e-49 },
    { 3.80039075485474342e-36, 1.74571580246525180e-52 },
    { -3.38715753552116180e-39, -5.09056148151084995e-56 },
    { 2.68822026628663633e-42, 5.35506116594333401e-59 },
    { -1.91196320504028195e-45, 2.78608221768831261e-62 },
    { 1.22561743912838585e-48, 6.03392734831560539e-68 },
};

static const int COS_COEFFS_COUNT = sizeof(cos_coeffs) / sizeof(cos_coeffs[0]);


// Вспомогательные операции с double
// точное сложение двух double с остатком
inline DDouble two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    double e = (a - (s - v)) + (b - v);
    return DDouble(s, e);
}

// умножение через FMA
inline DDouble two_mul(double a, double b) {
    double p = a * b;
    double e = std::fma(a, b, -p);
    return DDouble(p, e);
}

// быстрое сложение при |a|>=|b
inline DDouble quick_two_sum(double a, double b) {
    double s = a + b;
    double e = b - (s - a);
    return DDouble(s, e);
}

// быстрое сложение
inline std::pair<double, double> fast_two_sum(double a, double b) {
    double s = a + b;
    double z = s - a;
    double t = b - z;
    return {s, t};
}

// точное сложение с компенсацией
inline std::pair<double, double> is_two_sum(double a, double b) {
    double s = a + b;
    double a_ = s - b;
    double b_ = s - a_;
    double da = a - a_;
    double db = b - b_;
    double t = da + db;
    return {s, t};
}

// разделение Велткампа
inline std::pair<double, double> split(double x) {
    const double C = 134217729.0;
    double gamma = C * x;
    double delta = x - gamma;
    double x_hi = gamma + delta;
    double x_lo = x - x_hi;
    return {x_hi, x_lo};
}

// точное произведение Деккера
inline std::pair<double, double> dekker(double x, double y) {
    auto [x_hi, x_lo] = split(x);
    auto [y_hi, y_lo] = split(y);
    double p_hi = x * y;
    double t1 = -p_hi + x_hi * y_hi;
    double t2 = t1 + x_hi * y_lo;
    double t3 = t2 + x_lo * y_hi;
    double p_lo = t3 + x_lo * y_lo;
    return {p_hi, p_lo};
}

// Основные операции
// сложение двух double-double чисел
DDouble add_dd(const DDouble &x, const DDouble &y) {
    auto [sh, sl] = is_two_sum(x.hi, y.hi);
    auto [th, tl] = is_two_sum(x.lo, y.lo);
    double c = sl + th;
    auto [vh, vl] = fast_two_sum(sh, c);
    double w = tl + vl;
    auto [zh, zl] = fast_two_sum(vh, w);
    return DDouble(zh, zl);
}

// умножение двух double-double чисел
DDouble mul_dd(const DDouble &x, const DDouble &y) {
    auto [ph, pl] = dekker(x.hi, y.hi);
    pl += x.hi * y.lo + x.lo * y.hi;
    auto [zh, zl] = fast_two_sum(ph, pl);
    return DDouble(zh, zl);
}

// нормализация трёх чисел в нормализованную пару head+tail
DDouble normalize3(double a0, double a1, double a2) {
    auto [s, t2] = fast_two_sum(a1, a2);
    auto [t0, t1] = fast_two_sum(a0, s);
    auto [b0, s1] = fast_two_sum(t0, t1);
    auto [b0_final, b1] = fast_two_sum(b0, s1 + t2);
    return DDouble(b0_final, b1);
}

// смена знака double-double
DDouble neg_dd(const DDouble &x) {
    return DDouble(-x.hi, -x.lo);
}

// деление double-double на double-double
DDouble div_dd(const DDouble &x, const DDouble &y) {
    double q1 = x.hi / y.hi;
    DDouble p = mul_dd(DDouble(q1), y);
    DDouble r = add_dd(x, neg_dd(p));
    double q2 = r.hi / y.hi;
    DDouble res = quick_two_sum(q1, q2);
    return normalize3(res.hi, res.lo, 0.0);
}

// вычисление sin(r) через ряд Тейлора для малых r (редуцированных)
DDouble sin_poly(const DDouble &r) {
    DDouble r2 = mul_dd(r, r);
    DDouble term = r;
    DDouble sum = r;
    for (int i = 1; i < SIN_COEFFS_COUNT; ++i) {
        term = mul_dd(term, r2);
        DDouble delta = mul_dd(term, sin_coeffs[i]);
        sum = add_dd(sum, delta);
        sum = normalize3(sum.hi, sum.lo, 0.0);
        if (fabs(delta.hi) < fabs(sum.hi) * 1e-32) break;
    }
    return sum;
}

// вычисление cos(r) через ряд Тейлора для малых r (редуцированных)
DDouble cos_poly(const DDouble &r) {
    DDouble r2 = mul_dd(r, r);
    DDouble term = DDouble(1.0);
    DDouble sum = DDouble(1.0);
    for (int i = 1; i < COS_COEFFS_COUNT; ++i) {
        term = mul_dd(term, r2);
        DDouble delta = mul_dd(term, cos_coeffs[i]);
        sum = add_dd(sum, delta);
        sum = normalize3(sum.hi, sum.lo, 0.0);
        if (fabs(delta.hi) < fabs(sum.hi) * 1e-32) break;
    }
    return sum;
}

// приведение аргумента к [-pi/4, pi/4] с определением квадранта
int reduce_argument(const DDouble &x, DDouble &r) {
    double x_hi = x.hi, x_lo = x.lo;
    
    double k = std::round(x_hi * two_div_pi[0]);
    k += std::round(x_hi * two_div_pi[1] + x_lo * two_div_pi[0]);
    k += std::round(x_hi * two_div_pi[2] + x_lo * two_div_pi[1]);
    k += std::round(x_hi * two_div_pi[3] + x_lo * two_div_pi[2]);

    long long k_int = std::llround(k);
    double k1 = std::floor(static_cast<double>(k_int) * 1e-9) * 1e9;
    double k2 = static_cast<double>(k_int) - k1;
    
    DDouble r_val = x;
    double parts_k[2] = {k1, k2};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 5; ++j) {
            DDouble term = mul_dd(DDouble(parts_k[i]), DDouble(pi_div_two[j]));
            r_val = add_dd(r_val, neg_dd(term));
        }
    }
    
    const DDouble PI4(DD_PI_4.hi, DD_PI_4.lo);
    const DDouble HALF_PI(DD_PI_2.hi, DD_PI_2.lo);

    int quad = k_int & 3;

   while (r_val.hi > PI4.hi || (r_val.hi == PI4.hi && r_val.lo > PI4.lo)) {
        r_val = add_dd(r_val, neg_dd(HALF_PI));
        quad = (quad + 1) % 4;
    }
    while (r_val.hi < -PI4.hi || (r_val.hi == -PI4.hi && r_val.lo < -PI4.lo)) {
        r_val = add_dd(r_val, HALF_PI);
        quad = (quad + 3) % 4;
    }
    
    r = r_val;
    return quad;
}

// основная функция синуса с учётом знака и квадранта
DDouble dd_sin(const DDouble &x) {

    if (std::isnan(x.hi) || std::isinf(x.hi)) {
        return DDouble(std::nan(""), 0.0);
    }

    // запоминаем знак исходного аргумента
    bool neg = std::signbit(x.hi) || (x.hi == 0.0 && std::signbit(x.lo));
    DDouble abs_x = neg ? neg_dd(x) : x;
    
    DDouble r;
    int quad = reduce_argument(abs_x, r);
    DDouble sin_r = sin_poly(r);
    DDouble cos_r = cos_poly(r);
    
    DDouble res;
    switch (quad) {
        case 0: res = sin_r; break;
        case 1: res = cos_r; break;
        case 2: res = neg_dd(sin_r); break;
        case 3: res = neg_dd(cos_r); break;
        default: res = sin_r;
    }
    if (neg) res = neg_dd(res);
    return res;
}

// вспомогательный ввод
// НЕ поддерживает экспоненциальную форму
DDouble string_to_dd(const std::string &s) {
    size_t dot = s.find('.');
    if (dot == std::string::npos) return DDouble(std::stod(s), 0.0);    
    std::string int_part = s.substr(0, dot);
    std::string frac_part = s.substr(dot + 1);
    double num = std::stod(int_part + frac_part);
    double den = std::pow(10.0, static_cast<double>(frac_part.size()));
    return div_dd(DDouble(num), DDouble(den));
}

int main() {
    std::string input;
    std::cout << "Введите число: ";
    std::cin >> input;
    DDouble x = string_to_dd(input);
    DDouble res = dd_sin(x);
    double std_sin = std::sin(std::stod(input));

    std::cout << std::scientific << std::setprecision(17);
    std::cout << "x (hi, lo): " << x.hi << " , " << x.lo << "\n";
    std::cout << "sin(x) (hi, lo): " << res.hi << " , " << res.lo << "\n";
    std::cout << "sin(x) (реализованный): " << (res.hi + res.lo) << "\n";
    std::cout << "sin(x) (стандартный):   " << std_sin << "\n";
    std::cout << "Разница: " << (res.hi + res.lo - std_sin) << "\n";
    std::cout << "|lo/hi|: " << std::abs(res.lo / res.hi) << "\n";
    return 0;
}