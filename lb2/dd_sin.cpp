#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

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

// Предвычисленные коэффициенты рядов + разбиение констант для редукции аргумента в DDouble 
// Разбиение 2/pi на 4 части
static const double two_div_pi[] = {
    0.636619772367581343,
    5.912351290886566e-18,
    -6.046162238423234e-35,
    1.146313195287612e-51
};

// Разбиение pi/2 на 5 частей
static const double pi_div_two[] = {
    1.570796326794896558e+00, 6.123233995736766036e-17,
    2.022266248795950632e-33, 1.996755123683625902e-49,
    6.037002664545628121e-66
};

// Коэфициенты синуса для ряда Тейлора
static const DDouble sin_coeffs[]{
    {1.0, 0.0},
    {-1.666666666666666574e-01,-9.251858538542971261e-18},
    {8.3333333333333332176e-03, 1.1564823173178715840e-19},
    {-1.98412698412698412e-04,-1.72095582934207053e-22},
    {2.755731922398589251e-06,-1.85839327404647234e-22},
    {-2.50521083854417202e-08, 1.448814070935912168e-24},
    {1.605904383682161334e-10, 1.258529458875209764e-26},
    {-7.64716373181981640e-13, -7.03872877733453034e-30},
    {2.811457254345520598e-15, 1.650884273086143474e-31},
    {-8.22063524662432971e-18, -2.21418941196042675e-34},
    {1.957294106339126123e-20, -1.36435038300879101e-36},
    {-3.86817017063068412e-23, 8.843177655482344715e-40},
    {6.446950284384473589e-26, -1.93304042337034673e-42},
    {-9.18368986379554600e-29, -1.43031503967873239e-45},
    {1.130996288644771588e-31, 1.049801541295950759e-47},
    {-1.21612504155351789e-34, -5.58629056788880660e-51},
    {1.151633562077195089e-37, -6.09957445788454050e-54},
    {-9.67759295863189067e-41, -3.20229554864556237e-57},
    {7.265460179153071359e-44, -4.36409714935444624e-61},
    {-4.90246975651354351e-47, 1.213019100517928105e-63},
    {2.989310827142404614e-50, -1.040724770303315688e-66}
    };
    
static const int SIN_COEFFS_COUNT = sizeof(sin_coeffs) / sizeof(sin_coeffs[0]);

// Коэфициенты косинуса для ряда Тейлора
static const DDouble cos_coeffs[] = {
    { 1.0, 0.0 },
    { -5.000000000000000000e-01,  0.000000000000000000e+00 },
    {  4.166666666666666435e-02,  2.312964634635742664e-18 },
    { -1.388888888888888941e-03, -5.300543954373577065e-20 }, 
    {  2.480158730158730159e-05,  2.151194786677588162e-23 },
    { -2.755731922398588832e-07, -2.376771462225029730e-23 },
    {  2.087675698786810015e-09, -1.207345059113259965e-25 },
    { -1.147074559772972448e-11,  2.065551275283074542e-28 }, 
    {  4.779477332387385250e-14,  4.399205485834081263e-31 },
    { -1.554312234406214050e-16, -1.748367015502621006e-33 },
    {  4.116161026367464973e-19,  2.296734139988225573e-36 },
    { -8.866726349929881407e-22, -1.272023531649983177e-39 },
    {  1.606290998900340866e-24, -1.258670417726487532e-41 },
    { -2.471216921385139794e-27, -2.064434939763784110e-44 }, 
    {  3.268798837811031473e-30, -3.208151322045558000e-47 },
    { -3.757240043460955716e-33, -3.102604726588237581e-50 } 
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
inline void fast2Sum(double a, double b, double &s, double &t) {
    s = a + b;
    double z = s - a;
    t = b - z;
}

// точное сложение с компенсацией ошибки
inline void is2Sum(double a, double b, double &s, double &t) {
    s = a + b;
    double a_ = s - b;
    double b_ = s - a_;
    double da = a - a_;
    double db = b - b_;
    t = da + db;
}

// разделение double на две части методом Велткампа
inline void split(double x, double &x_hi, double &x_lo) {
    const double C = 134217729.0; // 2^27 + 1
    double gamma = C * x;
    double delta = x - gamma;
    x_hi = gamma + delta;
    x_lo = x - x_hi;
}

// точное произведение двух double
inline void dekker(double x, double y, double &p_hi, double &p_lo) {
    double x_hi, x_lo, y_hi, y_lo;
    split(x, x_hi, x_lo);
    split(y, y_hi, y_lo);
    p_hi = x * y;
    double t1 = -p_hi + x_hi * y_hi;
    double t2 = t1 + x_hi * y_lo;
    double t3 = t2 + x_lo * y_hi;
    p_lo = t3 + x_lo * y_lo;
}

// приведение трёх чисел к нормализованной паре head+tail
DDouble normalize3(double a0, double a1, double a2) {
    double s, t2;
    fast2Sum(a1, a2, s, t2);
    double t0, t1;
    fast2Sum(a0, s, t0, t1);
    double b0, b1;
    fast2Sum(t0, t1, b0, s);
    fast2Sum(b0, s + t2, b0, b1);
    return DDouble(b0, b1);
}

// Арифметика DDouble 
// сложение двух double-double чисел
DDouble add_dd(const DDouble &x, const DDouble &y) {
    double sh, sl, th, tl, c, vh, vl, wh, zh, zl;
    is2Sum(x.hi, y.hi, sh, sl);   // void, но sh, sl заполнятся
    is2Sum(x.lo, y.lo, th, tl);
    c = sl + th;
    fast2Sum(sh, c, vh, vl);
    wh = tl + vl;
    fast2Sum(vh, wh, zh, zl);
    return DDouble(zh, zl);
}

// умножение двух double-double чисел
DDouble mul_dd(const DDouble &x, const DDouble &y) {
    double ph, pl;
    dekker(x.hi, y.hi, ph, pl);
    pl += x.hi * y.lo + x.lo * y.hi;
    double zh, zl;
    fast2Sum(ph, pl, zh, zl);
    return DDouble(zh, zl);
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
    return quick_two_sum(q1, q2);
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
    int k_int = static_cast<int>(k);
        
    
    double k1 = std::floor(k_int * 1e-9) * 1e9;
    double k2 = k_int - k1;
    
    DDouble r_val = x;
    double parts_k[2] = {k1, k2};
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 5; ++j) {
            DDouble term = mul_dd(DDouble(parts_k[i]), DDouble(pi_div_two[j]));
            r_val = add_dd(r_val, neg_dd(term));
        }
    }
    
    
    const double PI4 = 0.7853981633974483;
    int quad = k_int & 3;
    DDouble half_pi(pi_div_two[0], pi_div_two[1]);
    while (r_val.hi > PI4) {
        r_val = add_dd(r_val, neg_dd(half_pi));
        quad = (quad + 1) % 4;
    }
    while (r_val.hi < -PI4) {
        r_val = add_dd(r_val, half_pi);
        quad = (quad + 3) % 4;
    }
    
    r = r_val;
    return quad;
}

// основная функция синуса с учётом знака и квадранта
DDouble dd_sin(const DDouble &x) {
    // Запоминаем знак исходного аргумента
    bool neg = (x.hi < 0.0) || (x.hi == 0.0 && x.lo < 0.0);
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

// Вспомогательный ввод
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