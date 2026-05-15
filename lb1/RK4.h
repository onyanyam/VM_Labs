#ifndef RK4_H
#define RK4_H

struct RK4 {
    static constexpr int stages = 4;
    // a[i][j] для i = 0..3, j = 0..i-1 (индексы с 0)
    static constexpr double a[stages][stages-1] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 0.0, 1.0}
    };
    static constexpr double b[stages] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    static constexpr double c[stages] = {0.0, 0.5, 0.5, 1.0};
};

#endif