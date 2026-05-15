#ifndef RK38_H
#define RK38_H

struct RK38 {
    static constexpr int stages = 4;
    static constexpr double a[stages][stages-1] = {
        {0.0, 0.0, 0.0},
        {1.0/3.0, 0.0, 0.0},
        {-1.0/3.0, 1.0, 0.0},
        {1.0, -1.0, 1.0}
    };
    static constexpr double b[stages] = {1.0/8.0, 3.0/8.0, 3.0/8.0, 1.0/8.0};
    static constexpr double c[stages] = {0.0, 1.0/3.0, 2.0/3.0, 1.0};
};

#endif