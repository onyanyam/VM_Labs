#ifndef RALSTON3_H
#define RALSTON3_H

struct Ralston3 {
    static constexpr int stages = 3;
    static constexpr double a[stages][stages-1] = {
        {0.0, 0.0},
        {0.5, 0.0},
        {0.0, 0.75}
    };
    static constexpr double b[stages] = {2.0/9.0, 1.0/3.0, 4.0/9.0};
    static constexpr double c[stages] = {0.0, 0.5, 0.75};
};

#endif