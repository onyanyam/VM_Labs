import mpmath
mpmath.mp.dps = 100

def split_to_hi_lo(value):
    hi = float(value)
    lo = float(value - mpmath.mpf(hi))
    return hi, lo

def print_dd_constant(name, value):
    hi, lo = split_to_hi_lo(value)
    print(f"static const DDouble {name}({hi:.17e}, {lo:.17e});")

def factorial_reciprocal(n):
    return mpmath.mpf(1) / mpmath.factorial(n)

def generate_sin_coeffs(max_n=41):
    coeffs = []
    for k in range(1, max_n + 1, 2):
        exact = factorial_reciprocal(k)
        hi, lo = split_to_hi_lo(exact)
        if (k // 2) % 2 == 1:
            hi = -hi
            lo = -lo
        coeffs.append((hi, lo))
    return coeffs

def generate_cos_coeffs(max_n=42):
    coeffs = []
    for k in range(0, max_n + 1, 2):
        if k == 0:
            exact = mpmath.mpf(1)
        else:
            exact = factorial_reciprocal(k)
        hi, lo = split_to_hi_lo(exact)
        if (k // 2) % 2 == 1:
            hi = -hi
            lo = -lo
        coeffs.append((hi, lo))
    return coeffs

def print_coeffs(coeffs, name):
    print(f"static const DDouble {name}[] = {{")
    for hi, lo in coeffs:
        print(f"    {{ {hi:.17e}, {lo:.17e} }},")
    print("};")


if __name__ == "__main__":
    sin_coeffs = generate_sin_coeffs(41)
    print("// Коэффициенты синуса (1, -1/3!, +1/5!, ..., до 1/41! )")
    print_coeffs(sin_coeffs, "sin_coeffs")
    print()

    cos_coeffs = generate_cos_coeffs(40)
    print("// Коэффициенты косинуса (1, -1/2!, +1/4!, ..., до 1/40! )")
    print_coeffs(cos_coeffs, "cos_coeffs")
    print()

    print("// Константы π/2, π/4")
    pi_2 = mpmath.pi / 2
    print_dd_constant("DD_PI_2", pi_2)

    pi_4 = mpmath.pi / 4
    print_dd_constant("DD_PI_4", pi_4)
    print()