import math
from decimal import Decimal, getcontext

getcontext().prec = 100

pi2 = math.pi / 2
pio2_h0 = 1.5707963267948966  # то же самое, что float(pi2)
pio2_h1 = pi2 - pio2_h0

print(f"PIO2_H0 = {pio2_h0:.16f}")
print(f"PIO2_H1 = {pio2_h1:.16e}")
print()

two_over_pi = 2 / math.pi
two_pi_hi = float(two_over_pi)
two_pi_lo = float(two_over_pi - two_pi_hi)

print(f"{two_pi_hi:.16f}")
print(f"{two_pi_lo:.16e}")

print()

for n in range(0, 9):
    coeff = (-1)**n / math.factorial(2*n + 1)
    print(f"{coeff:.17e}")

print()

for n in range(0, 9):
    coeff = (-1)**n / math.factorial(2*n)
    print(f"{coeff:.17e}")