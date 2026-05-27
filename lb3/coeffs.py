import numpy.polynomial.legendre as leg

def generate_cpp_table(n_points, precision=16):
    nodes, weights = leg.leggauss(n_points)
    print(f"static const int GAUSS_N = {n_points};")
    print(f"static const double GAUSS_X[{n_points}] = {{")
    for i, x in enumerate(nodes):
        comma = "," if i < n_points - 1 else ""
        print(f"    {x:.{precision}f}{comma}")
    print(f"}};")
    print()
    print(f"static const double GAUSS_W[{n_points}] = {{")
    for i, w in enumerate(weights):
        comma = "," if i < n_points - 1 else ""
        print(f"    {w:.{precision}f}{comma}")
    print(f"}};")


generate_cpp_table(100)