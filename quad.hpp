#pragma once
// ---------------------------------------------------------------------------
//  Квадратура Гаусса-Лежандра на отрезке [0, L].
//  Узлы -- корни полинома Лежандра P_n, находятся методом Ньютона.
//  n узлов интегрируют точно любой полином степени <= 2n-1.
// ---------------------------------------------------------------------------
#include <cmath>
#include <vector>

struct GaussLegendre {
    std::vector<double> x;   // узлы на [0, L]
    std::vector<double> w;   // веса

    GaussLegendre(int n, double L) : x(n), w(n) {
        const double PI = 3.14159265358979323846;
        for (int i = 0; i < n; ++i) {
            // начальное приближение (асимптотика Чебышёва)
            double t = std::cos(PI * (i + 0.75) / (n + 0.5));
            for (int it = 0; it < 100; ++it) {
                // рекуррентно считаем P_n(t) и P_{n-1}(t)
                double p0 = 1.0, p1 = 0.0;
                for (int k = 0; k < n; ++k) {
                    const double p2 = p1;
                    p1 = p0;
                    p0 = ((2.0 * k + 1.0) * t * p1 - k * p2) / (k + 1.0);
                }
                const double dp = n * (t * p0 - p1) / (t * t - 1.0);  // P_n'(t)
                const double dt = -p0 / dp;
                t += dt;
                if (std::fabs(dt) < 1e-15) break;
            }
            // пересчёт производной в найденном корне для веса
            double p0 = 1.0, p1 = 0.0;
            for (int k = 0; k < n; ++k) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * k + 1.0) * t * p1 - k * p2) / (k + 1.0);
            }
            const double dp = n * (t * p0 - p1) / (t * t - 1.0);

            // отображение [-1,1] -> [0,L]
            x[i] = 0.5 * L * (t + 1.0);
            w[i] = L / ((1.0 - t * t) * dp * dp);
        }
    }
};
