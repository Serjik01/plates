#pragma once
// ---------------------------------------------------------------------------
//  Метод Навье.
//  Прогиб раскладывается в двойной ряд по собственным функциям оператора:
//      w = sum_{m,n} A_mn * sin(m*pi*x/a) * sin(n*pi*y/b)
//  Коэффициенты нагрузки q_mn считаются квадратурой (годится для любой q),
//  коэффициенты прогиба -- по замкнутой формуле.
// ---------------------------------------------------------------------------
#include <cmath>
#include <vector>

#include "plate.hpp"
#include "quad.hpp"

struct NavierSolution {
    int                 K = 0;      // максимальный индекс по каждой координате
    double              a = 0, b = 0;
    std::vector<double> A;          // A[(m-1)*K + (n-1)]

    // производная порядка (dx, dy) от прогиба
    double deriv(int dx, int dy, double x, double y) const {
        const double PI = 3.14159265358979323846;
        double acc = 0.0;
        for (int m = 1; m <= K; ++m) {
            const double lx = m * PI / a;
            double fx;
            switch (dx) {
                case 0: fx =  std::sin(lx * x);           break;
                case 1: fx =  lx * std::cos(lx * x);      break;
                case 2: fx = -lx * lx * std::sin(lx * x); break;
                default: fx = 0.0;
            }
            if (fx == 0.0) continue;
            for (int n = 1; n <= K; ++n) {
                const double c = A[(m - 1) * K + (n - 1)];
                if (c == 0.0) continue;
                const double ly = n * PI / b;
                double fy;
                switch (dy) {
                    case 0: fy =  std::sin(ly * y);           break;
                    case 1: fy =  ly * std::cos(ly * y);      break;
                    case 2: fy = -ly * ly * std::sin(ly * y); break;
                    default: fy = 0.0;
                }
                acc += c * fx * fy;
            }
        }
        return acc;
    }

    double w(double x, double y) const { return deriv(0, 0, x, y); }

    // усечение до K' <= K членов (для исследования сходимости без пересчёта)
    NavierSolution truncated(int Kp) const {
        NavierSolution s;
        s.K = Kp; s.a = a; s.b = b;
        s.A.assign(std::size_t(Kp) * Kp, 0.0);
        for (int m = 1; m <= Kp; ++m)
            for (int n = 1; n <= Kp; ++n)
                s.A[(m - 1) * Kp + (n - 1)] = A[(m - 1) * K + (n - 1)];
        return s;
    }
};

inline NavierSolution solveNavier(const Plate& p, LoadType lt, int K, int nq = 0)
{
    const double PI = 3.14159265358979323846;
    const double d1 = p.delta1(), d2 = p.delta2(), D1 = p.D1();
    if (nq <= 0) nq = std::max(60, 5 * K);   // ~5 узлов на полуволну старшей гармоники

    NavierSolution sol;
    sol.K = K; sol.a = p.a; sol.b = p.b;
    sol.A.assign(std::size_t(K) * K, 0.0);

    GaussLegendre gx(nq, p.a), gy(nq, p.b);

    // таблицы синусов
    std::vector<double> SX(std::size_t(K) * nq), SY(std::size_t(K) * nq);
    for (int m = 1; m <= K; ++m)
        for (int i = 0; i < nq; ++i)
            SX[std::size_t(m - 1) * nq + i] = std::sin(m * PI * gx.x[i] / p.a);
    for (int n = 1; n <= K; ++n)
        for (int j = 0; j < nq; ++j)
            SY[std::size_t(n - 1) * nq + j] = std::sin(n * PI * gy.x[j] / p.b);

    // значения нагрузки в узлах
    std::vector<double> Q(std::size_t(nq) * nq);
    for (int i = 0; i < nq; ++i)
        for (int j = 0; j < nq; ++j)
            Q[std::size_t(i) * nq + j] = loadValue(p, lt, gx.x[i], gy.x[j]);

    // Sx[m][j] = sum_i w_i * q(x_i, y_j) * sin(m pi x_i / a)
    std::vector<double> Sx(std::size_t(K) * nq, 0.0);
    for (int m = 0; m < K; ++m)
        for (int i = 0; i < nq; ++i) {
            const double f = gx.w[i] * SX[std::size_t(m) * nq + i];
            if (f == 0.0) continue;
            for (int j = 0; j < nq; ++j)
                Sx[std::size_t(m) * nq + j] += f * Q[std::size_t(i) * nq + j];
        }

    for (int m = 1; m <= K; ++m) {
        const double lx = m * PI / p.a;
        for (int n = 1; n <= K; ++n) {
            double q_mn = 0.0;
            for (int j = 0; j < nq; ++j)
                q_mn += gy.w[j] * Sx[std::size_t(m - 1) * nq + j]
                                * SY[std::size_t(n - 1) * nq + j];
            q_mn *= 4.0 / (p.a * p.b);

            const double ly = n * PI / p.b;
            const double den = std::pow(lx, 4)
                             + 2.0 * d1 * lx * lx * ly * ly
                             + d2 * std::pow(ly, 4);
            sol.A[(m - 1) * K + (n - 1)] = q_mn / (D1 * den);
        }
    }
    return sol;
}
