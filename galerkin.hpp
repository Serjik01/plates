#pragma once
// ---------------------------------------------------------------------------
//  Метод Бубнова-Галёркина.
//
//  w_N = sum c_i * phi_i ,  невязка R = L[w_N] - q/D1 ортогональна всем phi_j:
//      K_ji * c_i = F_j ,   K_ji = int int L[phi_i]*phi_j dA
//                           F_j  = int int (q/D1)*phi_j dA
//
//  Матрица собирается ЧИСЛЕННО по общей формуле -- диагональность в случае
//  тригонометрического базиса не закладывается, а получается как результат.
// ---------------------------------------------------------------------------
#include <cmath>
#include <vector>

#include "basis.hpp"
#include "linalg.hpp"
#include "plate.hpp"
#include "quad.hpp"

struct GalerkinSolution {
    const TensorBasis*  basis = nullptr;
    std::vector<double> c;              // коэффициенты разложения
    double              pivotRatio = 0; // оценка обусловленности
    double              offDiagRatio = 0;   // max|K_ij|, i!=j, к max|K_ii|
    int                 N = 0;

    double deriv(int dx, int dy, double x, double y) const {
        double acc = 0.0;
        for (int i = 0; i < N; ++i) acc += c[i] * basis->val(i, dx, dy, x, y);
        return acc;
    }
    double w(double x, double y) const { return deriv(0, 0, x, y); }
};

inline GalerkinSolution solveGalerkin(const Plate& p, LoadType lt,
                                      const TensorBasis& basis, int nq = 80)
{
    const double d1 = p.delta1(), d2 = p.delta2(), D1 = p.D1();
    const int    N  = basis.size();

    GaussLegendre gx(nq, p.a), gy(nq, p.b);

    // --- предвычисление phi_i и L[phi_i] в узлах квадратуры ---
    const int G = nq * nq;
    std::vector<double> phi(std::size_t(N) * G), Lphi(std::size_t(N) * G);
    std::vector<double> wq(G), qv(G);

    for (int i = 0; i < nq; ++i)
        for (int j = 0; j < nq; ++j) {
            const int g = i * nq + j;
            wq[g] = gx.w[i] * gy.w[j];
            qv[g] = loadValue(p, lt, gx.x[i], gy.x[j]) / D1;
        }

    for (int k = 0; k < N; ++k)
        for (int i = 0; i < nq; ++i)
            for (int j = 0; j < nq; ++j) {
                const int g = i * nq + j;
                phi [std::size_t(k) * G + g] = basis.val (k, 0, 0, gx.x[i], gy.x[j]);
                Lphi[std::size_t(k) * G + g] = basis.oper(k, d1, d2, gx.x[i], gy.x[j]);
            }

    // --- сборка ---
    Matrix K(N, Vector(N, 0.0));
    Vector F(N, 0.0);

    for (int j = 0; j < N; ++j) {
        const double* pj = &phi[std::size_t(j) * G];
        double f = 0.0;
        for (int g = 0; g < G; ++g) f += wq[g] * qv[g] * pj[g];
        F[j] = f;
        for (int i = 0; i < N; ++i) {
            const double* li = &Lphi[std::size_t(i) * G];
            double s = 0.0;
            for (int g = 0; g < G; ++g) s += wq[g] * li[g] * pj[g];
            K[j][i] = s;
        }
    }

    // --- диагностика структуры матрицы ---
    double dmax = 0.0, omax = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            (i == j ? dmax : omax) = std::max(i == j ? dmax : omax, std::fabs(K[i][j]));

    GalerkinSolution sol;
    sol.basis        = &basis;
    sol.N            = N;
    sol.offDiagRatio = (dmax > 0.0 ? omax / dmax : 0.0);
    sol.c            = gaussSolve<double>(K, F, &sol.pivotRatio);
    return sol;
}

// ---------------------------------------------------------------------------
//  Демонстрация вырождения "наивного" базиса phi = x(a-x)*y(b-y).
//  У него phi,xxxx = phi,yyyy = 0, поэтому из оператора выпадают ОБЕ изгибные
//  жёсткости и остаётся только смешанный член -- результат нефизичен.
// ---------------------------------------------------------------------------
inline double naiveBasisCenterDeflection(const Plate& p, LoadType lt, int nq = 60)
{
    const double d1 = p.delta1(), D1 = p.D1();
    GaussLegendre gx(nq, p.a), gy(nq, p.b);

    auto phi = [&](double x, double y) { return x * (p.a - x) * y * (p.b - y); };
    // L[phi] = 0 + 2*d1*phi,xxyy + 0 = 2*d1*(-2)*(-2) = 8*d1
    const double Lphi = 8.0 * d1;

    double K = 0.0, F = 0.0;
    for (int i = 0; i < nq; ++i)
        for (int j = 0; j < nq; ++j) {
            const double x = gx.x[i], y = gy.x[j], w = gx.w[i] * gy.w[j];
            const double ph = phi(x, y);
            K += w * Lphi * ph;
            F += w * (loadValue(p, lt, x, y) / D1) * ph;
        }
    const double c = F / K;
    return c * phi(0.5 * p.a, 0.5 * p.b);
}
