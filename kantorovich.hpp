#pragma once
// ---------------------------------------------------------------------------
//  Метод Канторовича-Власова.
//
//  Нагрузка синусоидальна по y (нечётные варианты):  q = q0*sin(pi*y/b).
//  Решение ищется в виде   w(x,y) = f(x) * sin(s*y),  s = pi/b.
//  Множитель sin(s*y) тождественно удовлетворяет шарнирному опиранию на
//  кромках y = 0 и y = b. Подстановка в разрешающее уравнение даёт ОДУ:
//
//      f'''' - 2*delta1*s^2*f'' + delta2*s^4*f = q0/D1
//
//  Характеристическое уравнение  r^4 - 2*d1*s^2*r^2 + d2*s^4 = 0 решается
//  в комплексных числах -- одна ветвь кода покрывает все типы корней.
//  Экспоненты масштабируются (сдвиг аргумента), чтобы |exp| <= 1.
// ---------------------------------------------------------------------------
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "linalg.hpp"
#include "plate.hpp"

using cplx = std::complex<double>;

enum class EdgeBC { Clamped, Simple, Free };

inline const char* bcName(EdgeBC bc) {
    switch (bc) {
        case EdgeBC::Clamped: return "жёсткая заделка";
        case EdgeBC::Simple:  return "шарнирное опирание";
        case EdgeBC::Free:    return "свободный край";
    }
    return "?";
}

struct KantorovichSolution {
    double              a = 0, b = 0, s = 0;
    double              fp = 0;              // частное решение (константа)
    std::array<cplx, 4> r{};                 // корни характеристического уравнения
    std::array<double, 4> shift{};           // сдвиги аргумента экспонент
    std::array<cplx, 4> C{};                 // постоянные интегрирования

    // производная f порядка d
    double f(int d, double x) const {
        cplx acc(0.0, 0.0);
        for (int k = 0; k < 4; ++k)
            acc += C[k] * std::pow(r[k], d) * std::exp(r[k] * cplx(x - shift[k], 0.0));
        return acc.real() + (d == 0 ? fp : 0.0);
    }

    double deriv(int dx, int dy, double x, double y) const {
        const double sy = std::sin(s * y), cy = std::cos(s * y);
        double gy;
        switch (dy) {
            case 0: gy =  sy;         break;
            case 1: gy =  s * cy;     break;
            case 2: gy = -s * s * sy; break;
            default: gy = 0.0;
        }
        return f(dx, x) * gy;
    }
    double w(double x, double y) const { return deriv(0, 0, x, y); }
};

// bc0 -- условие на кромке x = 0, bcA -- на кромке x = a
inline KantorovichSolution solveKantorovich(const Plate& p, EdgeBC bc0, EdgeBC bcA)
{
    const double PI = 3.14159265358979323846;
    const double d1 = p.delta1(), d2 = p.delta2(), D1 = p.D1();
    const double s  = PI / p.b;
    const double nu2 = p.nu2, kap = p.kappa();

    KantorovichSolution sol;
    sol.a = p.a; sol.b = p.b; sol.s = s;
    sol.fp = p.q0 / (D1 * d2 * std::pow(s, 4));

    // --- корни: u = r^2/s^2 = d1 +- sqrt(d1^2 - d2) ---
    const cplx disc = std::sqrt(cplx(d1 * d1 - d2, 0.0));
    const cplx u1 = cplx(d1, 0.0) + disc;
    const cplx u2 = cplx(d1, 0.0) - disc;
    if (std::abs(u1 - u2) < 1e-12)
        throw std::runtime_error("кратные корни характеристического уравнения: "
                                 "нужен базис с множителем x");

    const cplx q1 = std::sqrt(u1) * s;
    const cplx q2 = std::sqrt(u2) * s;
    sol.r = { q1, -q1, q2, -q2 };

    for (int k = 0; k < 4; ++k)
        sol.shift[k] = (sol.r[k].real() > 0.0 ? p.a : 0.0);

    // --- граничные условия как линейные функционалы от (f, f', f'', f''') ---
    auto rows = [&](EdgeBC bc) {
        std::array<std::array<double, 4>, 2> R{};
        switch (bc) {
            case EdgeBC::Clamped:                       // f = 0, f' = 0
                R[0] = {1, 0, 0, 0};
                R[1] = {0, 1, 0, 0};
                break;
            case EdgeBC::Simple:                        // f = 0, M_x = 0
                R[0] = {1, 0, 0, 0};
                R[1] = {-nu2 * s * s, 0, 1, 0};
                break;
            case EdgeBC::Free:                          // M_x = 0, V_x = 0
                R[0] = {-nu2 * s * s, 0, 1, 0};
                R[1] = {0, -kap * s * s, 0, 1};
                break;
        }
        return R;
    };

    const auto R0 = rows(bc0);
    const auto RA = rows(bcA);

    std::vector<std::vector<cplx>> A(4, std::vector<cplx>(4));
    std::vector<cplx>              rhs(4);

    auto fill = [&](int row, const std::array<double, 4>& coef, double xe) {
        for (int k = 0; k < 4; ++k) {
            const cplx g = std::exp(sol.r[k] * cplx(xe - sol.shift[k], 0.0));
            cplx acc(0.0, 0.0);
            for (int d = 0; d < 4; ++d)
                if (coef[d] != 0.0) acc += coef[d] * std::pow(sol.r[k], d) * g;
            A[row][k] = acc;
        }
        // частное решение -- константа, вклад даёт только coef[0]
        rhs[row] = cplx(-coef[0] * sol.fp, 0.0);
    };

    fill(0, R0[0], 0.0);
    fill(1, R0[1], 0.0);
    fill(2, RA[0], p.a);
    fill(3, RA[1], p.a);

    const auto C = gaussSolve<cplx>(A, rhs);
    for (int k = 0; k < 4; ++k) sol.C[k] = C[k];
    return sol;
}
