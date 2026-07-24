#pragma once
// ---------------------------------------------------------------------------
//  Плотная линейная алгебра: метод Гаусса с частичным выбором главного
//  элемента. Шаблон работает и для double, и для std::complex<double>.
// ---------------------------------------------------------------------------
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

template <class T> inline double magnitude(const T& v)              { return std::fabs(v); }
template <>        inline double magnitude<std::complex<double>>(const std::complex<double>& v) { return std::abs(v); }

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// Решает A*x = rhs на месте. Возвращает x.
// Если pivotRatio != nullptr, туда кладётся отношение max|pivot| / min|pivot| --
// грубая, но дешёвая оценка обусловленности матрицы.
template <class T>
std::vector<T> gaussSolve(std::vector<std::vector<T>> A,
                          std::vector<T>              rhs,
                          double*                     pivotRatio = nullptr)
{
    const std::size_t n = rhs.size();
    if (A.size() != n) throw std::runtime_error("gaussSolve: несогласованные размеры");

    double pmax = 0.0, pmin = 0.0;
    bool   first = true;

    for (std::size_t k = 0; k < n; ++k) {
        // --- выбор главного элемента ---
        std::size_t piv = k;
        double      best = magnitude(A[k][k]);
        for (std::size_t i = k + 1; i < n; ++i) {
            const double m = magnitude(A[i][k]);
            if (m > best) { best = m; piv = i; }
        }
        if (best == 0.0) throw std::runtime_error("gaussSolve: вырожденная матрица");
        if (piv != k) { std::swap(A[piv], A[k]); std::swap(rhs[piv], rhs[k]); }

        if (first) { pmax = pmin = best; first = false; }
        else { if (best > pmax) pmax = best; if (best < pmin) pmin = best; }

        // --- исключение ---
        for (std::size_t i = k + 1; i < n; ++i) {
            const T f = A[i][k] / A[k][k];
            if (magnitude(f) == 0.0) continue;
            for (std::size_t j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            rhs[i] -= f * rhs[k];
        }
    }

    if (pivotRatio) *pivotRatio = (pmin > 0.0 ? pmax / pmin : 0.0);

    // --- обратный ход ---
    std::vector<T> x(n);
    for (std::size_t ii = n; ii-- > 0;) {
        T s = rhs[ii];
        for (std::size_t j = ii + 1; j < n; ++j) s -= A[ii][j] * x[j];
        x[ii] = s / A[ii][ii];
    }
    return x;
}
