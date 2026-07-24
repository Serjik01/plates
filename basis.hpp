#pragma once
// ---------------------------------------------------------------------------
//  Координатные (базисные) функции метода Бубнова-Галёркина.
//
//  Оба используемых базиса -- тензорные произведения одномерных наборов,
//  каждый из которых точно удовлетворяет условиям шарнирного опирания
//  на своих концах:   X(0) = X(L) = 0,   X''(0) = X''(L) = 0.
//
//  Все производные (до 4-го порядка) вычисляются АНАЛИТИЧЕСКИ -- никаких
//  конечных разностей, поэтому невязка ограничена только машинной точностью.
// ---------------------------------------------------------------------------
#include <cmath>
#include <string>
#include <vector>

// ============================ базовый интерфейс ============================
struct Set1D {
    virtual ~Set1D() = default;
    virtual int         size() const = 0;
    virtual std::string name() const = 0;
    // производная порядка `order` (0..4) k-й функции в точке x
    virtual double d(int k, int order, double x) const = 0;
};

// ====================== 1. тригонометрический базис ========================
//   X_k(x) = sin(k*pi*x/L),  k = 1..n
struct TrigSet : Set1D {
    double L;
    int    n;
    TrigSet(int n_, double L_) : L(L_), n(n_) {}

    int         size() const override { return n; }
    std::string name() const override { return "trig"; }

    double d(int k, int order, double x) const override {
        const double PI = 3.14159265358979323846;
        const double lam = (k + 1) * PI / L;          // k нумеруется с 0
        const double s   = std::sin(lam * x);
        const double c   = std::cos(lam * x);
        switch (order) {
            case 0: return  s;
            case 1: return  lam * c;
            case 2: return -lam * lam * s;
            case 3: return -lam * lam * lam * c;
            case 4: return  lam * lam * lam * lam * s;
        }
        return 0.0;
    }
};

// ========================= 2. полиномиальный базис =========================
//   Балочные полиномы: X_k -- решение X'''' = xi^(k-1) при X = X'' = 0 на концах.
//   В безразмерной координате xi = x/L:
//       X_k(xi) = xi^(k+3) - C*xi^3 + (C-1)*xi ,   C = (k+3)(k+2)/6
//   Первый из них (k=1):  xi^4 - 2*xi^3 + xi  -- классическая балочная функция.
//
//   Каждая функция нормируется на свой максимум по модулю: это не меняет
//   натянутое подпространство, но заметно улучшает обусловленность матрицы.
struct PolySet : Set1D {
    double                            L;
    int                               n;
    std::vector<std::vector<double>>  c;   // c[k][p] -- коэффициент при xi^p

    // скалярное произведение двух полиномов на [0,1]:  int xi^p * xi^q = 1/(p+q+1)
    static double dot(const std::vector<double>& u, const std::vector<double>& v) {
        double s = 0.0;
        for (std::size_t p = 0; p < u.size(); ++p) {
            if (u[p] == 0.0) continue;
            for (std::size_t q = 0; q < v.size(); ++q)
                if (v[q] != 0.0) s += u[p] * v[q] / double(p + q + 1);
        }
        return s;
    }

    PolySet(int n_, double L_) : L(L_), n(n_), c(n_) {
        // --- сырые балочные полиномы ---
        for (int k = 1; k <= n; ++k) {
            const double C = double(k + 3) * double(k + 2) / 6.0;
            std::vector<double> co(std::size_t(n) + 4, 0.0);   // общая длина -- удобно для МГШ
            co[k + 3] = 1.0;
            co[3]    += -C;
            co[1]    += (C - 1.0);
            c[k - 1] = std::move(co);
        }
        // --- модифицированный процесс Грама-Шмидта в L2[0,1] ---
        // Линейные комбинации функций, удовлетворяющих однородным граничным
        // условиям, удовлетворяют им же, поэтому базис остаётся допустимым,
        // а обусловленность матрицы Галёркина улучшается на порядки.
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < k; ++j) {
                const double pr = dot(c[k], c[j]);
                for (std::size_t p = 0; p < c[k].size(); ++p) c[k][p] -= pr * c[j][p];
            }
            const double nrm = std::sqrt(dot(c[k], c[k]));
            if (nrm > 0.0) for (double& v : c[k]) v /= nrm;
        }
    }

    int         size() const override { return n; }
    std::string name() const override { return "poly"; }

    double d(int k, int order, double x) const override {
        const double xi  = x / L;
        const auto&  co  = c[k];
        double       acc = 0.0;
        for (std::size_t p = order; p < co.size(); ++p) {
            if (co[p] == 0.0) continue;
            double fact = 1.0;                       // p!/(p-order)!
            for (int t = 0; t < order; ++t) fact *= double(p - t);
            acc += co[p] * fact * std::pow(xi, double(p - order));
        }
        return acc / std::pow(L, double(order));
    }
};

// ==================== тензорное произведение наборов =======================
struct TensorBasis {
    const Set1D& X;
    const Set1D& Y;

    TensorBasis(const Set1D& x, const Set1D& y) : X(x), Y(y) {}

    int size() const { return X.size() * Y.size(); }
    int ix(int i) const { return i / Y.size(); }
    int iy(int i) const { return i % Y.size(); }

    // d^(dx+dy) phi_i / dx^dx dy^dy
    double val(int i, int dx, int dy, double x, double y) const {
        return X.d(ix(i), dx, x) * Y.d(iy(i), dy, y);
    }

    // оператор задачи:  L[phi] = phi,xxxx + 2*d1*phi,xxyy + d2*phi,yyyy
    double oper(int i, double d1, double d2, double x, double y) const {
        return val(i, 4, 0, x, y)
             + 2.0 * d1 * val(i, 2, 2, x, y)
             + d2 * val(i, 0, 4, x, y);
    }

    std::string name() const { return X.name() + "x" + Y.name(); }
};
