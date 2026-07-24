#pragma once
// ---------------------------------------------------------------------------
//  Геометрия, материал и нагрузки тонкой ортотропной пластины
//  Разрешающее уравнение (нормировано на D1):
//      w,xxxx + 2*delta1*w,xxyy + delta2*w,yyyy = q(x,y) / D1
// ---------------------------------------------------------------------------
#include <cmath>
#include <stdexcept>
#include <string>

struct Plate {
    // --- геометрия ---
    double a  = 0.75;      // размер вдоль x, м
    double b  = 0.50;      // размер вдоль y, м
    double h  = 0.01;      // толщина, м

    // --- материал (ортотропный) ---
    double E1  = 1.2e11;   // Па, модуль вдоль x
    double E2  = 0.6e11;   // Па, модуль вдоль y
    double nu1 = 0.071;    // коэффициент Пуассона
    double nu2 = 0.036;
    double G   = 0.7e10;   // Па, модуль сдвига в плоскости пластины

    // --- нагрузка ---
    double q0 = 1.0e4;     // Па (10 кПа)

    // --- цилиндрические жёсткости ---
    double D1()  const { return E1 * h * h * h / (12.0 * (1.0 - nu1 * nu2)); }
    double D2()  const { return E2 * h * h * h / (12.0 * (1.0 - nu1 * nu2)); }
    double D12() const { return G  * h * h * h / 12.0; }

    // delta1 = nu2 + 2*G*(1-nu1*nu2)/E1  == nu2 + 2*D12/D1
    double delta1() const { return nu2 + 2.0 * D12() / D1(); }

    // delta2 = D2/D1 = E2/E1  (условие взаимности даёт также nu2/nu1,
    // но в исходных данных nu2 округлено, поэтому берём отношение модулей)
    double delta2() const { return D2() / D1(); }

    // множитель в условии свободного края для обобщённой силы Кирхгофа:
    //   V_x = -D1*w,xxx + (nu2 + 4*D12/D1)*D1*w,xyy  ,   nu2 + 4*D12/D1 = 2*delta1 - nu2
    double kappa() const { return 2.0 * delta1() - nu2; }

    void check() const {
        if (h <= 0 || a <= 0 || b <= 0) throw std::runtime_error("bad geometry");
        if (std::fabs(nu1 * E2 - nu2 * E1) / (nu1 * E2) > 0.05)
            throw std::runtime_error("нарушено условие взаимности E1*nu2 = E2*nu1");
    }
};

// --- типы поперечной нагрузки ---
enum class LoadType {
    Uniform,   // q = q0                          (пункт а)
    SinXY,     // q = q0*sin(pi x/a)*sin(pi y/b)  (пункт б)
    SinY       // q = q0*sin(pi y/b)              (пункт в, нечётный вариант)
};

inline const char* loadName(LoadType t) {
    switch (t) {
        case LoadType::Uniform: return "uniform";
        case LoadType::SinXY:   return "sin_xy";
        case LoadType::SinY:    return "sin_y";
    }
    return "?";
}

inline double loadValue(const Plate& p, LoadType t, double x, double y) {
    const double PI = 3.14159265358979323846;
    switch (t) {
        case LoadType::Uniform: return p.q0;
        case LoadType::SinXY:   return p.q0 * std::sin(PI * x / p.a) * std::sin(PI * y / p.b);
        case LoadType::SinY:    return p.q0 * std::sin(PI * y / p.b);
    }
    return 0.0;
}
