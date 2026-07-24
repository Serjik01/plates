// ---------------------------------------------------------------------------
//  Изгиб тонкой ортотропной прямоугольной пластины. Вариант 3.
//
//    а) шарнирное опирание по контуру, равномерная нагрузка
//       -- метод Навье и метод Бубнова-Галёркина (два базиса), сходимость
//    б) то же при синусоидальной по x и y нагрузке
//    в) нагрузка q0*sin(pi*y/b); кромки y=0,b шарнирно оперты,
//       x=0 -- жёсткая заделка, x=a -- свободный край
//       -- метод Канторовича-Власова
//
//  Сборка:  g++ -std=c++17 -O2 main.cpp -o plate
// ---------------------------------------------------------------------------
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "basis.hpp"
#include "galerkin.hpp"
#include "kantorovich.hpp"
#include "navier.hpp"
#include "output.hpp"
#include "plate.hpp"

namespace fs = std::filesystem;

static const std::string OUT = "out";
static const int GX = 121, GY = 81;      // сетка вывода полей

static void rule(const char* title) {
    std::cout << "\n=====================================================\n"
              << "  " << title << "\n"
              << "=====================================================\n";
}

// центральная точка пластины
template <class S> static double wCenter(const Plate& p, const S& s) {
    return s.deriv(0, 0, 0.5 * p.a, 0.5 * p.b);
}

// ---------------------------------------------------------------------------
static void printParams(const Plate& p)
{
    rule("Исходные данные и жёсткостные характеристики");
    std::cout << std::fixed;
    std::printf("  a = %.3f м,  b = %.3f м,  h = %.4f м,  q0 = %.1f Па\n",
                p.a, p.b, p.h, p.q0);
    std::printf("  E1 = %.3e Па,  E2 = %.3e Па,  G = %.3e Па\n", p.E1, p.E2, p.G);
    std::printf("  nu1 = %.4f,  nu2 = %.4f\n", p.nu1, p.nu2);
    std::printf("  D1  = %.4f Н*м,  D2 = %.4f Н*м,  D12 = %.4f Н*м\n",
                p.D1(), p.D2(), p.D12());
    std::printf("  delta1 = %.6f,  delta2 = %.6f\n", p.delta1(), p.delta2());
    std::printf("  delta1/sqrt(delta2) = %.4f   (изотропия дала бы 1)\n",
                p.delta1() / std::sqrt(p.delta2()));
    std::printf("  nu2/nu1 = %.6f  против  E2/E1 = %.6f  -> расхождение %.2f %%\n",
                p.nu2 / p.nu1, p.E2 / p.E1,
                100.0 * std::fabs(p.nu2 / p.nu1 - p.E2 / p.E1) / (p.E2 / p.E1));
    std::printf("  h/min(a,b) = %.4f\n", p.h / std::min(p.a, p.b));
}

// ---------------------------------------------------------------------------
//  Проверка: изотропный предел против табличного решения Тимошенко
//  для квадратной шарнирно опёртой пластины: w_max = 0.00406 * q * a^4 / D
// ---------------------------------------------------------------------------
static void isotropicCheck()
{
    rule("Верификация: изотропный предел (Тимошенко)");
    Plate iso;
    iso.a = 1.0; iso.b = 1.0; iso.h = 0.01;
    iso.E1 = iso.E2 = 2.0e11;
    iso.nu1 = iso.nu2 = 0.3;
    iso.G = iso.E1 / (2.0 * (1.0 + iso.nu1));
    iso.q0 = 1.0e4;

    const double D = iso.D1();
    std::printf("  delta1 = %.10f,  delta2 = %.10f  (должны быть 1)\n",
                iso.delta1(), iso.delta2());

    const auto nav = solveNavier(iso, LoadType::Uniform, 41);
    const double wnum = wCenter(iso, nav);
    const double wref = 0.00406 * iso.q0 * std::pow(iso.a, 4) / D;
    std::printf("  w_max (Навье, K=41) = %.6f мм\n", wnum * 1e3);
    std::printf("  w_max (табличное)   = %.6f мм\n", wref * 1e3);
    std::printf("  расхождение         = %.3f %%  (табличный коэффициент округлён)\n",
                100.0 * std::fabs(wnum - wref) / wref);
}

// ---------------------------------------------------------------------------
static void taskA(const Plate& p)
{
    rule("Задача А: равномерная нагрузка, шарнирное опирание по контуру");

    const int KREF = 41;
    const auto navRef = solveNavier(p, LoadType::Uniform, KREF);
    const double xc = 0.5 * p.a, yc = 0.5 * p.b;

    // изгибающий момент в центре по прогибу
    auto momentX = [&](auto&& s) {
        return -p.D1() * (s.deriv(2, 0, xc, yc) + p.nu2 * s.deriv(0, 2, xc, yc));
    };

    const double wRef  = wCenter(p, navRef);
    const double MxRef = momentX(navRef);
    std::printf("  Эталон: метод Навье, K = %d\n", KREF);
    std::printf("    w_center  = %.6f мм\n", wRef * 1e3);
    std::printf("    Mx_center = %.4f Н*м/м\n\n", MxRef);

    const int KMAX_TRIG = 15, KMAX_POLY = 10;

    std::ofstream csv(OUT + "/convergence_A.csv");
    csv << "K,N,w_trig_mm,err_w_trig_pct,Mx_trig,err_Mx_trig_pct,"
           "w_poly_mm,err_w_poly_pct,Mx_poly,err_Mx_poly_pct,pivot_poly,"
           "offdiag_trig,galerkin_minus_navier_rel\n";
    csv << std::setprecision(10);

    std::cout <<
      "        |        Галёркин, тригонометрический       |        Галёркин, полиномиальный          | Галёркин\n"
      "   K   N|   w, мм    погр.%    Mx, Н*м/м   погр.%   |   w, мм    погр.%    Mx, Н*м/м   погр.%   | vs Навье\n"
      "  ------+------------------------------------------+------------------------------------------+---------\n";

    for (int K = 1; K <= KMAX_TRIG; ++K) {
        const int nq = std::max(60, 5 * K);

        TrigSet tx(K, p.a), ty(K, p.b);
        TensorBasis bt(tx, ty);
        const auto gt = solveGalerkin(p, LoadType::Uniform, bt, nq);
        const double wt = wCenter(p, gt), Mt = momentX(gt);

        // то же усечение ряда Навье -- для прямого сопоставления методов
        const auto nvK = navRef.truncated(K);
        const double wn = wCenter(p, nvK);
        const double dGN = std::fabs(wt - wn) / std::fabs(wRef);

        double wp = 0.0, Mp = 0.0, pivp = 0.0;
        if (K <= KMAX_POLY) {
            PolySet px(K, p.a), py(K, p.b);
            TensorBasis bp(px, py);
            const auto gp = solveGalerkin(p, LoadType::Uniform, bp, nq);
            wp = wCenter(p, gp); Mp = momentX(gp); pivp = gp.pivotRatio;
        }

        const double ewt = 100.0 * std::fabs(wt - wRef)  / std::fabs(wRef);
        const double eMt = 100.0 * std::fabs(Mt - MxRef) / std::fabs(MxRef);
        const double ewp = 100.0 * std::fabs(wp - wRef)  / std::fabs(wRef);
        const double eMp = 100.0 * std::fabs(Mp - MxRef) / std::fabs(MxRef);

        std::printf("  %3d%4d| %9.6f %8.3f %10.4f %8.3f   |", K, K * K,
                    wt * 1e3, ewt, Mt, eMt);
        if (K <= KMAX_POLY)
            std::printf(" %9.6f %8.3f %10.4f %8.3f   |", wp * 1e3, ewp, Mp, eMp);
        else
            std::printf("     --       --        --       --      |");
        std::printf(" %8.1e\n", dGN);

        csv << K << ',' << K * K << ',' << wt * 1e3 << ',' << ewt << ','
            << Mt << ',' << eMt << ',';
        if (K <= KMAX_POLY) csv << wp * 1e3 << ',' << ewp << ',' << Mp << ',' << eMp << ',' << pivp << ',';
        else                csv << ",,,,,";
        csv << gt.offDiagRatio << ',' << dGN << '\n';
    }

    std::cout << "\n  Последний столбец -- относительная разность между методом\n"
                 "  Бубнова-Галёркина и усечённым на том же K рядом Навье.\n";

    // диагональность матрицы Галёркина на тригонометрическом базисе
    {
        TrigSet tx(6, p.a), ty(6, p.b);
        TensorBasis bt(tx, ty);
        const auto g = solveGalerkin(p, LoadType::Uniform, bt, 70);
        std::printf("\n  Тригонометрический базис (K=6): max|K_ij, i!=j| / max|K_ii| = %.3e\n",
                    g.offDiagRatio);
        std::cout << "  -> матрица диагональна с точностью до ошибок округления,\n"
                     "     метод Бубнова-Галёркина вырождается в метод Навье.\n";
    }

    // "наивный" базис
    {
        const double wn = naiveBasisCenterDeflection(p, LoadType::Uniform);
        std::printf("\n  Наивный базис x(a-x)y(b-y): w_center = %.6f мм  (эталон %.6f мм,"
                    " завышение в %.1f раза)\n", wn * 1e3, wRef * 1e3, wn / wRef);
        std::cout << "  -> phi,xxxx = phi,yyyy = 0, поэтому из оператора выпадают ОБЕ\n"
                     "     изгибные жёсткости и остаётся только смешанный член:\n"
                     "         c = q0 / (8*delta1*D1)\n"
                     "     -- прогиб не зависит ни от a, ни от b, ни от delta2.\n"
                     "     Функция не удовлетворяет условию M_n = 0 на контуре, и\n"
                     "     отброшенные при этом граничные члены как раз и несут\n"
                     "     потерянную жёсткость. Отсюда -- балочные полиномы выше.\n";

        // подтверждение независимости от размеров: удвоим сторону a
        Plate p2 = p; p2.a = 2.0 * p.a;
        const double wn2 = naiveBasisCenterDeflection(p2, LoadType::Uniform);
        std::printf("     контроль: при a = %.2f м тот же базис даёт c*phi_max, а сам\n"
                    "     коэффициент c не меняется (w = %.6f мм только из-за роста phi)\n",
                    p2.a, wn2 * 1e3);
    }

    // поля
    const auto F = sampleField(p, navRef, GX, GY);
    writeVTU(F, OUT + "/task_A_uniform.vtu");
    writeSectionCSV(F, OUT + "/task_A_sections.csv");

    double wmax = 0.0, Mxmax = 0.0, Mymax = 0.0;
    for (std::size_t g = 0; g < F.w.size(); ++g) {
        wmax  = std::max(wmax,  std::fabs(F.w[g]));
        Mxmax = std::max(Mxmax, std::fabs(F.Mx[g]));
        Mymax = std::max(Mymax, std::fabs(F.My[g]));
    }
    std::printf("\n  w_max = %.6f мм,  w_max/h = %.4f\n", wmax * 1e3, wmax / p.h);
    std::printf("  |Mx|_max = %.3f Н*м/м,  |My|_max = %.3f Н*м/м\n", Mxmax, Mymax);
}

// ---------------------------------------------------------------------------
static void taskB(const Plate& p)
{
    rule("Задача Б: синусоидальная нагрузка q0*sin(pi x/a)*sin(pi y/b)");

    const double PI = 3.14159265358979323846;
    const double lx = PI / p.a, ly = PI / p.b;
    const double den = std::pow(lx, 4) + 2.0 * p.delta1() * lx * lx * ly * ly
                     + p.delta2() * std::pow(ly, 4);
    const double wExact = p.q0 / (p.D1() * den);
    std::printf("  Точное решение (один член ряда): w_center = %.9f мм\n", wExact * 1e3);

    const auto nav = solveNavier(p, LoadType::SinXY, 8);
    std::printf("  Метод Навье (K=8):               w_center = %.9f мм   [откл. %.2e %%]\n",
                wCenter(p, nav) * 1e3,
                100.0 * std::fabs(wCenter(p, nav) - wExact) / wExact);

    for (int K : {2, 4, 6, 8}) {
        TrigSet tx(K, p.a), ty(K, p.b);
        TensorBasis bt(tx, ty);
        const auto g = solveGalerkin(p, LoadType::SinXY, bt, 70);
        std::printf("  Галёркин, тригон., K=%d:          w_center = %.9f мм   [откл. %.2e %%]\n",
                    K, wCenter(p, g) * 1e3,
                    100.0 * std::fabs(wCenter(p, g) - wExact) / wExact);
    }

    std::cout << "\n  Сходимость полиномиального базиса против ТОЧНОГО решения:\n";
    std::cout << "   K    N   |   w, мм        отн. погрешность   оценка обусл.\n";
    for (int K = 1; K <= 8; ++K) {
        PolySet px(K, p.a), py(K, p.b);
        TensorBasis bp(px, py);
        const auto g = solveGalerkin(p, LoadType::SinXY, bp, 70);
        const double w = wCenter(p, g);
        std::printf("  %3d %4d  |  %10.7f     %10.3e        %9.2e\n",
                    K, K * K, w * 1e3,
                    100.0 * std::fabs(w - wExact) / wExact, g.pivotRatio);
    }

    const auto F = sampleField(p, nav, GX, GY);
    writeVTU(F, OUT + "/task_B_sinxy.vtu");
    writeSectionCSV(F, OUT + "/task_B_sections.csv");
}

// ---------------------------------------------------------------------------
static void taskC(const Plate& p)
{
    rule("Задача В: метод Канторовича-Власова, q = q0*sin(pi y/b)");

    std::printf("  Кромки y=0, y=b: шарнирное опирание (обеспечено множителем sin)\n");
    std::printf("  Кромка  x=0:     %s\n", bcName(EdgeBC::Clamped));
    std::printf("  Кромка  x=a:     %s\n\n", bcName(EdgeBC::Free));

    // --- контрольный прогон: шарнир-шарнир должен совпасть с методом Навье ---
    {
        const auto kv  = solveKantorovich(p, EdgeBC::Simple, EdgeBC::Simple);
        const auto nav = solveNavier(p, LoadType::SinY, 61);
        const double wk = wCenter(p, kv), wn = wCenter(p, nav);
        std::printf("  Контроль (обе кромки шарнирные):\n");
        std::printf("    Канторович-Власов (замкнутая форма) = %.9f мм\n", wk * 1e3);
        std::printf("    Навье (K=61)                        = %.9f мм\n", wn * 1e3);
        std::printf("    расхождение                         = %.3e %%\n\n",
                    100.0 * std::fabs(wk - wn) / std::fabs(wn));
    }

    const auto kv = solveKantorovich(p, EdgeBC::Clamped, EdgeBC::Free);

    std::printf("  Корни характеристического уравнения (1/м):\n");
    for (int k = 0; k < 4; ++k)
        std::printf("    r%d = %+10.5f %+10.5f i\n", k + 1, kv.r[k].real(), kv.r[k].imag());
    std::printf("  Частное решение f_p = %.9f мм\n\n", kv.fp * 1e3);

    // невязка граничных условий
    const double s2 = kv.s * kv.s;
    std::printf("  Проверка граничных условий:\n");
    std::printf("    f(0)                        = %+.3e\n", kv.f(0, 0.0));
    std::printf("    f'(0)                       = %+.3e\n", kv.f(1, 0.0));
    std::printf("    f''(a) - nu2*s^2*f(a)       = %+.3e\n",
                kv.f(2, p.a) - p.nu2 * s2 * kv.f(0, p.a));
    std::printf("    f'''(a) - kappa*s^2*f'(a)   = %+.3e\n\n",
                kv.f(3, p.a) - p.kappa() * s2 * kv.f(1, p.a));

    // невязка самого ОДУ в нескольких точках
    double rmax = 0.0;
    for (int i = 0; i <= 100; ++i) {
        const double x = p.a * i / 100.0;
        const double R = kv.f(4, x) - 2.0 * p.delta1() * s2 * kv.f(2, x)
                       + p.delta2() * s2 * s2 * kv.f(0, x) - p.q0 / p.D1();
        rmax = std::max(rmax, std::fabs(R));
    }
    std::printf("  Максимальная невязка ОДУ по пластине: %.3e (правая часть %.3e)\n\n",
                rmax, p.q0 / p.D1());

    // профиль f(x)
    std::cout << "  Профиль прогиба вдоль x при y = b/2:\n";
    std::cout << "    x/a      w, мм\n";
    for (int i = 0; i <= 10; ++i) {
        const double x = p.a * i / 10.0;
        std::printf("    %.2f   %9.6f\n", i / 10.0, kv.deriv(0, 0, x, 0.5 * p.b) * 1e3);
    }

    const auto F = sampleField(p, kv, GX, GY);
    writeVTU(F, OUT + "/task_C_kantorovich.vtu");
    writeSectionCSV(F, OUT + "/task_C_sections.csv");

    double wmax = 0.0;
    for (double v : F.w) wmax = std::max(wmax, std::fabs(v));
    std::printf("\n  w_max = %.6f мм (на свободной кромке x = a)\n", wmax * 1e3);
}

// ---------------------------------------------------------------------------
int main()
{
    fs::create_directories(OUT);

    Plate p;              // параметры варианта 3 заданы по умолчанию в plate.hpp
    p.check();

    printParams(p);
    isotropicCheck();
    taskA(p);
    taskB(p);
    taskC(p);

    rule("Файлы записаны");
    for (const auto& e : fs::directory_iterator(OUT))
        std::cout << "  " << e.path().string() << '\n';
    return 0;
}
