#pragma once
// ---------------------------------------------------------------------------
//  Выборка решения на регулярной сетке и вывод в файлы:
//    *.vtu -- XML VTK UnstructuredGrid (ячейки VTK_QUAD) для ParaView
//    *.csv -- таблицы сходимости и сечения
// ---------------------------------------------------------------------------
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "plate.hpp"

struct Field {
    int                 nx = 0, ny = 0;
    std::vector<double> x, y;         // координаты узлов
    std::vector<double> w, Mx, My, Mxy;
};

// Sol должен предоставлять deriv(dx, dy, x, y) для (0,0), (2,0), (0,2), (1,1)
template <class Sol>
Field sampleField(const Plate& p, const Sol& s, int nx, int ny)
{
    Field F;
    F.nx = nx; F.ny = ny;
    F.x.resize(nx); F.y.resize(ny);
    for (int i = 0; i < nx; ++i) F.x[i] = p.a * i / (nx - 1.0);
    for (int j = 0; j < ny; ++j) F.y[j] = p.b * j / (ny - 1.0);

    const std::size_t np = std::size_t(nx) * ny;
    F.w.resize(np); F.Mx.resize(np); F.My.resize(np); F.Mxy.resize(np);

    const double D1 = p.D1(), D2 = p.D2(), D12 = p.D12();
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const std::size_t g = std::size_t(j) * nx + i;
            const double X = F.x[i], Y = F.y[j];
            const double wv  = s.deriv(0, 0, X, Y);
            const double wxx = s.deriv(2, 0, X, Y);
            const double wyy = s.deriv(0, 2, X, Y);
            const double wxy = s.deriv(1, 1, X, Y);
            F.w  [g] = wv;
            F.Mx [g] = -D1 * (wxx + p.nu2 * wyy);
            F.My [g] = -D2 * (wyy + p.nu1 * wxx);
            F.Mxy[g] = -2.0 * D12 * wxy;
        }
    return F;
}

inline void writeVTU(const Field& F, const std::string& path)
{
    std::ofstream f(path);
    f << std::scientific << std::setprecision(10);
    const int nx = F.nx, ny = F.ny;
    const long np = long(nx) * ny;
    const long nc = long(nx - 1) * (ny - 1);

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << nc << "\">\n";

    f << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            f << F.x[i] << ' ' << F.y[j] << " 0\n";
    f << "        </DataArray>\n      </Points>\n";

    f << "      <Cells>\n"
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int j = 0; j < ny - 1; ++j)
        for (int i = 0; i < nx - 1; ++i) {
            const long p0 = long(j) * nx + i;
            f << p0 << ' ' << p0 + 1 << ' ' << p0 + nx + 1 << ' ' << p0 + nx << '\n';
        }
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (long c = 1; c <= nc; ++c) f << 4 * c << (c % 20 ? ' ' : '\n');
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (long c = 0; c < nc; ++c) f << "9" << ((c + 1) % 40 ? ' ' : '\n');
    f << "\n        </DataArray>\n      </Cells>\n";

    auto array = [&](const char* name, const std::vector<double>& v, double scale) {
        f << "        <DataArray type=\"Float64\" Name=\"" << name
          << "\" format=\"ascii\">\n";
        for (long g = 0; g < np; ++g) f << v[g] * scale << '\n';
        f << "        </DataArray>\n";
    };

    f << "      <PointData Scalars=\"w_mm\">\n";
    array("w_mm",   F.w,   1000.0);   // прогиб в мм
    array("w_m",    F.w,   1.0);      // прогиб в м (для Warp By Scalar)
    array("Mx",     F.Mx,  1.0);      // Н*м/м
    array("My",     F.My,  1.0);
    array("Mxy",    F.Mxy, 1.0);
    f << "      </PointData>\n";

    f << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
}

// сечение вдоль x при y = yFix (и наоборот) -> CSV
inline void writeSectionCSV(const Field& F, const std::string& path)
{
    std::ofstream f(path);
    f << std::setprecision(10);
    f << "x,w_mm_at_ymid,y,w_mm_at_xmid\n";
    const int jm = F.ny / 2, im = F.nx / 2;
    const int n  = std::max(F.nx, F.ny);
    for (int k = 0; k < n; ++k) {
        if (k < F.nx) f << F.x[k] << ',' << F.w[std::size_t(jm) * F.nx + k] * 1000.0;
        else          f << ",";
        if (k < F.ny) f << ',' << F.y[k] << ',' << F.w[std::size_t(k) * F.nx + im] * 1000.0;
        else          f << ",,";
        f << '\n';
    }
}
