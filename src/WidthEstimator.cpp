// [[Rcpp::plugins(cpp14)]]

#include<vector>
#include<cmath>
#include<iostream>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>

namespace py = pybind11;
using Shifts = std::vector<std::pair<int, int>>;

// Bresenham's algorithm by X
Shifts get_shifts_byX(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int yi = 1;
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }

    int D = (2 * dy) - dx;
    int y = y0;

    Shifts s;
    for (auto x = x0; x <= x1; x++) {
        s.emplace_back(x, y);
        if (D > 0) {
            y += yi;
            D += 2 * (dy - dx);
        } else {
            D += 2 * dy;
        }
    }

    return s;
}

// Bresenham's algorithm by Y
Shifts get_shifts_byY(int x0, int y0, int x1, int y1) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int xi = 1;
    if (dx < 0) {
        xi = -1;
        dx = -dx;
    }

    int D = (2 * dx) - dy;
    int x = x0;

    Shifts s;
    for (auto y = y0; y <= y1; y++) {
        s.emplace_back(x, y);
        if (D > 0) {
            x += xi;
            D += 2 * (dx - dy);
        } else {
            D += 2 * dx;
        }
    }

    return s;
}

// Bresenham's algorithm in any direction
Shifts get_shifts(int x0, int y0, int x1, int y1) {
    if (abs(y1 - y0) < abs(x1 - x0)) {
        if (x0 > x1)
            return get_shifts_byX(x1, y1, x0, y0);
        else
            return get_shifts_byX(x0, y0, x1, y1);
    } else {
        if (y0 > y1)
            return get_shifts_byY(x1, y1, x0, y0);
        else
            return get_shifts_byY(x0, y0, x1, y1);
    }
}

std::pair<int, int> operator +(const std::pair<int, int>& x, const std::pair<int, int>& y) {
    return std::make_pair(x.first + y.first, x.second + y.second);
}

double dist(int dx, int dy) {
    return sqrt(dx * dx + dy * dy);
}

bool is_within(int i, int j, int imax, int jmax) {
    return (i >= 0) and (j >= 0) and (i < imax) and (j < jmax);
}

// Region length calculation
py::array_t<double> estimate_length(py::array_t<double> cobst, py::array_t<double> clength,
                                    double cellsize, double nodata, int ndirs, double radius) {

    std::vector<Shifts> shifts;
    auto r = int(ceil(radius/cellsize));
    auto angle = 2 * M_PI / ndirs;
    double a = 0;
    int i, j;

    for (auto k = 0; k < ndirs; k++) {
        i = -r * cos(a);
        j =  r * sin(a);
        shifts.push_back(std::move(get_shifts(0, 0, j, i)));
        a += angle;
    }

    auto buf_obst = cobst.request(), buf_length = clength.request();

    if (buf_obst.size != buf_length.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_obst.shape[0];
    auto ncol = buf_obst.shape[1];

    auto *obst = (double *) buf_obst.ptr,
       *length = (double *) buf_length.ptr;

    int ndirs2 = ndirs / 2;

    std::vector<double> d1(ndirs2); // reusable vector;
    std::vector<double> d2(ndirs2); // reusable vector;
    std::vector<double> d12(ndirs2); // reusable vector;

    std::pair<int, int> ij;
    double d;
    int ik, jk;

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            for (auto k = 0; k < ndirs2; k++) {
                d1[k] = r;
                for (auto s: shifts[k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (obst[ik * ncol + jk] == 1) {
                            d1[k] = dist(s.first, s.second);
                            break;
                        }
                    }
                }

                d2[k] = r;
                for (auto s: shifts[ndirs2 + k]) {
                    ik = i + s.first;
                    jk = j + s.second;
                    if (is_within(ik, jk, nrow, ncol)) {
                        if (obst[ik * ncol + jk] == 1) {
                            d2[k] = dist(s.first, s.second);
                            break;
                        }
                    }
                }

                d12[k] = d1[k] + d2[k];
            }

            length[i * ncol + j] = *std::max_element(d12.begin(), d12.end());

        }
    }

    return clength;
}

// Region width calculation
py::array_t<double> estimate_width(py::array_t<double> ceuc, py::array_t<double> cwidth, double cellsize, double nodata) {


    auto buf_euc = ceuc.request(), buf_width = cwidth.request();

    if (buf_euc.size != buf_width.size)
        throw std::runtime_error("Input shapes must match");

    auto nrow = buf_euc.shape[0];
    auto ncol = buf_euc.shape[1];

    auto *euc = (double *) buf_euc.ptr,
         *width = (double *) buf_width.ptr;

    for (int i = 0; i < nrow; i++) {

        for (int j = 0; j < ncol; j++) {

            auto radius = euc[i * ncol + j];

            if (radius == nodata) {
                width[i * ncol + j] = nodata;
                continue;
            }

            auto w = int(ceil(radius / cellsize));

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k * k + l * l > w * w)
                        continue;

                    auto ik = i + k;
                    auto jl = j + l;

                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
                        continue;

                    if (width[ik * ncol + jl] < 2 * radius)
                        width[ik * ncol + jl] = 2 * radius;
                }
            }
        }
    }

    return cwidth;
}

PYBIND11_MODULE(WidthEstimator3, m) {
    m.doc() = R"pbdoc(
        C++ plugin for region width estimation
        -----------------------

        .. currentmodule:: WidthEstimator3

        .. autosummary::
           :toctree: _generate

           estimate_width, estimate_length

    )pbdoc";

    m.def("estimate_width", &estimate_width, R"pbdoc(
        Estimate local region width
    )pbdoc");

    m.def("estimate_length", &estimate_length, R"pbdoc(
        Estimate local region length
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}