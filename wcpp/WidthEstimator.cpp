// [[Rcpp::plugins(cpp14)]]

#include<vector>
#include<cmath>
#include<iostream>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>

#include "ndarray.h"

namespace py = pybind11;

py::array_t<double> EstimateWidth2(py::array_t<double> ceuc, py::array_t<double> cwidth, double cellsize, double nodata) {


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

void EstimateWidth(numpyArray<double> ceuc, numpyArray<double> cwidth, double cellsize, double nodata) {

    Ndarray<double,2> euc(ceuc);
    Ndarray<double,2> width(cwidth);

    auto nrow = euc.getShape(0);
    auto ncol = euc.getShape(1);

    for (int i = 0; i < nrow; i++) {

        for (int j = 0; j < ncol; j++) {
            auto radius = euc[i][j];

            if (radius == nodata) {
                width[i][j] = nodata;
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

                    if (width[ik][jl] < 2 * radius)
                        width[ik][jl] = 2 * radius;
                }
            }
        }
    }
}

PYBIND11_MODULE(WidthEstimator, m) {
    m.doc() = R"pbdoc(
        C++ plugin for region width estimation
        -----------------------

        .. currentmodule:: WidthEstimator

        .. autosummary::
           :toctree: _generate

           EstimateWidth

    )pbdoc";

    m.def("EstimateWidth2", &EstimateWidth2, R"pbdoc(
        Estimate local region width

        See the following article for details: Samsonov T., Koshel S, Walther D., Jenny B. Automated placement of supplementary contour lines // International Journal of Geographical Information Science. — 2019. — Vol. 33. — DOI: 10.1080/13658816.2019.1610965
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}