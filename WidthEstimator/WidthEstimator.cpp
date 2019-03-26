// [[Rcpp::plugins(cpp14)]]

#include<vector>
#include<cmath>

//#include <Rcpp.h>
//using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> EstimateWidth(std::vector<double> euc, double cellsize, int nrow, int ncol, double nodata) {

    std::vector<double> width(nrow * ncol, 0);

    for (int i = 0; i < nrow; i++){
        for (int j = 0; j < ncol; j++){
            auto radius = euc[i + j*nrow];

            if (radius == nodata){
                width[i + j*nrow] = nodata;
                continue;
            }


            auto w = int(ceil(radius/cellsize));

            for (int k = -w; k <= w; ++k) {
                for (int l = -w; l <= w; ++l) {

                    if (k*k + l*l > w*w)
                        continue;

                    auto ik = i + k;
                    auto jl = j + l;

                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
                        continue;

                    if (width[ik + jl*nrow] < 2*radius)
                        width[ik + jl*nrow] = 2*radius;
                }
            }
        }
    }

    return width;
}

//NumericMatrix EstimateWidth(NumericMatrix euc, double cellsize, double nodata) {
//
//    auto nrow = euc.nrow();
//    auto ncol = euc.ncol();
//
//    NumericMatrix width(nrow, ncol);
//
//    std::fill(width.begin(), width.end(), 0);
//
//    for (int i = 0; i < nrow; i++){
//        for (int j = 0; j < ncol; j++){
//            auto radius = euc(i, j);
//
//            if (radius == nodata){
//                width(i, j) = nodata;
//                continue;
//            }
//
//
//            auto w = int(ceil(radius/cellsize));
//
//            std::vector<int> idx;
//            std::vector<int> jdx;
//
//            for (int k = -w; k <= w; ++k) {
//                for (int l = -w; l <= w; ++l) {
//
//                    if (k*k + l*l > w*w)
//                        continue;
//
//                    auto ik = i + k;
//                    auto jl = j + l;
//
//                    if (ik < 0 || ik >= nrow || jl < 0 || jl >= ncol)
//                        continue;
//
//                    if (width(ik, jl) < 2*radius)
//                        width(ik, jl) = 2*radius;
//                }
//            }
//        }
//    }
//
//    return width;
//}