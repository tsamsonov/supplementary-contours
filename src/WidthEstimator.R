library(Rcpp)
library(stars)

Rcpp::sourceCpp('WidthEstimator.cpp')

euc = read_stars('euc.tif')
plot(euc)

# width = EstimateWidth(euc[[1]], 0.133008, -1)

width = EstimateWidth(euc[[1]], 0.133008, nrow(euc), ncol(euc), -1) %>% 
  matrix(nrow(euc), ncol(euc))

rwidth = st_as_stars(width, dimensions = attr(euc, 'dimensions'))
plot(rwidth)

write_stars(rwidth, 'width.tif')
