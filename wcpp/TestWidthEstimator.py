import numpy
import ctypes
import sys

sys.path.insert(0, 'C:\Pmodules')

import WidthEstimator

from numpyctypes import c_ndarray

from osgeo import gdal
from osgeo.gdalconst import *


def CalculateWidth(npeuc, cellsize, nodata):
    nrow = npeuc.shape[0]
    ncol = npeuc.shape[1]

    npwid = numpy.full((nrow, ncol), 0)

    c_npeuc = c_ndarray(npeuc, dtype=numpy.double, ndim=2, shape=(nrow, ncol))
    c_npwid = c_ndarray(npwid, dtype=numpy.double, ndim=2, shape=(nrow, ncol))

    return WidthEstimator.EstimateWidth2(npeuc, npwid, float(cellsize), float(nodata))

    # return (numpy.ndarray(
    #     shape=(nrow, ncol),
    #     buffer=(ctypes.c_double * nrow * ncol).from_address(ctypes.addressof(c_npwid.data.contents))
    # ))


test = gdal.Open('euc5.tif')

npeuc = numpy.array(test.GetRasterBand(1).ReadAsArray())

output = CalculateWidth(npeuc, 20, -1)

driver = test.GetDriver()
outDs = driver.Create("width5.tif", output.shape[1], output.shape[0], 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)

# write the data
outBand.WriteArray(output, 0, 0)

# flush data to disk, set the NoData value and calculate stats
outBand.FlushCache()
outBand.SetNoDataValue(-99)

# georeference the image and set the projection
outDs.SetGeoTransform(test.GetGeoTransform())
outDs.SetProjection(test.GetProjection())
