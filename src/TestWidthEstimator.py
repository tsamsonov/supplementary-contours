import numpy
import arcpy
import sys
import os
# from osgeo import gdal
# from osgeo.gdalconst import *

# print(os.path.dirname(os.path.abspath(__file__)))
# sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/modules")

import WidthEstimator as WidthEstimator


def CalculateWidth(npeuc, cellsize, nodata):
    nrow = npeuc.shape[0]
    ncol = npeuc.shape[1]

    npwid = numpy.full((nrow, ncol), 0)

    # c_npeuc = c_ndarray(npeuc, dtype=numpy.double, ndim=2, shape=(nrow, ncol))
    # c_npwid = c_ndarray(npwid, dtype=numpy.double, ndim=2, shape=(nrow, ncol))

    return WidthEstimator.estimate_width(npeuc, npwid, float(cellsize), float(nodata))

    # return (numpy.ndarray(
    #     shape=(nrow, ncol),
    #     buffer=(ctypes.c_double * nrow * ncol).from_address(ctypes.addressof(c_npwid.data.contents))
    # ))

cell_size = 0.01

test = arcpy.Raster('euc2.tif')
lowerLeft = arcpy.Point(test.extent.XMin, test.extent.YMin)

# npeuc = numpy.array(test.GetRasterBand(1).ReadAsArray())

npeuc = arcpy.RasterToNumPyArray(test)

npwidth = CalculateWidth(npeuc, cell_size, -1)

width = arcpy.NumPyArrayToRaster(npwidth, lowerLeft, cell_size, cell_size, -1)
arcpy.CopyRaster_management(width, 'Y:/GitHub/supplementary-contours/src/width.tif')


#
# driver = test.GetDriver()
# outDs = driver.Create("width.tif", output.shape[1], output.shape[0], 1, GDT_Float32)
# outBand = outDs.GetRasterBand(1)
#
# # write the data
# outBand.WriteArray(output, 0, 0)
#
# # flush data to disk, set the NoData value and calculate stats
# outBand.FlushCache()
# outBand.SetNoDataValue(-99)
#
# # georeference the image and set the projection
# outDs.SetGeoTransform(test.GetGeoTransform())
# outDs.SetProjection(test.GetProjection())
