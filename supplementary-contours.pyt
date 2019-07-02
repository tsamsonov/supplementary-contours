# -*- coding: utf-8 -*-
# Automated generation of supplementary contours
# 2017-2019 Timofey Samsonov and Dmitry Walther, Lomonosov Moscow State University
#
# Naming convention: function arguments with underscore prefix (_) indicate arcpy Raster objects:
#  width_raster - this is a string path to a raster (e.g. 'in_memory/raster' or 'D:/Folder/raster.tif')
# _width_raster - this is a Raster wrapper around the same object, possibly obtained by arcpy.Raster(width_raster)

import arcpy, numpy, os, sys
import math
from arcpy.sa import *

can_use_cpp = True

if sys.version_info[:2] == (2, 7): # ArcGIS for Desktop 10.3+, Python 2.7 (32 Bit)
    try:
        import WidthEstimator
    except:
        can_use_cpp = False
elif sys.version_info[:2] == (3, 6): # ArcGIS Pro, Python 3.6 (64 Bit)
    try:
        import WidthEstimator3 as WidthEstimator
    except:
        can_use_cpp = False


class Toolbox(object):
    def __init__(self):

        self.label = "Supplementary Contours"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Centrality, RegionWidth, RegionBorders,
                      WidthCentralityMask, SupplementaryContours, SupplementaryContoursFull]

class Centrality(object):
    def __init__(self):

        self.label = "Centrality"
        self.description = "Calculate centrality raster"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_lines = arcpy.Parameter(
            displayName="Input lines",
            name="in_lines",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output centrality raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        snap_raster = arcpy.Parameter(
            displayName="Snap raster",
            name="snap_raster",
            datatype=["GPRasterLayer"],
            parameterType="Optional",
            direction="Input")

        parameters = [in_lines, out_raster, cell_size, snap_raster]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def calculate_centrality(self, in_lines, out_raster, cell_size, _snap_raster):

        def_snap_raster = arcpy.env.snapRaster
        def_extent = arcpy.env.extent

        if (_snap_raster is not None):
            arcpy.env.snapRaster = _snap_raster
            arcpy.env.extent = _snap_raster.extent

        in_dist = EucDistance(in_lines, '', cell_size)
        arcpy.env.snapRaster = in_dist

        # compute geometric parameters
        lleft = arcpy.Point(in_dist.extent.XMin, in_dist.extent.YMin)
        crs = in_dist.spatialReference


        id_field = arcpy.ListFields(in_lines)[0].name
        # TODO: check whether using the first field is reliable

        alloc = EucAllocation(in_lines, cell_size=cell_size, source_field=id_field)

        alloc_regions = "in_memory/alloc_regions"
        arcpy.RasterToPolygon_conversion(alloc, alloc_regions)

        alloc_lines = "in_memory/alloc_lines"
        arcpy.PolygonToLine_management(alloc_regions, alloc_lines)

        alloc_lines_lyr = "alloc_lines_lyr"
        arcpy.MakeFeatureLayer_management(alloc_lines, alloc_lines_lyr, '"LEFT_FID" <> -1')

        in_lines_r =  "in_memory/in_lines_r"

        arcpy.PolylineToRaster_conversion(in_lines,
                                          id_field,
                                          in_lines_r,
                                          cellsize=cell_size)

        in_lines_0 = Con(IsNull(in_lines_r), -1, in_lines_r)

        max = arcpy.GetRasterProperties_management(in_lines_r, "MAXIMUM")
        cost = Reclassify(in_lines_0,
                          reclass_field="Value",
                          remap=RemapRange([[-1, -1, 1],
                                            [0, max, "NODATA"]]))

        centr_dist = CostDistance(alloc_lines_lyr, cost)

        centr_null = in_dist / (in_dist + centr_dist)

        output = CreateConstantRaster(0, "Float", cell_size, centr_null)

        # TODO: put mosaic directly into the folder
        arcpy.Mosaic_management(centr_null, output, "MAXIMUM")

        arcpy.CopyRaster_management(output, out_raster)

        arcpy.env.snapRaster = def_snap_raster
        arcpy.env.extent = def_extent

    def execute(self, parameters, messages):
        in_lines = parameters[0].valueAsText
        out_raster = parameters[1].valueAsText
        cell_size = float(parameters[2].valueAsText.replace(",","."))
        snap_raster = parameters[3].valueAsText

        self.calculate_centrality(in_lines, out_raster, cell_size, arcpy.Raster(snap_raster))

# TODO: overlay width inside closed empty supplementary contours
class RegionWidth(object):
    def __init__(self):

        self.label = "Region width"
        self.description = "Calculate raster of region width"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_features = arcpy.Parameter(
            displayName="Input features",
            name="in_features",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output region width raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        snap_raster = arcpy.Parameter(
            displayName="Snap raster",
            name="snap_raster",
            datatype=["GPRasterLayer"],
            parameterType="Optional",
            direction="Input")

        interest = arcpy.Parameter(
            displayName="Regions of interest",
            name="interest",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        interest.value = 'ALL'
        interest.filter.list = ['ALL', 'INSIDE POLYGONS', 'OUTSIDE POLYGONS']

        mode = arcpy.Parameter(
            displayName="Computation mode",
            name="mode",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        mode.filter.list = ['CPP', 'PYTHON']
        if (can_use_cpp):
            mode.value = 'CPP'
        else:
            mode.value = 'PYTHON'
            mode.enabled = False


        parameters = [in_features, out_raster, cell_size, snap_raster, interest, mode]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    # this function performs at least 10 times faster than the one below
    def calculate_width_circles_cpp(self, npdist, cell_size, nodata = -1):

            nrow = npdist.shape[0]
            ncol = npdist.shape[1]

            npwid = numpy.full((nrow, ncol), 0)

            return WidthEstimator.estimate_width(npdist, npwid, cell_size, nodata)

    def calculate_width_circles(self, npdist, cell_size, nodata = -1):

        N = npdist.shape[0]  # row number
        M = npdist.shape[1]  # column number
        output = numpy.full((N, M), 0)

        for i in range(N):
            for j in range(M):

                radius = npdist[i, j]

                if (radius == nodata):
                    continue

                w = int(math.ceil(radius/cell_size)) # calculate kernel radius (rounded)

                k = range(-w, w+1) # calculate kernel indices

                idx = numpy.meshgrid(k, k) # generate coordinate matrices

                flt = (idx[0]**2 + idx[1]**2 <= w**2) # filter by distance

                x = idx[0][flt] + i
                y = idx[1][flt] + j

                flt_xy = (x >= 0) * (x < N) * (y >= 0) * (y < M) # filter by domain

                x = x[flt_xy]
                y = y[flt_xy]

                flt_distance = output[x, y] < radius*2 # filter by values

                x = x[flt_distance]
                y = y[flt_distance]

                output[x, y] = radius * 2

        return output

    def call(self, in_features, out_raster, cell_size, _snap_raster, interest, mode):

        obstacles = "in_memory/obstacles"

        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(in_features, frame, "ENVELOPE", "ALL")

        buffers = 'in_memory/buffers'

        desc = arcpy.Describe(in_features)

        if desc.shapeType == 'Polygon':
            arcpy.PolygonToLine_management(in_features, obstacles)
            if interest == 'INSIDE POLYGONS':
                arcpy.Buffer_analysis(in_features, buffers, cell_size)
                arcpy.env.mask = buffers
            elif interest == 'OUTSIDE POLYGONS':
                arcpy.Buffer_analysis(in_features, buffers, -cell_size)
                inv_buffers = 'in_memory/inv_buffers'
                arcpy.Erase_analysis(frame, buffers, inv_buffers)
                arcpy.env.mask = inv_buffers
        elif desc.shapeType == 'Polyline':
            arcpy.CopyFeatures_management(in_features, obstacles)
        elif desc.shapeType == 'Point':
            arcpy.Buffer_analysis(in_features, buffers, 0.5 * cell_size)
            arcpy.PolygonToLine_management(buffers, obstacles)
        else:
            return

        if interest != 'INSIDE_POLYGONS':
            frame_lines = "in_memory/lines"
            arcpy.PolygonToLine_management(frame, frame_lines)
            arcpy.Append_management(frame_lines, obstacles, schema_type='NO_TEST')

        # calculate distance raster
        def_snap_raster = arcpy.env.snapRaster
        def_extent = arcpy.env.extent

        if (_snap_raster is not None):
            arcpy.env.snapRaster = _snap_raster
            arcpy.env.extent = _snap_raster.extent

        dist = EucDistance(obstacles, '', cell_size)

        npdist = arcpy.RasterToNumPyArray(dist, nodata_to_value=-1)

        # execute width calculation
        width = self.calculate_width_circles_cpp(npdist, cell_size, -1) if mode == 'CPP' \
            else self.calculate_width_circles(npdist, cell_size, -1)

        # convert to georeferenced raster

        lleft = arcpy.Point(dist.extent.XMin, dist.extent.YMin)
        out = arcpy.NumPyArrayToRaster(width, lleft, cell_size, cell_size, -1)
        arcpy.DefineProjection_management(out, dist.spatialReference)

        if interest == 'INSIDE POLYGONS':
            ExtractByMask(out, in_features).save(out_raster)
        elif interest == 'OUTSIDE POLYGONS':
            inv_features = 'in_memory/inv_features'
            arcpy.Erase_analysis(frame, in_features, inv_features)
            ExtractByMask(out, inv_features).save(out_raster)
        else:
            out.save(out_raster)

        arcpy.env.snapRaster = def_snap_raster
        arcpy.env.extent = def_extent

    def execute(self, parameters, messages):

        in_features = parameters[0].valueAsText
        out_raster = parameters[1].valueAsText
        cell_size = float(parameters[2].valueAsText.replace(",","."))
        snap_raster = parameters[3].valueAsText
        interest = parameters[4].valueAsText
        mode = parameters[5].valueAsText

        self.call(in_features, out_raster, cell_size, snap_raster, interest, mode)

class RegionBorders(object):
    def __init__(self):

        self.label = "Region borders"
        self.description = "Generate region borders"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input elevation raster",
            name="in_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        out_features = arcpy.Parameter(
            displayName="Output region borders feature class",
            name="out_features",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        contour_interval=arcpy.Parameter(
            displayName="Contour interval",
            name="contour_interval",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        base_contour=arcpy.Parameter(
            displayName="Base contour level",
            name="base_contour",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        base_contour.value = 0.0

        parameters = [in_raster, out_features, contour_interval, base_contour]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
        in_raster = parameters[0].valueAsText
        out_features = parameters[1].valueAsText
        contour_interval = float(parameters[2].valueAsText.replace(",","."))
        base_contour = float(parameters[3].valueAsText.replace(",","."))

        main_contours = "in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        # Add fields
        arcpy.AddField_management(main_contours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(main_contours, "Type", "'Regular'", "PYTHON", "")

        addbaselevel = contour_interval / 2
        addcontours = "in_memory/additional_contours"
        Contour(in_raster, addcontours, contour_interval, addbaselevel, 1)

        # Add fields
        arcpy.AddField_management(addcontours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(addcontours, "Type", "'Supplementary'", "PYTHON", "")

        addclosed = "in_memory/addclosed"
        arcpy.FeatureToPolygon_management(addcontours, addclosed, "", "ATTRIBUTES", "")


        addclosedlayer = "addclosedlayer"
        arcpy.MakeFeatureLayer_management(addclosed, addclosedlayer)

        arcpy.SelectLayerByLocation_management(addclosedlayer, 'CONTAINS', main_contours,
                                               "#", "NEW_SELECTION", "INVERT")

        arcpy.SelectLayerByLocation_management(addclosedlayer, 'COMPLETELY_CONTAINS', addcontours,
                                               "#", "SUBSET_SELECTION", "INVERT")

        seladdclosed = "in_memory/selclosed"
        arcpy.CopyFeatures_management(addclosedlayer, seladdclosed)

        addlayer = "addlayer"
        arcpy.MakeFeatureLayer_management(addcontours, addlayer)

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed)

        # generate obstacles
        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(main_contours, frame, "ENVELOPE", "ALL")

        arcpy.PolygonToLine_management(frame, out_features)

        # Add fields
        arcpy.AddField_management(out_features, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(out_features, "Type", "'Boundary'", "PYTHON", "")

        arcpy.Append_management([main_contours, addlayer], out_features, schema_type='NO_TEST')

class WidthCentralityMask(object):
    def __init__(self):

        self.label = "Width-centrality mask"
        self.description = "Generate width-centrality mask"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_width_raster = arcpy.Parameter(
            displayName="Input region width raster",
            name="in_width_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        in_centrality_raster = arcpy.Parameter(
            displayName="Input centrality raster",
            name="in_centrality_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output width-centrality mask raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        width_min = arcpy.Parameter(
            displayName="Region width (minimal)",
            name="width_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_min.value = 0.25

        width = arcpy.Parameter(
            displayName="Region width (optimal)",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width.value = 0.5

        width_max = arcpy.Parameter(
            displayName="Region width (maximal)",
            name="width_max",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_max.value = 0.75

        centrality_min = arcpy.Parameter(
            displayName="Centrality (minimal)",
            name="centrality_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality_min.value = 0.4

        centrality = arcpy.Parameter(
            displayName="Centrality (optimal)",
            name="centrality",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality.value = 0.8

        parameters = [in_width_raster, in_centrality_raster, out_raster,
                      width_min, width, width_max, centrality_min, centrality]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def calculate_mask(self, npwidth, npcentr, wmin, wopt, wmax, cmin, copt):

        N = npwidth.shape[0]  # row number
        M = npwidth.shape[1]  # column number
        output = numpy.zeros((N, M))

        k1 = (copt - cmin) / (wopt - wmin)
        b1 = cmin - wmin * k1

        k2 = (1 - copt) / (wmax - wopt)
        b2 = copt - wopt * k2

        for i in range(N):
            for j in range(M):

                w = npwidth[i, j]
                c = npcentr[i, j]

                flag = 0

                if w < wmin:
                    flag = 0
                elif w < wopt:
                    flag = 1 if c <= w * k1 + b1 else 0
                elif w < wmax:
                    flag = 2 if c <= w * k2 + b2 else 0
                else:
                    flag = 3

                output[i, j] = flag

        return output

    def execute(self, parameters, messages):

        in_width_raster = parameters[0].valueAsText
        in_centrality_raster = parameters[1].valueAsText

        out_raster = parameters[2].valueAsText

        rwmin = float(parameters[3].valueAsText.replace(",", "."))
        rwopt = float(parameters[4].valueAsText.replace(",","."))
        rwmax = float(parameters[5].valueAsText.replace(",","."))

        cmin = float(parameters[6].valueAsText.replace(",","."))
        copt = float(parameters[7].valueAsText.replace(",","."))

        width = arcpy.Raster(in_width_raster)
        centr = arcpy.Raster(in_centrality_raster)

        npwidth = arcpy.RasterToNumPyArray(width)
        npcentr = arcpy.RasterToNumPyArray(centr)

        maxwidth = numpy.amax(npwidth)

        wopt = rwopt * maxwidth
        wmin = rwmin * maxwidth
        wmax = rwmax * maxwidth

        npmask = self.calculate_mask(npwidth, npcentr, wmin, wopt, wmax, cmin, copt)

        lleft = arcpy.Point(width.extent.XMin, width.extent.YMin)
        cell_size = width.meanCellWidth
        arcpy.env.snapRaster = width
        arcpy.env.extent = width.extent

        out = arcpy.NumPyArrayToRaster(npmask, lleft, cell_size)
        arcpy.DefineProjection_management(out, width.spatialReference)

        out.save(out_raster)

class SupplementaryContours(object):
    def __init__(self):

        self.label = "Supplementary contours"
        self.description = "Generate regular, index and supplementary contours from elevation, region width and centrality rasters."
        self.canRunInBackground = True
        self.params = arcpy.GetParameterInfo()

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input elevation raster",
            name="in_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        in_width_raster = arcpy.Parameter(
            displayName="Input region width raster",
            name="in_width_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        in_centrality_raster = arcpy.Parameter(
            displayName="Input centrality raster",
            name="in_centrality_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        out_features = arcpy.Parameter(
            displayName="Output contours feature class",
            name="out_features",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        contour_interval=arcpy.Parameter(
            displayName="Contour interval",
            name="contour_interval",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        base_contour=arcpy.Parameter(
            displayName="Base contour level",
            name="base_contour",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        base_contour.value = 0.0

        index_contour = arcpy.Parameter(
            displayName="Index contour label (each)",
            name="index_contour",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        index_contour.value = 5
        index_contour.filter.type = "Range"
        index_contour.filter.list = [2, 10]

        closed_width_avg = arcpy.Parameter(
            displayName="Closed contour width (average)",
            name="closed_width_avg",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        closed_width_avg.value = 0.125

        width_min = arcpy.Parameter(
            displayName="Region width (minimal)",
            name="width_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_min.value = 0.25

        width=arcpy.Parameter(
            displayName="Region width (optimal)",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width.value = 0.5

        width_max = arcpy.Parameter(
            displayName="Region width (maximal)",
            name="width_max",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_max.value = 0.75

        centrality_min = arcpy.Parameter(
            displayName="Centrality (minimal)",
            name="centrality_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality_min.value = 0.4

        centrality=arcpy.Parameter(
            displayName="Centrality (optimal)",
            name="centrality",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality.value = 0.8

        centrality_ext = arcpy.Parameter(
            displayName="Centrality (extension)",
            name="centrality_ext",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality_ext.value = 0.8

        min_gap=arcpy.Parameter(
            displayName="Gap length (maximal)",
            name="min_gap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_gap.value = 0.5

        min_len=arcpy.Parameter(
            displayName="Segment length (minimal)",
            name="min_len",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_len.value = 1.0

        ext_len = arcpy.Parameter(
            displayName="Extension length",
            name="ext_len",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        ext_len.value = 0.5

        extend = arcpy.Parameter(
            displayName="Extend to defined centrality",
            name="extend",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        extend.value = 'true'

        absolute = arcpy.Parameter(
            displayName="Set parameter values in projection units (absolute values)",
            name="absolute",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        absolute.value = 'false'

        mode = arcpy.Parameter(
            displayName="Region width computation mode",
            name="mode",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        mode.filter.list = ['CPP', 'PYTHON']
        if (can_use_cpp):
            mode.value = 'CPP'
        else:
            mode.value = 'PYTHON'
            mode.enabled = False

        parameters = [in_raster, in_width_raster, in_centrality_raster, out_features,
                      contour_interval, base_contour, index_contour,
                      closed_width_avg, width_min, width, width_max,
                      centrality_min, centrality, centrality_ext,
                      min_gap, min_len, ext_len, extend, absolute, mode]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def calculate_length(self, list1, flaglist):
        list2 = []
        n = 1
        s = 0

        idx = range(len(flaglist) - 1)

        for i in idx:
            if flaglist[i] == flaglist[i + 1]:
                n = n + 1
                s = s + list1[i + 1]
            else:
                s = s + list1[i + 1] * 0.5
                while n > 0:
                    list2.append(s)
                    n = n - 1
                s = list1[i + 1] * 0.5
                n = 1
        while n > 0:
            list2.append(s)
            n = n - 1

        return list2

    def extend_segments(self, list1, flaglist, ext_len, coordinateslist, lowerLeft, cell_size,
                        npcentr, npwidth, width_min, centrality_ext):

        N = len(coordinateslist)
        ni = npwidth.shape[0]

        i1 = 0
        i2 = 0
        while i2 < N - 1:
            if flaglist[i2] == flaglist[i2 + 1]:
                i2 += 1
            elif flaglist[i2] == 1:

                d = 0
                i0 = i1 - 1

                imax = i0
                cmax = 0

                while d <= ext_len:

                    ci = ni - int((coordinateslist[i0][1] - lowerLeft.Y) / cell_size) - 1
                    cj = int((coordinateslist[i0][0] - lowerLeft.X) / cell_size)

                    c = npcentr[ci, cj]
                    w = npwidth[ci, cj]

                    if w <= width_min or i0 == 0:
                        cmax = 1
                        imax = i0
                        break
                    else:
                        # find the highest centrality in length
                        if c >= centrality_ext and c > cmax:
                            imax = i0
                            cmax = c
                        i0 -= 1
                        d += list1[i0]

                if cmax > 0:
                    while imax < i1:
                        flaglist[imax] = 1
                        imax += 1

                d = 0
                i3 = i2 + 1

                imax = i0
                cmax = 0

                while d <= ext_len:
                    ci = ni - int((coordinateslist[i3][1] - lowerLeft.Y) / cell_size) - 1
                    cj = int((coordinateslist[i3][0] - lowerLeft.X) / cell_size)

                    c = npcentr[ci, cj]
                    w = npwidth[ci, cj]

                    if w <= width_min or i3 == N - 1:
                        cmax = 1
                        imax = i3
                        break
                    else:
                        # find the highest centrality in length
                        if c >= centrality_ext and c > cmax:
                            imax = i3
                            cmax = c
                        i3 += 1
                        d += list1[i3]

                if cmax > 0:
                    while imax > i2:
                        flaglist[imax] = 1
                        imax -= 1
                    i2 = i3

                i2 += 1
                i1 = i2

            else:
                i2 += 1
                i1 = i2

        return flaglist

    # TODO: do not filter or extend, if parameters are set to 1 or 0
    # TODO: prohibit filling of the terminal gaps in non-closed contours (?)
    def filter_vertices(self, addlayer, _widthRaster, _centrRaster, width, width_min, width_max,
                        closed_width_avg, min_gap, min_len, ext_len, centrality, centrality_min, centrality_ext, extend):

        cell_size = _widthRaster.meanCellWidth
        lowerLeft = arcpy.Point(_widthRaster.extent.XMin, _widthRaster.extent.YMin)
        crs = _widthRaster.spatialReference  # TODO: remove?

        npwidth = arcpy.RasterToNumPyArray(_widthRaster)
        npcentr = arcpy.RasterToNumPyArray(_centrRaster)

        k1 = (centrality - centrality_min) / (width - width_min)
        b1 = centrality_min - width_min * k1

        k2 = (1 - centrality) / (width_max - width)
        b2 = centrality - width * k2

        cursor = arcpy.da.SearchCursor(addlayer, ['SHAPE@', 'Id', 'Contour'])

        fc_list = []  # Feature coordinates list
        ff_list = []  # Feature flag list
        fi_list = []  # Feature id list
        fh_list = []  # Feature height list

        ni = npwidth.shape[0]

        # filter by width and centrality
        for row in cursor:
            if row[0] != None:
                for part in row[0].getPart():
                    flaglist = []
                    coordinateslist = []
                    for pnt in part:
                        coordinateslist.append([pnt.X, pnt.Y])

                        # Nearest neighbor interpolation
                        i = ni - int((pnt.Y - lowerLeft.Y) / cell_size) - 1
                        j = int((pnt.X - lowerLeft.X) / cell_size)

                        w = npwidth[i, j]
                        c = npcentr[i, j]

                        flag = 0

                        if w < width_min:
                            flag = 0
                        elif w < width:
                            flag = int(c <= w * k1 + b1)
                        elif w < width_max:
                            flag = int(c <= w * k2 + b2)
                        else:
                            flag = 1

                        flaglist.append(flag)

                    fc_list.append(coordinateslist)
                    ff_list.append(flaglist)
                    fi_list.append(row[1])
                    fh_list.append(row[2])

        feature_info = []
        feature_Show = []
        feature_id = []
        feature_height = []

        # filter gaps and segments, extend

        for coordinateslist, flaglist, id, height in zip(fc_list, ff_list, fi_list, fh_list):
            # Fill the cumulative distance
            list1 = [0]
            N = len(coordinateslist)

            idx = range(N - 1)

            idx0 = range(N)

            for i in idx0:
                if i == 0:
                    continue
                dist = (pow(coordinateslist[i][0] - coordinateslist[i - 1][0], 2) +
                        pow(coordinateslist[i][1] - coordinateslist[i - 1][1], 2)) ** 0.5
                list1.append(dist)

            # Calculate segment lengths
            list2 = self.calculate_length(list1, flaglist)

            # Fill gaps
            if numpy.sum(flaglist) > 0.5:  # if this is not a single gap-line
                for i in idx0:
                    if flaglist[i] == 0 and list2[i] <= min_gap:
                        flaglist[i] = 1

            # Recalculate list2 after filling gaps
            list2 = self.calculate_length(list1, flaglist)

            # Remove short segments (switch between 1 and 0)
            for i in idx0:
                if flaglist[i] == 1 and list2[i] <= min_len:
                    flaglist[i] = 0

            # Recalculate list2 after removing short segments
            list2 = self.calculate_length(list1, flaglist)

            # Extend lines
            if extend == 'true' and numpy.sum(flaglist) > 0.5:

                flaglist = self.extend_segments(list1, flaglist, ext_len, coordinateslist,
                                                lowerLeft, cell_size, npcentr, npwidth,
                                                width_min, centrality_ext)

                # Recalculate list2 after extending lines
                list2 = self.calculate_length(list1, flaglist)

                # Fill gaps between extended lines
                for i in idx0:
                    if flaglist[i] == 0 and list2[i] <= min_gap:
                        flaglist[i] = 1

            # Below there is no need to recalculate the length of the segments

            feature = []

            flag = flaglist[0]
            for i in idx0:
                f = flaglist[i]
                p = arcpy.Point(*coordinateslist[i])
                feature.append(p)
                if f != flag or i == N - 1:
                    feature_info.append(arcpy.Polyline(arcpy.Array(feature)))
                    feature_Show.append(flag)
                    feature_id.append(id)
                    feature_height.append(height)
                    feature = [p]
                    flag = f

        return feature_info, feature_Show, feature_id, feature_height, closed_width_avg, width_min

    def process_contours(self, main_contours, addcontours, _width_raster, _centrality_raster, out_features,
                         contour_interval, base_contour, index_contour,
                         closed_width_avg, width_min, width, width_max,
                         centrality_min, centrality, centrality_ext,
                         min_gap, min_len, ext_len, extend, absolute, mode):

        # Add fields
        arcpy.AddField_management(main_contours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(main_contours, "Type", "'Regular'", "PYTHON", "")

        arcpy.AddField_management(main_contours, "Show", "SHORT")
        arcpy.CalculateField_management(main_contours, "Show", 1, "PYTHON", "")

        arcpy.AddField_management(addcontours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(addcontours, "Type", "'Supplementary'", "PYTHON_9.3")

        arcpy.AddField_management(addcontours, "Show", "SHORT")
        arcpy.CalculateField_management(addcontours, "Show", 1, "PYTHON_9.3")

        addlayer = "addlayer"
        arcpy.MakeFeatureLayer_management(addcontours, addlayer)

        addclosed = "in_memory/addclosed"
        arcpy.FeatureToPolygon_management(addcontours, addclosed, "", "ATTRIBUTES", "")

        seladdclosed = "in_memory/selclosed"
        closed_processing = False

        if (int(arcpy.GetCount_management(addclosed).getOutput(0)) > 0):

            arcpy.AddField_management(addclosed, "Show", "SHORT")
            arcpy.CalculateField_management(addclosed, "Show", 1, "PYTHON", "")

            addclosedlayer = "addclosedlayer"
            arcpy.MakeFeatureLayer_management(addclosed, addclosedlayer)

            arcpy.SelectLayerByLocation_management(addclosedlayer, 'CONTAINS', main_contours)
            arcpy.CalculateField_management(addclosedlayer, "Show", 0, "PYTHON", "")

            arcpy.SelectLayerByLocation_management(addclosedlayer, 'COMPLETELY_CONTAINS', addcontours)
            arcpy.CalculateField_management(addclosedlayer, "Show", 0, "PYTHON", "")

            arcpy.SelectLayerByAttribute_management(addclosedlayer, "NEW_SELECTION", '"Show" = 1')
            seladdclosed = "in_memory/selclosed"
            arcpy.CopyFeatures_management(addclosedlayer, seladdclosed)

            if (int(arcpy.GetCount_management(seladdclosed).getOutput(0)) > 0):

                # filter small closed contours by length and average width
                arcpy.AddMessage('PROCESSING CLOSED SUPPLEMENTARY CONTOURS...')
                closed_processing = True
                # cursor = arcpy.da.UpdateCursor(seladdclosed, ['OID@', 'SHAPE@', 'Show'])
                #
                # for row in cursor:
                #     partnum = 0
                #     for part in row[1]:
                #         for pnt in part:
                #             if pnt:
                #                 None
                #             else:
                #                 row[2] = 0
                #                 cursor.updateRow(row)
                #         partnum += 1

                arcpy.AddGeometryAttributes_management(seladdclosed, "PERIMETER_LENGTH")

                arcpy.AddMessage('-- Estimating width...')
                inside_width = 'in_memory/winside'
                widthCalculator = RegionWidth()
                widthCalculator.call(seladdclosed, inside_width, _width_raster.meanCellHeight, _width_raster,
                                     'INSIDE POLYGONS', mode)

                arcpy.AddMessage('-- Filtering lines...')
                zonal_stats = 'in_memory/zonalstats'
                ZonalStatisticsAsTable(seladdclosed, "OBJECTID", inside_width,
                                       zonal_stats, "NODATA", "MEAN")

                arcpy.JoinField_management(seladdclosed, 'OBJECTID', zonal_stats, 'OBJECTID_1', "MEAN")

                seladdclosed_layer = "selected_add_closed_layer"
                arcpy.MakeFeatureLayer_management(seladdclosed, seladdclosed_layer)

                # SMALL
                arcpy.SelectLayerByAttribute_management(seladdclosed_layer, "NEW_SELECTION",
                                                        ' "PERIMETER" <= ' + str(min_len))

                seladdclosed_small = "in_memory/seladdclosed_small"
                arcpy.CopyFeatures_management(seladdclosed_layer, seladdclosed_small)

                arcpy.SelectLayerByLocation_management(addlayer,
                                                       'SHARE_A_LINE_SEGMENT_WITH',
                                                       seladdclosed_small, "",
                                                       "NEW_SELECTION")

                arcpy.CalculateField_management(addlayer, "Show", 0, "PYTHON", "")

                # NARROW
                arcpy.SelectLayerByAttribute_management(seladdclosed_layer, "NEW_SELECTION",
                                                        ' "MEAN" < ' + str(closed_width_avg))

                seladdclosed_narrow = "in_memory/seladdclosed_narrow"
                arcpy.CopyFeatures_management(seladdclosed_layer, seladdclosed_narrow)

                arcpy.SelectLayerByLocation_management(addlayer,
                                                       'SHARE_A_LINE_SEGMENT_WITH',
                                                       seladdclosed_narrow, "",
                                                       "NEW_SELECTION")

                arcpy.CalculateField_management(addlayer, "Show", 0, "PYTHON", "")

                # TODO: apply width criteria to border-attached contours too

                arcpy.SelectLayerByLocation_management(addlayer,
                                                       'SHARE_A_LINE_SEGMENT_WITH',
                                                       seladdclosed, "",
                                                       "NEW_SELECTION",
                                                       "INVERT")

        arcpy.AddMessage('PROCESSING OPEN SUPPLEMENTARY CONTOURS...')
        arcpy.AddMessage('-- Filtering vertices...')

        # TODO: more elegant return with single value
        feature_info, feature_Show, feature_id, feature_height, min_area, width_min = \
            self.filter_vertices(addlayer, _width_raster, _centrality_raster, width,
                                 width_min, width_max, closed_width_avg, min_gap, min_len, ext_len,
                                 centrality, centrality_min, centrality_ext, extend)

        arcpy.AddMessage('-- Reconstructing lines...')

        cursor = arcpy.da.InsertCursor(main_contours, ["SHAPE@", "Type", "Show", "Id", "Contour"])

        for feature, Show, id, height in zip(feature_info, feature_Show, feature_id, feature_height):
            cursor.insertRow([feature, "Supplementary", Show, id, height])


        arcpy.AddMessage('MERGING RESULTS...')

        # Merge all contour types
        if closed_processing:
            arcpy.SelectLayerByLocation_management(addlayer,
                                                   'SHARE_A_LINE_SEGMENT_WITH',
                                                   seladdclosed, "",
                                                   "NEW_SELECTION")
            arcpy.Merge_management([main_contours, addlayer], out_features)
        else:
            arcpy.CopyFeatures_management(main_contours, out_features)

        arcpy.AddField_management(out_features, "Index", "SHORT")

        arcpy.CalculateField_management(out_features, "Index",
                                        "abs(!Contour! - " + str(base_contour) + ") % " +
                                        str(index_contour * contour_interval) + " < " + str(contour_interval / (index_contour + 1)),
                                        "PYTHON", "")

        out_layer = 'out_layer'
        arcpy.MakeFeatureLayer_management(out_features, out_layer)
        arcpy.SelectLayerByAttribute_management(out_layer, 'NEW_SELECTION', '"Index" = 1')
        arcpy.CalculateField_management(out_layer, 'Type', '"Index"')
        arcpy.SelectLayerByAttribute_management(out_layer, 'CLEAR_SELECTION')

        arcpy.DeleteField_management(out_features, 'Index')

        return

    def execute(self, parameters, messages):

        in_raster = parameters[0].valueAsText
        width_raster = parameters[1].valueAsText
        centrality_raster = parameters[2].valueAsText

        out_features = parameters[3].valueAsText

        contour_interval = float(parameters[4].valueAsText.replace(",","."))
        base_contour = float(parameters[5].valueAsText.replace(",","."))
        index_contour = int(parameters[6].valueAsText)

        rclosed_width_avg = float(parameters[7].valueAsText.replace(",", "."))

        rwidth_min = float(parameters[8].valueAsText.replace(",", "."))
        rwidth = float(parameters[9].valueAsText.replace(",","."))
        rwidth_max = float(parameters[10].valueAsText.replace(",","."))

        centrality_min = float(parameters[11].valueAsText.replace(",","."))
        centrality = float(parameters[12].valueAsText.replace(",","."))
        centrality_ext = float(parameters[13].valueAsText.replace(",","."))

        rmin_gap = float(parameters[14].valueAsText.replace(",","."))
        rmin_len = float(parameters[15].valueAsText.replace(",","."))
        rext_len = float(parameters[16].valueAsText.replace(",","."))
        extend = parameters[17].valueAsText
        absolute = parameters[18].valueAsText
        mode = parameters[19].valueAsText

        arcpy.AddMessage('PREPARING CONTOURS...')

        main_contours = "in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        addbaselevel = contour_interval / 2
        addcontours = "in_memory/additional_contours"
        Contour(in_raster, addcontours, contour_interval, addbaselevel, 1)

        # NORMALIZE THRESHOLD VALUES
        wmax = 1
        if absolute == 'false':
            wmax = float(arcpy.GetRasterProperties_management(width_raster, "MAXIMUM").getOutput(0))

        width = rwidth * wmax
        width_min = rwidth_min * wmax
        width_max = rwidth_max * wmax
        closed_width_avg = rclosed_width_avg * wmax

        w = 1
        if absolute == 'false':
            w = width
            arcpy.AddMessage('NORMALIZING THREHOLDS ON\nmax(W) = ' + str(wmax) + ',\nWopt = ' + str(w))

        min_gap = rmin_gap * w
        min_len = rmin_len * w
        ext_len = rext_len * w

        self.process_contours(main_contours, addcontours,
                              arcpy.Raster(width_raster), arcpy.Raster(centrality_raster),
                              out_features,
                              contour_interval, base_contour, index_contour,
                              closed_width_avg, width_min, width, width_max,
                              centrality_min, centrality, centrality_ext,
                              min_gap, min_len, ext_len,
                              extend, absolute, mode)

        return

class SupplementaryContoursFull(object):
    def __init__(self):

        self.label = "Supplementary contours (full)"
        self.description = "Create regular, index and supplementary contours from raster DEM"
        self.canRunInBackground = True
        self.params = arcpy.GetParameterInfo()

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input elevation raster",
            name="in_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        contour_interval=arcpy.Parameter(
            displayName="Contour interval",
            name="contour_interval",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        base_contour=arcpy.Parameter(
            displayName="Base contour level",
            name="base_contour",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        base_contour.value = 0.0

        index_contour = arcpy.Parameter(
            displayName="Index contour label (each)",
            name="index_contour",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        index_contour.value = 5
        index_contour.filter.type = "Range"
        index_contour.filter.list = [2, 10]

        closed_width_avg = arcpy.Parameter(
            displayName="Closed contour width (average)",
            name="closed_width_avg",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        closed_width_avg.value = 0.1

        width_min = arcpy.Parameter(
            displayName="Region width (minimal)",
            name="width_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_min.value = 0.25

        width=arcpy.Parameter(
            displayName="Region width (optimal)",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width.value = 0.5

        width_max = arcpy.Parameter(
            displayName="Region width (maximal)",
            name="width_max",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width_max.value = 0.75

        centrality_min = arcpy.Parameter(
            displayName="Centrality (minimal)",
            name="centrality_min",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality_min.value = 0.4

        centrality=arcpy.Parameter(
            displayName="Centrality (optimal)",
            name="centrality",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality.value = 0.8

        centrality_ext = arcpy.Parameter(
            displayName="Centrality (extension)",
            name="centrality_ext",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality_ext.value = 0.8

        min_gap=arcpy.Parameter(
            displayName="Gap length (maximal)",
            name="min_gap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_gap.value = 0.5

        min_len=arcpy.Parameter(
            displayName="Segment length (minimal)",
            name="min_len",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_len.value = 1.0

        ext_len = arcpy.Parameter(
            displayName="Extension length",
            name="ext_len",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        ext_len.value = 0.5

        out_features = arcpy.Parameter(
            displayName="Output contours feature class",
            name="out_features",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        extend = arcpy.Parameter(
            displayName="Extend to defined centrality",
            name="extend",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        extend.value = 'true'

        absolute = arcpy.Parameter(
            displayName="Set parameter values in projection units (absolute values)",
            name="absolute",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        absolute.value = 'false'

        mode = arcpy.Parameter(
            displayName="Region width computation mode",
            name="mode",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        mode.filter.list = ['CPP', 'PYTHON']
        if (can_use_cpp):
            mode.value = 'CPP'
        else:
            mode.value = 'PYTHON'
            mode.enabled = False

        parameters = [in_raster, out_features, cell_size, contour_interval, base_contour, index_contour, closed_width_avg, width_min, width, width_max,
                      centrality_min, centrality, centrality_ext, min_gap, min_len, ext_len, extend, absolute, mode]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        in_raster = parameters[0].valueAsText
        out_features = parameters[1].valueAsText

        cell_size = float(parameters[2].valueAsText.replace(",","."))
        contour_interval = float(parameters[3].valueAsText.replace(",","."))
        base_contour = float(parameters[4].valueAsText.replace(",","."))
        index_contour = int(parameters[5].valueAsText)

        rclosed_width_avg = float(parameters[6].valueAsText.replace(",", "."))

        rwidth_min = float(parameters[7].valueAsText.replace(",", "."))
        rwidth = float(parameters[8].valueAsText.replace(",","."))
        rwidth_max = float(parameters[9].valueAsText.replace(",","."))

        centrality_min = float(parameters[10].valueAsText.replace(",","."))
        centrality = float(parameters[11].valueAsText.replace(",","."))
        centrality_ext = float(parameters[12].valueAsText.replace(",","."))

        rmin_gap = float(parameters[13].valueAsText.replace(",","."))
        rmin_len = float(parameters[14].valueAsText.replace(",","."))
        rext_len = float(parameters[15].valueAsText.replace(",","."))
        extend = parameters[16].valueAsText
        absolute = parameters[17].valueAsText
        mode = parameters[18].valueAsText

        arcpy.AddMessage('PREPARING CONTOURS...')

        main_contours="in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        addbaselevel = contour_interval / 2
        addcontours = "in_memory/additional_contours"
        Contour(in_raster, addcontours, contour_interval, addbaselevel, 1)

        # generate obstacles
        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(main_contours, frame, "ENVELOPE", "ALL")

        lines = "in_memory/lines"
        arcpy.PolygonToLine_management(frame, lines)
        arcpy.Append_management(main_contours, lines, schema_type='NO_TEST')

        arcpy.AddMessage('ESTIMATING WIDTH...')

        _in_raster = arcpy.Raster(in_raster)
        lowerLeft = arcpy.Point(_in_raster.extent.XMin, _in_raster.extent.YMin)
        crs = _in_raster.spatialReference # TODO: remove?

        arcpy.env.snapRaster = _in_raster
        arcpy.env.extent = _in_raster.extent

        # calculate distance raster
        dist = EucDistance(lines, '', cell_size)
        npdist = arcpy.RasterToNumPyArray(dist)

        # calculate width
        widthCalculator = RegionWidth()

        npwidth = widthCalculator.calculate_width_circles_cpp(npdist, cell_size, -1) if mode == 'CPP' \
            else widthCalculator.calculate_width_circles(npdist, cell_size, -1)

        _width_raster = arcpy.NumPyArrayToRaster(npwidth, lowerLeft, cell_size, cell_size, -1)
        arcpy.DefineProjection_management(_width_raster, dist.spatialReference)

        arcpy.AddMessage('ESTIMATING CENTRALITY...')

        # calculate centrality
        centralityCalculator = Centrality()
        centrality_raster = "in_memory/centr"

        centralityCalculator.calculate_centrality(main_contours, centrality_raster, cell_size, _width_raster)

        # NORMALIZE THRESHOLD VALUES
        wmax = 1
        if absolute == 'false':
            wmax = float(arcpy.GetRasterProperties_management(_width_raster, "MAXIMUM").getOutput(0))

        width = rwidth * wmax
        width_min = rwidth_min * wmax
        width_max = rwidth_max * wmax
        closed_width_avg = rclosed_width_avg * wmax

        w = 1
        if absolute == 'false':
            w = width
            arcpy.AddMessage('NORMALIZING THREHOLDS ON\n-- max(W) = ' + str(wmax) + ',\n-- Wopt = ' + str(w))

        min_gap = rmin_gap * w
        min_len = rmin_len * w
        ext_len = rext_len * w

        mainProcessor = SupplementaryContours()

        mainProcessor.process_contours(main_contours, addcontours,
                                       _width_raster, arcpy.Raster(centrality_raster),
                                       out_features,
                                       contour_interval, base_contour, index_contour,
                                       closed_width_avg, width_min, width, width_max,
                                       centrality_min, centrality, centrality_ext,
                                       min_gap, min_len, ext_len,
                                       extend, absolute, mode)

        return