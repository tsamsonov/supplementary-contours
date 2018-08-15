# -*- coding: utf-8 -*-
import arcpy, numpy, os
import math

from arcpy.sa import *

class Toolbox(object):
    def __init__(self):

        self.label = "Supplementary Contours"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [CalculateCentrality, CalculateWidth, GenerateWidthCentralityMask, SupplContours, SupplContoursFull]

class CalculateCentrality(object):
    def __init__(self):

        self.label = "Local centrality"
        self.description = "Create raster of local centrality"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_lines = arcpy.Parameter(
            displayName="Input lines",
            name="in_lines",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        snap_raster = arcpy.Parameter(
            displayName="Snap raster",
            name="snap_raster",
            datatype=["GPRasterLayer"],
            parameterType="Optional",
            direction="Input")

        parameters = [in_lines, cell_size, out_raster, snap_raster]

        return parameters

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def calculate_centrality(self, in_lines, cell_size, out_raster, snap_raster):

        # calculate distance raster

        if (snap_raster):
            snapRaster = arcpy.Raster(snap_raster)
            arcpy.env.snapRaster = snapRaster
            arcpy.env.extent = snapRaster.extent

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

    def execute(self, parameters, messages):
        in_lines = parameters[0].valueAsText
        cell_size = float(parameters[1].valueAsText.replace(",","."))
        out_raster = parameters[2].valueAsText
        snap_raster = parameters[3].valueAsText

        self.calculate_centrality(in_lines, cell_size, out_raster, snap_raster)
        """The source code of the tool."""

class CalculateWidth(object):
    def __init__(self):

        self.label = "Local width"
        self.description = "Create raster of local region width"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_lines = arcpy.Parameter(
            displayName="Input lines",
            name="in_lines",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="cell_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        out_raster = arcpy.Parameter(
            displayName="Output raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        snap_raster = arcpy.Parameter(
            displayName="Snap raster",
            name="snap_raster",
            datatype=["GPRasterLayer"],
            parameterType="Optional",
            direction="Input")

        parameters = [in_lines, cell_size, out_raster, snap_raster]

        return parameters

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def calculate_width_circles(self, npdist, cell_size):

        N = npdist.shape[0]  # row number
        M = npdist.shape[1]  # column number
        output = numpy.zeros((N, M))

        for i in range(N):
            for j in range(M):
                radius = npdist[i, j]

                w = int(math.ceil(radius/cell_size)) # calculate kernel radius (rounded)

                k = range(-w, w+1) # calculate kernel indices

                idx = numpy.meshgrid(k, k) # generate coordinate matrices

                flt = (idx[0]**2 + idx[1]**2 <= w**2) # filter by distance

                x = idx[0][flt] + i
                y = idx[1][flt] + j

                flt_xy = (x >= 0) * (x < N) * (y >= 0) * (y < M) # filter by domain

                x = x[flt_xy]
                y = y[flt_xy]

                flt_distance = (output[x, y] < radius*2) # filter by values

                x = x[flt_distance]
                y = y[flt_distance]

                output[x, y] = radius * 2

        return output

    def execute(self, parameters, messages):
        in_lines = parameters[0].valueAsText
        cell_size = float(parameters[1].valueAsText.replace(",","."))
        out_raster = parameters[2].valueAsText
        snap_raster = parameters[3].valueAsText

        arcpy.AddMessage('Generating obstacles...')
        # generate obstacles
        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(in_lines, frame, "ENVELOPE", "ALL")

        lines = "in_memory/lines"
        arcpy.PolygonToLine_management(frame, lines)
        arcpy.Append_management(in_lines, lines, schema_type='NO_TEST')

        # calculate distance raster
        if(snap_raster):
            snapRaster = arcpy.Raster(snap_raster)
            arcpy.env.snapRaster = snapRaster
            arcpy.env.extent = snapRaster.extent

        dist = EucDistance(lines, '', cell_size)
        npdist = arcpy.RasterToNumPyArray(dist)

        # execute width calculation
        arcpy.AddMessage('Estimating width...')
        width = self.calculate_width_circles(npdist, cell_size)

        # convert to georeferenced raster
        arcpy.AddMessage('Saving width raster...')
        lleft = arcpy.Point(dist.extent.XMin, dist.extent.YMin)
        out = arcpy.NumPyArrayToRaster(width, lleft, cell_size)
        arcpy.DefineProjection_management(out, dist.spatialReference)

        out.save(out_raster)

class GenerateWidthCentralityMask(object):
    def __init__(self):

        self.label = "Width-centrality mask"
        self.description = "Generate width-centrality mask"
        self.canRunInBackground = True

    def getParameterInfo(self):

        in_width_raster = arcpy.Parameter(
            displayName="Input width raster",
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

        width_min = arcpy.Parameter(
            displayName="Region width (minimum)",
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
            displayName="Region width (maximum)",
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

        out_raster = arcpy.Parameter(
            displayName="Output raster",
            name="out_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Output")

        parameters = [in_width_raster, in_centrality_raster,
                      width_min, width, width_max, centrality_min, centrality,
                      out_raster]

        return parameters

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
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

        rwmin = float(parameters[2].valueAsText.replace(",", "."))
        rwopt = float(parameters[3].valueAsText.replace(",","."))
        rwmax = float(parameters[4].valueAsText.replace(",","."))

        cmin = float(parameters[5].valueAsText.replace(",","."))
        copt = float(parameters[6].valueAsText.replace(",","."))

        out_raster = parameters[7].valueAsText

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

class SupplContours(object):
    def __init__(self):

        self.label = "Supplementary contours"
        self.description = "Create main and additional contours in the most suitable places of map. " \
                           "This version requires width and centrality rasters to be calculated beforehand"
        self.canRunInBackground = True
        self.params = arcpy.GetParameterInfo()

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input raster",
            name="in_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        in_width_raster = arcpy.Parameter(
            displayName="Input width raster",
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

        min_area = arcpy.Parameter(
            displayName="Closed contour area (minimum)",
            name="min_area",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_area.value = 0.1

        width_min = arcpy.Parameter(
            displayName="Region width (minimum)",
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
            displayName="Region width (maximum)",
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
            displayName="Gap length",
            name="min_gap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_gap.value = 0.5

        min_len=arcpy.Parameter(
            displayName="Segment length",
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

        # Выходной параметр — дополнительные горизонтали
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
            parameterType="Optional",
            direction="Input")
        extend.value = 'true'

        parameters = [in_raster, in_width_raster, in_centrality_raster, out_features, contour_interval, base_contour, min_area, width_min, width, width_max,
                      centrality_min, centrality, centrality_ext, min_gap, min_len, ext_len, extend]

        return parameters

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
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
    def filter_vertices(self, addlayer, in_width_raster, in_centrality_raster, rwidth, rwidth_min, rwidth_max,
               rmin_area, rmin_gap, rmin_len, rext_len, centrality, centrality_min, centrality_ext, extend):

        widthRaster = arcpy.Raster(in_width_raster)
        centrRaster = arcpy.Raster(in_centrality_raster)
        cell_size = widthRaster.meanCellWidth
        lowerLeft = arcpy.Point(widthRaster.extent.XMin, widthRaster.extent.YMin)
        crs = widthRaster.spatialReference  # TODO: remove?

        npwidth = arcpy.RasterToNumPyArray(widthRaster)
        npcentr = arcpy.RasterToNumPyArray(centrRaster)

        maxwidth = numpy.amax(npwidth)

        width = rwidth * maxwidth
        width_min = rwidth_min * maxwidth
        width_max = rwidth_max * maxwidth

        min_area = rmin_area * (width ** 2)

        min_gap = rmin_gap * width
        min_len = rmin_len * width
        ext_len = rext_len * width

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
        feature_show = []
        feature_id = []
        feature_height = []

        # filter gaps and segments, extend

        for coordinateslist, flaglist, id, height in zip(fc_list, ff_list, fi_list, fh_list):
            # Заполняем список расстояний от предыдущей точки
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

            # Calculating segment lengths
            list2 = self.calculate_length(list1, flaglist)

            # Filling gaps
            if numpy.sum(flaglist) > 0.5:  # if this is not a single gap-line
                for i in idx0:
                    if flaglist[i] == 0 and list2[i] <= min_gap:
                        flaglist[i] = 1

            # Recalculate list2 after filling gaps
            list2 = self.calculate_length(list1, flaglist)

            # Removing short segments (замена 1 на 0)
            for i in idx0:
                if flaglist[i] == 1 and list2[i] <= min_len:
                    flaglist[i] = 0

            # Recalculate list2 after removing short segments
            list2 = self.calculate_length(list1, flaglist)

            # Extending lines

            if extend == 'true' and numpy.sum(flaglist) > 0.5:

                flaglist = self.extend_segments(list1, flaglist, ext_len, coordinateslist,
                                                lowerLeft, cell_size, npcentr, npwidth,
                                                width_min, centrality_ext)

                # Recalculate list2 after extending lines
                list2 = self.calculate_length(list1, flaglist)

                # Filling gaps between extended lines
                for i in idx0:
                    if flaglist[i] == 0 and list2[i] <= min_gap:
                        flaglist[i] = 1

            # Here there is no need to recalculate the length of the segments

            feature = []

            flag = flaglist[0]
            for i in idx0:
                f = flaglist[i]
                p = arcpy.Point(*coordinateslist[i])
                feature.append(p)
                if f != flag or i == N - 1:
                    feature_info.append(arcpy.Polyline(arcpy.Array(feature)))
                    feature_show.append(flag)
                    feature_id.append(id)
                    feature_height.append(height)
                    feature = [p]
                    flag = f

        return feature_info, feature_show, feature_id, feature_height, min_area

    def process_contours(self, in_raster, in_width_raster, in_centrality_raster, out_features,
                        contour_interval, base_contour, rmin_area, rwidth_min, rwidth, rwidth_max,
                        centrality_min, centrality, centrality_ext, rmin_gap, rmin_len, rext_len, extend):

        arcpy.AddMessage('Preparing contours...')

        main_contours = "in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        addbaselevel = contour_interval / 2
        addcontours = "in_memory/additional_contours"
        Contour(in_raster, addcontours, contour_interval, addbaselevel, 1)

        addclosed = "in_memory/addclosed"
        arcpy.FeatureToPolygon_management(addcontours, addclosed, "", "ATTRIBUTES", "")

        # generate obstacles
        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(main_contours, frame, "ENVELOPE", "ALL")

        lines = "in_memory/lines"
        arcpy.PolygonToLine_management(frame, lines)
        arcpy.Append_management(main_contours, lines, schema_type='NO_TEST')

        # Add fields
        arcpy.AddField_management(main_contours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(main_contours, "Type", "'Regular'", "PYTHON", "")

        arcpy.AddField_management(main_contours, "SHOW", "SHORT")
        arcpy.CalculateField_management(main_contours, "SHOW", 1, "PYTHON", "")

        arcpy.AddField_management(addcontours, "Type", "TEXT", field_length=13)
        arcpy.CalculateField_management(addcontours, "Type", "'Supplementary'", "PYTHON_9.3")

        arcpy.AddField_management(addcontours, "SHOW", "SHORT")
        arcpy.CalculateField_management(addcontours, "SHOW", 1, "PYTHON_9.3")

        arcpy.AddField_management(addclosed, "SHOW", "SHORT")
        arcpy.CalculateField_management(addclosed, "SHOW", 1, "PYTHON", "")

        cursor = arcpy.da.UpdateCursor(addclosed, ['OID@', 'SHAPE@', 'SHOW'])

        for row in cursor:
            partnum = 0
            for part in row[1]:
                for pnt in part:
                    if pnt:
                        None
                    else:
                        row[2] = 0
                        cursor.updateRow(row)
                partnum += 1

        addclosedlayer = "addclosedlayer"
        arcpy.MakeFeatureLayer_management(addclosed, addclosedlayer)

        arcpy.SelectLayerByLocation_management(addclosedlayer, 'CONTAINS', main_contours)
        arcpy.CalculateField_management(addclosedlayer, "SHOW", 0, "PYTHON", "")

        arcpy.SelectLayerByLocation_management(addclosedlayer, 'COMPLETELY_CONTAINS', addcontours)
        arcpy.CalculateField_management(addclosedlayer, "SHOW", 0, "PYTHON", "")

        arcpy.SelectLayerByAttribute_management(addclosedlayer, "NEW_SELECTION", '"SHOW" = 1')
        seladdclosed = "in_memory/selclosed"
        arcpy.CopyFeatures_management(addclosedlayer, seladdclosed)

        addlayer = "addlayer"
        arcpy.MakeFeatureLayer_management(addcontours, addlayer)

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed, "",
                                               "NEW_SELECTION",
                                               "INVERT")

        arcpy.AddMessage('Filtering vertices...')

        # TODO: more elegant return with single value
        feature_info, feature_show, feature_id, feature_height, min_area = \
            self.filter_vertices(addlayer, in_width_raster, in_centrality_raster, rwidth,
                                 rwidth_min, rwidth_max, rmin_area, rmin_gap, rmin_len, rext_len,
                                 centrality, centrality_min, centrality_ext, extend)

        arcpy.AddMessage('Combining results...')

        cursor = arcpy.da.InsertCursor(main_contours, ["SHAPE@", "Type", "SHOW", "Id", "Contour"])

        for feature, show, id, height in zip(feature_info, feature_show, feature_id, feature_height):
            cursor.insertRow([feature, "Supplementary", show, id, height])

        # filter small closed contours by area

        arcpy.AddGeometryAttributes_management(seladdclosed, "AREA")
        seladdclosed_layer = "selected_add_closed_layer"
        arcpy.MakeFeatureLayer_management(seladdclosed, seladdclosed_layer)

        arcpy.SelectLayerByAttribute_management(seladdclosed_layer, "NEW_SELECTION", ' "POLY_AREA" <= ' + str(min_area))

        seladdclosed_small = "in_memory/seladdclosed_small"
        arcpy.CopyFeatures_management(seladdclosed_layer, seladdclosed_small)

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed_small, "",
                                               "NEW_SELECTION")

        arcpy.CalculateField_management(addlayer, "SHOW", 0, "PYTHON", "")

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed, "",
                                               "NEW_SELECTION")

        arcpy.AddMessage('Saving output...')

        # Объединение всех видов горизонталей
        arcpy.Merge_management([main_contours, addlayer], out_features)

        arcpy.AddField_management(out_features, "Index", "SHORT")

        arcpy.CalculateField_management(out_features, "Index",
                                        "abs(!Contour! - " + str(base_contour) + ") % " + str(
                                            5 * contour_interval) + " < " + str(0.25 * contour_interval),
                                        "PYTHON", "")
        return

    def execute(self, parameters, messages):

        in_raster = parameters[0].valueAsText
        in_width_raster = parameters[1].valueAsText
        in_centrality_raster = parameters[2].valueAsText

        out_features = parameters[3].valueAsText

        contour_interval = float(parameters[4].valueAsText.replace(",","."))
        base_contour = float(parameters[5].valueAsText.replace(",","."))

        rmin_area = float(parameters[6].valueAsText.replace(",", "."))

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

        self.process_contours(in_raster, in_width_raster, in_centrality_raster, out_features,
                        contour_interval, base_contour, rmin_area, rwidth_min, rwidth, rwidth_max,
                        centrality_min, centrality, centrality_ext, rmin_gap, rmin_len, rext_len, extend)

        return

class SupplContoursFull(object):
    def __init__(self):

        self.label = "Supplementary contours (full)"
        self.description = "Create main and additional contours in the most suitable places of map"
        self.canRunInBackground = True
        self.params = arcpy.GetParameterInfo()

    def getParameterInfo(self):

        in_raster = arcpy.Parameter(
            displayName="Input raster",
            name="in_raster",
            datatype=["GPRasterLayer"],
            parameterType="Required",
            direction="Input")

        cell_size = arcpy.Parameter(
            displayName="Cell size",
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

        min_area = arcpy.Parameter(
            displayName="Closed contour area (minimum)",
            name="min_area",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_area.value = 0.1

        width_min = arcpy.Parameter(
            displayName="Region width (minimum)",
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
            displayName="Region width (maximum)",
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
            displayName="Gap length",
            name="min_gap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_gap.value = 0.5

        min_len=arcpy.Parameter(
            displayName="Segment length",
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

        # Выходной параметр — дополнительные горизонтали
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
            parameterType="Optional",
            direction="Input")
        extend.value = 'true'

        parameters = [in_raster, out_features, cell_size, contour_interval, base_contour, min_area, width_min, width, width_max,
                      centrality_min, centrality, centrality_ext, min_gap, min_len, ext_len, extend]

        return parameters

    def isLicensed(self):
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
            return False  # tool cannot be executed
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

        rmin_area = float(parameters[5].valueAsText.replace(",", "."))

        rwidth_min = float(parameters[6].valueAsText.replace(",", "."))
        rwidth = float(parameters[7].valueAsText.replace(",","."))
        rwidth_max = float(parameters[8].valueAsText.replace(",","."))

        centrality_min = float(parameters[9].valueAsText.replace(",","."))
        centrality = float(parameters[10].valueAsText.replace(",","."))
        centrality_ext = float(parameters[11].valueAsText.replace(",","."))

        rmin_gap = float(parameters[12].valueAsText.replace(",","."))
        rmin_len = float(parameters[13].valueAsText.replace(",","."))
        rext_len = float(parameters[14].valueAsText.replace(",","."))
        extend = parameters[15].valueAsText

        inRaster = arcpy.Raster(in_raster)
        lowerLeft = arcpy.Point(inRaster.extent.XMin, inRaster.extent.YMin)
        crs = inRaster.spatialReference # TODO: remove?

        arcpy.AddMessage('Preparing contours...')

        main_contours="in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        # generate obstacles
        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(main_contours, frame, "ENVELOPE", "ALL")

        lines = "in_memory/lines"
        arcpy.PolygonToLine_management(frame, lines)
        arcpy.Append_management(main_contours, lines, schema_type='NO_TEST')

        arcpy.AddMessage('Estimating region width...')

        arcpy.env.snapRaster = inRaster
        arcpy.env.extent = inRaster.extent

        # calculate distance raster
        dist = EucDistance(lines, '', cell_size)
        npdist = arcpy.RasterToNumPyArray(dist)


        # calculate width
        widthCalculator = CalculateWidth()
        npwidth = widthCalculator.calculate_width_circles(npdist, cell_size)
        in_width_raster = arcpy.NumPyArrayToRaster(npwidth, lowerLeft, cell_size)

        arcpy.AddMessage('Estimating centrality...')

        # calculate centrality
        centralityCalculator = CalculateCentrality()
        in_centrality_raster = "in_memory/centr"

        centralityCalculator.calculate_centrality(main_contours, cell_size, centr, in_width_raster)

        mainProcessor = SupplContours()

        mainProcessor.execute(in_raster, in_width_raster, in_centrality_raster, out_features,
                        contour_interval, base_contour, rmin_area, rwidth_min, rwidth, rwidth_max,
                        centrality_min, centrality, centrality_ext, rmin_gap, rmin_len, rext_len, extend)

        return