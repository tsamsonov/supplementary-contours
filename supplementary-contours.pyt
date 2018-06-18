# -*- coding: utf-8 -*-
import arcpy, numpy, os
import math

from arcpy.sa import *

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Supplementary Contours"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [CalculateCentrality, CalculateWidth, SupplContours]

class CalculateCentrality(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Local centrality"
        self.description = "Create raster of local centrality"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

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

        parameters = [in_lines, cell_size, out_raster]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def calculate_centrality(self, in_lines, cell_size, out_raster):

        # calculate distance raster
        in_dist = EucDistance(in_lines, '', cell_size)
        arcpy.env.snapRaster = in_dist

        # compute geometric parameters
        lleft = arcpy.Point(in_dist.extent.XMin, in_dist.extent.YMin)
        crs = in_dist.spatialReference

        alloc = EucAllocation(in_lines, cell_size=cell_size, source_field="OBJECTID")

        alloc_regions = "in_memory/alloc_regions"
        arcpy.RasterToPolygon_conversion(alloc, alloc_regions)

        alloc_lines = "in_memory/alloc_lines"
        arcpy.PolygonToLine_management(alloc_regions, alloc_lines)

        alloc_lines_lyr = "alloc_lines_lyr"
        arcpy.MakeFeatureLayer_management(alloc_lines, alloc_lines_lyr, '"LEFT_FID" <> -1')

        in_lines_r =  "in_memory/in_lines_r"
        # arcpy.AddMessage([field.name for field in arcpy.ListFields(in_lines)])

        # TODO: replace FID with calculated OID value
        arcpy.PolylineToRaster_conversion(in_lines, "FID", in_lines_r)

        in_lines_0 = Con(IsNull(in_lines_r), -1, in_lines_r)

        max = arcpy.GetRasterProperties_management(in_lines_r, "MAXIMUM")
        cost = Reclassify(in_lines_0,
                          reclass_field="Value",
                          remap=RemapRange([[-1, -1, 1],
                                            [0, max, "NODATA"]]))

        centr_dist = CostDistance(alloc_lines_lyr, cost)

        centr_null = in_dist / (in_dist + centr_dist)

        output = CreateConstantRaster(0, "Float", cell_size)
        arcpy.Mosaic_management(centr_null, output, "MAXIMUM")

        output.save(out_raster)

    def execute(self, parameters, messages):
        in_lines = parameters[0].valueAsText
        cell_size = float(parameters[1].valueAsText.replace(",","."))
        out_raster = parameters[2].valueAsText

        self.calculate_centrality(in_lines, cell_size, out_raster)
        """The source code of the tool."""

class CalculateWidth(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Local width"
        self.description = "Create raster of local region width"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

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

        parameters = [in_lines, cell_size, out_raster]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def calculate_width_circles(self, in_lines, cell_size, out_raster):

        frame = "in_memory/frame"
        arcpy.MinimumBoundingGeometry_management(in_lines, frame, "ENVELOPE", "ALL")

        lines = "in_memory/lines"
        arcpy.PolygonToLine_management(frame, lines)
        arcpy.Append_management(in_lines, lines, schema_type='NO_TEST')

        # calculate distance raster
        dist = EucDistance(lines, '', cell_size)

        # compute geometric parameters
        lleft = arcpy.Point(dist.extent.XMin, dist.extent.YMin)
        crs = dist.spatialReference

        npdist = arcpy.RasterToNumPyArray(dist)

        N = npdist.shape[0]  # число строк
        M = npdist.shape[1]  # число столбцов
        output = numpy.zeros((N, M))

        for i in range(N):
            for j in range(M):
                radius = npdist[i, j]

                w = int(math.ceil(radius/cell_size)) # calculate kernel radius (rounded)

                arcpy.AddMessage(w)
                arcpy.AddMessage(cell_size)
                break

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

        out = arcpy.NumPyArrayToRaster(output, lleft, cell_size)

        arcpy.DefineProjection_management(out, crs)

        out.save(out_raster)

    def execute(self, parameters, messages):
        in_lines = parameters[0].valueAsText
        cell_size = float(parameters[1].valueAsText.replace(",","."))
        out_raster = parameters[2].valueAsText

        self.calculate_width_circles(in_lines, cell_size, out_raster)
        """The source code of the tool."""

class SupplContours(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Supplementary contours"
        self.description = "Create main and additional contours in the most suitable places of map"
        self.canRunInBackground = True
        self.params = arcpy.GetParameterInfo()

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Входной парметр — ЦМР
        in_raster = arcpy.Parameter(
            displayName="Input raster",
            name="in_raster",
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

        width=arcpy.Parameter(
            displayName="Region width (minimum)",
            name="width",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        width.value = 0.0

        centrality=arcpy.Parameter(
            displayName="Centrality (maximum)",
            name="centrality",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        centrality.value = 1.0

        min_gap=arcpy.Parameter(
            displayName="Gap length (minimum)",
            name="min_gap",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_gap.value = 0.0

        min_len=arcpy.Parameter(
            displayName="Segment length (minimum)",
            name="min_len",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_len.value = 0.0

        # Входной параметр — минимальная длина дополнительных горизонталей
        min_area = arcpy.Parameter(
            displayName="Closed contour area (minimum)",
            name="min_area",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        min_area.value = 0.0

        # Выходной параметр — дополнительные горизонтали
        out_features = arcpy.Parameter(
            displayName="Output contours feature class",
            name="out_features",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        parameters = [in_raster, contour_interval, base_contour, width, centrality, min_gap, min_len, min_area, out_features]

        return parameters

    def isLicensed(self):
        # try:
        #     if arcpy.CheckExtension("Spatial") != "Available":
        #         raise Exception
        # except Exception:
        #     return False  # tool cannot be executed
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        in_raster = parameters[0].valueAsText
        contour_interval = float(parameters[1].valueAsText.replace(",","."))
        base_contour = float(parameters[2].valueAsText.replace(",","."))
        width = float(parameters[3].valueAsText.replace(",","."))
        centrality = float(parameters[4].valueAsText.replace(",","."))
        min_gap = float(parameters[5].valueAsText.replace(",","."))
        min_len = float(parameters[6].valueAsText.replace(",","."))
        min_area = float(parameters[7].valueAsText.replace(",","."))
        out_features = parameters[8].valueAsText

        inRaster = arcpy.Raster(in_raster)
        lowerLeft = arcpy.Point(inRaster.extent.XMin, inRaster.extent.YMin)
        cellSize = inRaster.meanCellWidth
        crs = inRaster.spatialReference

        # Построение основных горизонталей
        main_contours="in_memory/main_contours"
        Contour(in_raster, main_contours, contour_interval, base_contour, 1)

        # Добавление поля Type к основным горизонталям
        arcpy.AddField_management(main_contours, "Type", "TEXT")
        arcpy.CalculateField_management(main_contours, "Type", "'Main'", "PYTHON", "")

        arcpy.AddField_management(main_contours, "SHOW", "SHORT")
        arcpy.CalculateField_management(main_contours, "SHOW", 1, "PYTHON", "")

        width_raster = "in_memory/width_raster"
        widthCalculator = CalculateWidth()
        widthCalculator.calculate_width_circles(main_contours, cellSize, width_raster)

        # Построение дополнительных горизонталей
        addbaselevel = contour_interval/2
        addcontours="in_memory/aditional_contours"
        Contour(in_raster, addcontours, contour_interval,  addbaselevel, 1)

        # Добавление поля Type к дополнительным горизонталям
        arcpy.AddField_management(addcontours, "Type", "TEXT")
        arcpy.CalculateField_management(addcontours, "Type", "'Add'", "PYTHON_9.3")

        arcpy.AddField_management(addcontours, "SHOW", "SHORT")
        arcpy.CalculateField_management(addcontours, "SHOW", 1, "PYTHON_9.3")

        addclosed = "in_memory/addclosed"
        arcpy.FeatureToPolygon_management(addcontours, addclosed, "", "ATTRIBUTES", "")

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

        arcpy.SelectLayerByAttribute_management(addclosedlayer, "NEW_SELECTION", '"SHOW" = 1')
        seladdclosed = "in_memory/sellines"
        arcpy.CopyFeatures_management(addclosedlayer, seladdclosed)

        seladdclosed_layer = "selected_add_closed_layer"
        arcpy.MakeFeatureLayer_management(seladdclosed, seladdclosed_layer)

        addlayer = "addlayer"
        arcpy.MakeFeatureLayer_management(addcontours, addlayer)

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed_layer, "",
                                               "NEW_SELECTION",
                                               "INVERT")

        addcontours_interp = "in_memory/addcontours_interp"

        arcpy.InterpolateShape_3d(width_raster, addlayer, addcontours_interp)

        cursor = arcpy.da.SearchCursor(addcontours_interp, ['SHAPE@'])

        fc_list = [] # Feature coordinates list
        ff_list = [] # Festure flag list

        # заполняем  список флагов (0,1) и координат
        for row in cursor:
            for part in row[0].getPart():
                flaglist = []
                coordinateslist = []
                for pnt in part:
                    coordinateslist.append([pnt.X, pnt.Y])
                    if pnt.Z >= width:
                        flaglist.append(1)
                    else:
                        flaglist.append(0)
                fc_list.append(coordinateslist)
                ff_list.append(flaglist)

        feature_info = []
        feature_show = []
        for coordinateslist, flaglist in zip(fc_list, ff_list):
            # Заполняем список расстояний от предыдущей точки
            list1 = [0]
            N = len(coordinateslist)

            # arcpy.AddMessage(N)

            idx = range(N-1)

            idx0 = range(N)

            for i in idx0:
                if i == 0:
                    continue
                dist = (pow(coordinateslist[i][0] - coordinateslist[i-1][0], 2) +
                        pow(coordinateslist[i][1] - coordinateslist[i-1][1], 2)) ** 0.5
                list1.append(dist)

            # Заполняем список длин однородных участков
            list2 = []
            n = 1
            s = 0

            for i in idx:

                if flaglist[i] == flaglist[i+1]:
                    n = n + 1
                    s = s + list1[i + 1]

                else:
                    s = s + list1[i + 1] * 0.5
                    while n > 0:
                        list2.append(s)
                        n = n-1
                    s = list1[i + 1] * 0.5
                    n = 1

            # Обработка последней точки, если она отличается от предпоследней
            while n > 0:
                list2.append(s)
                n = n - 1

            # Удаление дыр (замена 0 на 1)
            for i in idx:
                if flaglist[i] == 0 and list2[i] <= min_gap:
                    flaglist[i] = 1

            # Пересчёт list2 после удалений
            list2 = []
            n = 1
            s = 0

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

            # Удаление коротких отрезков (замена 1 на 0)
            for i in idx:
                if flaglist[i] == 1 and list2[i] <= min_len:
                    flaglist[i] = 0

            # Пересчёт list2 после удалений
            list2 = []
            n = 1
            s = 0

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

            feature = []

            flag = flaglist[0]
            for i in idx0:
                f = flaglist[i]
                p = arcpy.Point(*coordinateslist[i])
                feature.append(p)
                if f != flag or i == N-1:
                    feature_info.append(arcpy.Polyline(arcpy.Array(feature)))
                    feature_show.append(flag)
                    feature = [p]
                    flag = f

        cursor = arcpy.da.InsertCursor(main_contours, ["SHAPE@", "Type", "SHOW"])

        for feature, show in zip(feature_info, feature_show):
            cursor.insertRow([feature, "Add", show])

        # Выделение всех замкнутых дополнительных горизонталей

        arcpy.AddGeometryAttributes_management(seladdclosed_layer, "AREA")

        # arcpy.AddMessage([f.name for f in arcpy.ListFields(seladdclosed_layer)])

        arcpy.SelectLayerByAttribute_management(seladdclosed_layer, "NEW_SELECTION", ' "POLY_AREA" >= ' + str(min_area))

        seladdclosed_layer_copy = "in_memory/seladdclosed_layer_copy"
        arcpy.CopyFeatures_management(seladdclosed_layer, seladdclosed_layer_copy)

        arcpy.SelectLayerByLocation_management(addlayer,
                                               'SHARE_A_LINE_SEGMENT_WITH',
                                               seladdclosed_layer_copy, "",
                                               "NEW_SELECTION")

        # Объединение всех видов горизонталей
        arcpy.Merge_management([main_contours, addlayer], out_features)

        return