# Supplementary contours

Contours are widely used for graphical representation of continuous fields such as elevation, temperature and gravity. Contour lines are dense in areas of high gradient of the field (where the value grows rapidly) and are sparse in other areas. Therefore, contours depict the areas of high gradient in more detail. To reveal the surface features that are hidden between sparsely placed contours, additional supplementary contours are used. These contours may start and end at
arbitrary points and does not have to fully trace the corresponding isoline. Supplementary contours are drawn only in the areas where their presence gives valuable information to map user and does not produce additional clutter in the image.


This repository contains ArcGIS Python toolbox for automated supplementary contours generation. You need to have ArcGIS 10.3+ to use it, since Python toolboxes were introduced in 10.3 version.

The **processing workflow** is as follows:

1. Regular contours are generated at the selected height interval. These contours subdivide map area into regions.
2. Supplementary contours are generated in each region
3. Supplementray contour vertices are filtered based on width and centrality constraints. Only their sections located in wide regions and that are close to regular contours are kept.
4. Supplementray contour vertices are filtered based on continuity constraints. Line sections separated by short gaps are joined, and remaining short line sections are then removed.
5. Supplementary contours are extended in forward and backward direction until they reach central position in region.

Five tools are contained in the repository:

1. **Local centrality** tool calculates the centrality raster field. The value of each pixel of this raster provides estimation of how close it is located to the central line of the region between contours. Central pixels have *C = 1*,
pixels under the regular contours have *C = 0*.

2. **Local width** tool calculates the region width raster field. The value of each pixel of this raster provides estimation of local width of the region. It is modeled using the circle-based approach: the width is equal to the diameter of a largest circle covering the center of the pixel and located completely within the current region.

3. **Width-centrality mask** tool combines width and centrality rasters and produces the raster which masks the suitable areas for supplementary contours generation.

4. **Supplementary contours** tool takes width, centrality and elevation rasters and produces supplementary contours guided by the set of constraints.

5. **Supplementary contours (full)** tool works similarly to the previous one, but does not require width and centrality rasters as input parameters. Instead it calculates them during processing. But since the calculation of width raster
is very computationally expensive, it is recommended to prepare width and centrality rasters once and then use them to experiment with supplementary contours generation parameters in **Supplementary contours** tool.

Example digital elevation models are contained in DEM.zip archive available from `data` folder of this repository.

(c) Timofey Samsonov & Dmitry Walther, Lomonosov MSU Faculty of Geography, 2017-2019
