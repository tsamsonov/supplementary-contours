# Supplementary contours

**Supplementary contour lines** are placed between regular contour lines to visualize small but important terrain forms that regular contour lines are unable to show. On topographic maps, typical forms are hillcrests, depressions, saddles, terraces, banks, and levees. 

*Contour types drawn by Eduard Imhof (Cartographic Relief Representation, 1982):*

![Supplementary contours by Eduard Imhof (Cartographic Relief Representation, 1982)](img/contour_types.png)

**supplementary-contours** ArcGIS Python toolbox provides the complete set of tools for generation of supplementary, regular and index contours. You need ArcGIS Pro or ArcGIS 10.3+ to use it.

Processing workflow includes the following steps:

1. Extraction of regular and supplementary contours.
2. Subdivision of the area into the set of regions bordered by regular contours and the boundary of interpolation area.
3. Calculation of region width and centrality rasters.
4. Filtering of supplementary contours' vertices based on width and centrality criteria.
5. Short gap filling, small segment removal and extension of supplementary contour segments.
6. Filtering of closed supplementary contours based on length and average width criteria.
7. Merging the results and flagging of index contours

Six tools are contained in the repository:

1. **Centrality** tool calculates the centrality raster field. The value of each pixel of this raster provides estimation of how close it is located to the central line of the region between contours. Central pixels have *C = 1*,
pixels under the regular contours have *C = 0*.

2. **Region borders** tool subdivides the space into the set of regions formed by regular contours and the boundary of interpolation area. Edges between the resulting regions are returned as a result.

3. **Region width** tool calculates the region width raster field. The value of each pixel of this raster provides estimation of local width of the region. It is modeled using the circle-based approach: the width is equal to the diameter of a largest circle covering the center of the pixel and located completely within the current region.

4. **Supplementary contours** tool takes width, centrality and elevation rasters and produces supplementary contours guided by the set of constraints.

5. **Supplementary contours (full)** tool works similarly to the previous one, but does not require width and centrality rasters as input parameters. It is the main tool of the toolbox which combines all stages of supplementary contours placement.

6. **Width-centrality mask** tool combines width and centrality rasters and produces the raster which masks the suitable areas for supplementary contours generation.

Example digital elevation models are contained in DEM.zip archive available from `data` folder of this repository.

A detailed description of the method can be found in the following paper:

*Samsonov T., Koshel S, Walther D., Jenny B.* Automated placement of supplementary contour lines // **International Journal of Geographical Information Science**. — 2019. — Vol. 33. — DOI: 10.1080/13658816.2019.1610965

(c) Timofey Samsonov & Dmitry Walther, Lomonosov MSU Faculty of Geography, 2017-2019
