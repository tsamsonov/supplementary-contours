cmake_minimum_required(VERSION 2.8.12)
project(WidthEstimator3)

#SET (PYBIND11_PYTHON_VERSION 3.7 CACHE STRING "")
#SET (PYTHON_EXECUTABLE "C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/python.exe")
#add_subdirectory(pybind11)
#pybind11_add_module(WidthEstimator3 WidthEstimator.cpp)

# macOS
SET (PYBIND11_PYTHON_VERSION 3.7 CACHE STRING "")
SET (PYTHON_EXECUTABLE "/usr/local/Cellar/python/3.7.2_1/bin/python3.7")
add_subdirectory(pybind11)
pybind11_add_module(WidthEstimator3 WidthEstimator.cpp)


#SET (PYBIND11_PYTHON_VERSION 2.7 CACHE STRING "")
#SET (PYTHON_EXECUTABLE "C:/Python27/ArcGIS10.6/python.exe")
#add_subdirectory(pybind11)
#pybind11_add_module(WidthEstimator WidthEstimator.cpp)