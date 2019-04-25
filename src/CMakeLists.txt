cmake_minimum_required(VERSION 3.6)
project(WidthEstimator)

SET (PYBIND11_PYTHON_VERSION 3.6 CACHE STRING "")
SET (PYTHON_EXECUTABLE "C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/python.exe")

add_subdirectory(pybind11)
pybind11_add_module(WidthEstimator WidthEstimator.cpp)