CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(ellipseC)

ADD_EXECUTABLE(ellipse-params centers.h image.h LM.h matrix.h pgm_io.h image.cpp main.cpp matrix.cpp vector.cpp misc.h spline.h spline.cpp)
