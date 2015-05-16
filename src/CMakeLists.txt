CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(ellipseC)

ADD_EXECUTABLE(ellipse-params centers.h image.h LM.h matrix.h pgm_io.h pgm_io.cpp image.cpp main.cpp matrix.cpp vector.cpp misc.h misc.cpp spline.h spline.cpp)
