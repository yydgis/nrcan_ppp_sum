cmake_minimum_required(VERSION 3.16)
project(nrcan_ppp)

set(CMAKE_CXX_STANDARD 17)

add_executable(nrcan_ppp
nrcan_ppp.cpp
itrf.h
itrf.c
    )

target_link_libraries(nrcan_ppp m)