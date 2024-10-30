# CMakefile for xdrfile-1.1.1
# author: Benoist LAURENT
# version: 1.0
# date: 08/02/2011


set(XDRFILE_SRCS
    src/xdrfile.c
    src/xdrfile_trr.c
    src/xdrfile_xtc.c
    )

add_library(xdrfile STATIC ${XDRFILE_SRCS})

target_include_directories(xdrfile INTERFACE src/)
