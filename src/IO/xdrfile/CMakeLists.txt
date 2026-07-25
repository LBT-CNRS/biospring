add_library(biospring-xdrfile STATIC
    src/xdrfile.c
    src/xdrfile_trr.c
    src/xdrfile_xtc.c
)
add_library(BioSpring::xdrfile ALIAS biospring-xdrfile)

target_include_directories(biospring-xdrfile PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

if(MSVC)
    target_compile_definitions(biospring-xdrfile PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
