add_library(biospring-viewer-support STATIC
    SpringNetworkViewer.cpp
    TrackBall.cpp
    Quaternion.cpp
)
add_library(BioSpring::viewer ALIAS biospring-viewer-support)

target_compile_definitions(biospring-viewer-support PUBLIC OPENGL_SUPPORT)
target_include_directories(biospring-viewer-support PUBLIC
    "${PROJECT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(biospring-viewer-support PUBLIC
    BioSpring::core
    OpenGL::GL
    OpenGL::GLU
    GLUT::GLUT
    GLEW::GLEW
)

if(MSVC)
    target_compile_definitions(biospring-viewer-support PRIVATE NOMINMAX _USE_MATH_DEFINES)
endif()
