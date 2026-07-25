target_sources(biospring-core PRIVATE SpringNetworkOpenCL.cpp)
target_include_directories(biospring-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(biospring-core PUBLIC OPENCL_SUPPORT)
target_compile_definitions(biospring-core PRIVATE
    CL_TARGET_OPENCL_VERSION=120
    CL_USE_DEPRECATED_OPENCL_1_2_APIS
)
target_link_libraries(biospring-core PUBLIC OpenCL::OpenCL OpenGL::GL)

configure_file(biospring.cl "${PROJECT_BINARY_DIR}/biospring.cl" COPYONLY)

if(OPENGL_SUPPORT)
    target_compile_definitions(biospring-core PRIVATE OPENGL_SUPPORT)
    target_link_libraries(biospring-core PRIVATE
        OpenGL::GLU GLUT::GLUT GLEW::GLEW)
endif()
