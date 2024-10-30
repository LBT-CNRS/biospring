# Enable/disable OpenCL support.
option(OPENCL_SUPPORT "OpenCL support for parallelism on GPUs" OFF)

# Search for OpenCL if required.
if(OPENCL_SUPPORT)
    find_package(OpenCL REQUIRED)
endif(OPENCL_SUPPORT)

if(OpenCL_FOUND)
    target_compile_definitions(biospring-core PUBLIC -DOPENCL_SUPPORT)
    target_link_libraries(biospring-core PUBLIC OpenCL::OpenCL)
endif()


# install OpenCL extension
install(FILES src/biospring.cl DESTINATION ${BIN_INSTALL_DIR})
