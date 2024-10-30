# FindMDDriver.cmake v2.0 07/2023 Benoist LAURENT & Hubert Santuz
#
# This cmake module aims to find the MDDriver depencies in the context
# of the BioSpring project.
#
# Sets:
#     - MDDRIVER_LIBRARY
#     - MDDRIVER_INCLUDE_DIR
# Exports also the target:
#     - MDDriver::MDDriver

include(FindPackageHandleStandardArgs)

find_package(MDDriver NO_MODULE)

find_package_handle_standard_args(MDDriver
    REQUIRED_VARS MDDRIVER_LIBRARY MDDRIVER_INCLUDE_DIR
    FAIL_MESSAGE "MDDriver not found. Set the MDDriver_DIR cmake cache entry to the directory containing MDDriverConfig.cmake."
)

if (MDDriver_FOUND)
    message(STATUS "MDDriver found")
    message(STATUS "MDDriver library: ${MDDRIVER_LIBRARY}")
    message(STATUS "MDDriver include directory: ${MDDRIVER_INCLUDE_DIR}")
endif()

if (MDDriver_FOUND AND NOT TARGET MDDriver::MDDriver)
  add_library(MDDriver::MDDriver STATIC IMPORTED)
  set_property(TARGET MDDriver::MDDriver PROPERTY IMPORTED_LOCATION ${MDDRIVER_LIBRARY})
  target_include_directories(MDDriver::MDDriver INTERFACE ${MDDRIVER_INCLUDE_DIR})
endif()
