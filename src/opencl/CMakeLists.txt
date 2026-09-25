target_sources(biospring-core PRIVATE SpringNetworkOpenCL.cpp)
target_include_directories(biospring-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(biospring-core PUBLIC OPENCL_SUPPORT)
target_compile_definitions(biospring-core PRIVATE
    CL_TARGET_OPENCL_VERSION=120
    CL_USE_DEPRECATED_OPENCL_1_2_APIS
)
target_link_libraries(biospring-core PUBLIC OpenCL::OpenCL OpenGL::GL)
if(OPENCL_ICD_LOADER)
    # Tells the bundled cl.hpp to take the portable <CL/...> headers rather than
    # Apple's <OpenCL/...>, which is the whole difference on this platform.
    target_compile_definitions(biospring-core PUBLIC BIOSPRING_OPENCL_ICD)
endif()

# The kernel source is compiled into the binary. It used to be copied to the
# top of the build tree and opened at run time through a bare relative path,
# which made the OpenCL backend work only when biospring was launched from
# that one directory. See KernelSource.h.
set(BIOSPRING_CL_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/biospring.cl")
set(BIOSPRING_CL_EMBEDDED "${CMAKE_CURRENT_BINARY_DIR}/BiospringKernelSource.cpp")

# Generated at configure time rather than by a custom command: biospring-core
# is defined in another directory, and a custom command's output is only
# visible to targets declared alongside it. Listing the kernel as a configure
# dependency makes CMake re-run itself when it is edited, so the generated
# file cannot go stale.
#
# The kernel text is the shared force field headers followed by biospring.cl.
# Those headers are compiled a second time here, as OpenCL C, so a formula the
# CPU uses and a formula the GPU uses are the same characters rather than two
# transcriptions that drift. Add a header here when its formulas move into
# src/forcefield/shared/.
set(BIOSPRING_CL_SHARED_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/cellgrid_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/spring_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/electrostatic_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/steric_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/hydrophobic_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/torsion_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/imp_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/hbond_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/random_shared.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/../forcefield/shared/langevin_shared.h"
)

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${BIOSPRING_CL_SOURCE}" ${BIOSPRING_CL_SHARED_HEADERS})

set(BIOSPRING_CL_TEXT "")
foreach(_shared IN LISTS BIOSPRING_CL_SHARED_HEADERS)
    file(READ "${_shared}" _shared_text)
    get_filename_component(_shared_name "${_shared}" NAME)
    string(APPEND BIOSPRING_CL_TEXT
        "// ---- shared with the CPU: ${_shared_name} ----\n${_shared_text}\n")
endforeach()
file(READ "${BIOSPRING_CL_SOURCE}" _kernel_text)
string(APPEND BIOSPRING_CL_TEXT "${_kernel_text}")
configure_file(KernelSource.cpp.in "${BIOSPRING_CL_EMBEDDED}" @ONLY)
target_sources(biospring-core PRIVATE "${BIOSPRING_CL_EMBEDDED}")

if(OPENGL_SUPPORT)
    target_compile_definitions(biospring-core PRIVATE OPENGL_SUPPORT)
    # SpringNetworkOpenCL wraps the viewer's OpenGL VBO in a cl::BufferGL for
    # zero-copy interop when both OpenCL and OpenGL support are enabled; the
    # viewer library must therefore be linked here, not just built.
    target_link_libraries(biospring-core PRIVATE
        OpenGL::GLU GLUT::GLUT GLEW::GLEW BioSpring::viewer)

    if(APPLE)
        # biospring-core now links BioSpring::viewer (above), which itself
        # links BioSpring::core: a genuine mutual dependency between the two
        # static libraries. Apple's linker already resolves this correctly
        # (one copy of the archive is kept), it just warns loudly about it
        # for every executable that links biospring-core. This is Apple
        # ld's own sanctioned flag for exactly this benign case.
        target_link_options(biospring-core PUBLIC "LINKER:-no_warn_duplicate_libraries")
    endif()
endif()
