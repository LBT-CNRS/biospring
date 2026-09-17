#ifndef __BIOSPRING_OPENCL_KERNEL_SOURCE_H__
#define __BIOSPRING_OPENCL_KERNEL_SOURCE_H__

namespace biospring
{
namespace opencl
{

// The text of biospring.cl, compiled into the binary by EmbedKernel.cmake.
//
// It used to be read at run time with ifstream("biospring.cl"): a bare
// relative path, resolved against whatever directory biospring was launched
// from. CMake copies the kernel to the top of the build tree, so that worked
// from there and nowhere else -- not from an example directory, not from an
// installed prefix. Carrying the source in the binary removes the lookup, and
// with it the question of where the file has to be.
extern const char * const KERNEL_SOURCE;

} // namespace opencl
} // namespace biospring

#endif // __BIOSPRING_OPENCL_KERNEL_SOURCE_H__
