#ifndef __TESTS_TEMPPATH_H__
#define __TESTS_TEMPPATH_H__

#include <atomic>
#include <filesystem>
#include <random>
#include <string>

namespace biospring
{
namespace tests
{

// A unique path in the system temporary directory, ending in `suffix`.
//
// WHY NOT std::tmpnam. It is deprecated on every current toolchain, and the
// reason is real: the name it returns is not reserved, so anything may claim
// the path between the call and the open. mkstemp fixes that by creating the
// file, but it cannot place anything AFTER its random part, and these tables
// are chosen by their extension -- a reader keys on .hbond, .rbody or .psite --
// so the suffix is not decoration.
//
// Uniqueness here comes from a random 64-bit draw plus a counter, so repeated
// calls within one test cannot collide with each other the way some tmpnam
// implementations do when called more than once.
inline std::string temp_path(const std::string & suffix)
{
    static std::atomic<unsigned long> counter{0};
    static std::random_device device;
    static const unsigned long salt = (static_cast<unsigned long>(device()) << 32) ^ device();
    return (std::filesystem::temp_directory_path() /
            ("biospring-test-" + std::to_string(salt) + "-" + std::to_string(counter++)))
               .string() +
           suffix;
}

} // namespace tests
} // namespace biospring

#endif // __TESTS_TEMPPATH_H__
