#ifndef __TOPOLOGY_WRITER_BASE_HPP__
#define __TOPOLOGY_WRITER_BASE_HPP__

#include "WriterBase.hpp"

#include <string>

// Forward-declared, not included: SpringNetwork.h includes IO/modern.hpp,
// which lands back here, and this header only needs to hold a reference.
// The cycle used to be broken by luck -- reduce/Reduce.h happened to pull
// SpringNetwork.h in first -- and removing that legacy header exposed it.
namespace biospring
{
namespace spn
{
class SpringNetwork;
}
namespace io
{
namespace modern
{

class TopologyWriterBase : public WriterBase
{
  protected:
    const spn::SpringNetwork & _topology;

  public:
    TopologyWriterBase(const std::string & path, const spn::SpringNetwork & topology)
        : WriterBase(path), _topology(topology)
    {
    }

    virtual void write() = 0;
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __TOPOLOGY_WRITER_BASE_HPP__