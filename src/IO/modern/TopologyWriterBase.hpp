#ifndef __TOPOLOGY_WRITER_BASE_HPP__
#define __TOPOLOGY_WRITER_BASE_HPP__

#include "SpringNetwork.h"
#include "WriterBase.hpp"

#include <string>

namespace biospring
{
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