#ifndef __TRAJECTORY_WRITER_BASE_HPP__
#define __TRAJECTORY_WRITER_BASE_HPP__

#include "WriterBase.hpp"
#include <string>

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

class TrajectoryWriterBase : public WriterBase
{
  protected:
    const spn::SpringNetwork & _topology;
    size_t _write_frequency; // Number of steps between frames.
    size_t _current_frame;   // Internal counter of the current frame number.

  public:
    TrajectoryWriterBase(const std::string & path, const spn::SpringNetwork & topology, size_t write_frequency = 1)
        : WriterBase(path), _topology(topology), _write_frequency(write_frequency), _current_frame(0)
    {
        safe_open();
    }

    size_t write_frequency() const { return _write_frequency; }

    virtual ~TrajectoryWriterBase() {}

    virtual void write_step(){};
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __TRAJECTORY_WRITER_BASE_HPP__