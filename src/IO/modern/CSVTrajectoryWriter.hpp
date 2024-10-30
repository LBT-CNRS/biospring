#ifndef __CSV_TRAJECTORY_WRITER_HPP__
#define __CSV_TRAJECTORY_WRITER_HPP__

#include "TrajectoryWriterBase.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

class CSVTrajectoryWriter : public TrajectoryWriterBase
{
  public:
    CSVTrajectoryWriter(const std::string & path, const spn::SpringNetwork & topology, size_t write_frequency = 1)
        : TrajectoryWriterBase(path, topology, write_frequency)
    {
        write_header();
    }

    virtual void write_step();

  protected:
    void write_header();
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __CSV_TRAJECTORY_WRITER_HPP__