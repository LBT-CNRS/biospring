#ifndef __PDB_TRAJECTORY_WRITER_HPP__
#define __PDB_TRAJECTORY_WRITER_HPP__

#include "TrajectoryWriterBase.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

class PDBTrajectoryWriter : public TrajectoryWriterBase
{
  public:
    using TrajectoryWriterBase::TrajectoryWriterBase;
    virtual void write_step();
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __PDB_TRAJECTORY_WRITER_HPP__