#ifndef __PDBTRAJECTORYWRITER_H__
#define __PDBTRAJECTORYWRITER_H__

#include "PDBWriter.h"

namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class PDBTrajectoryWriter : public PDBWriter
{
  public:
    PDBTrajectoryWriter(const std::string & path, const biospring::spn::SpringNetwork * const spn)
        : PDBWriter(path, spn), _frameid(1)
    {
        _timestep = 100;
    }
    PDBTrajectoryWriter(const char * const path, const biospring::spn::SpringNetwork * const spn)
        : PDBWriter(path, spn), _frameid(1)
    {
        _timestep = 100;
    }
    PDBTrajectoryWriter() : PDBTrajectoryWriter("", nullptr) {}

    void write() override;

  protected:
    size_t _frameid;
};

#endif // __PDBTRAJECTORYWRITER_H__
