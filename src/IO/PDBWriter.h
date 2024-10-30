#ifndef __PDBWRITER_H__
#define __PDBWRITER_H__

#include "WriterBase.h"

namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class PDBWriter : public TrajectoryWriterBase
{
  public:
    using TrajectoryWriterBase::TrajectoryWriterBase;

    virtual void write() override;
    virtual void writeModel(size_t modelid);

    void setIsConnect(bool b) { _isconnect = b; }

  protected:
    bool _isconnect;

  private:
};

#endif // __PDBWRITER_H__
