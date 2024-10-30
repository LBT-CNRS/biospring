#ifndef __XTC_TRAJECTORY_WRITER_HPP__
#define __XTC_TRAJECTORY_WRITER_HPP__

#include "../xdrfile/src/xdrfile.h"
#include "TrajectoryWriterBase.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

class XTCTrajectoryWriter : public TrajectoryWriterBase
{
  protected:
    // Cannot use modern C++ here, like std::unique_ptr, because the xdrfile library
    // prevents it by design.
    XDRFILE * _xdr;

  public:
    using TrajectoryWriterBase::TrajectoryWriterBase;
    void safe_open();
    virtual void write_step();

    ~XTCTrajectoryWriter();
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __XTC_TRAJECTORY_WRITER_HPP__