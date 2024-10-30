

#ifndef __XTCTRAJWRITER_H__
#define __XTCTRAJWRITER_H__

#include "WriterBase.h"
#include "xdrfile.h"
#include <string>

namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class XTCTrajWriter : public TrajectoryWriterBase
{
  public:
    XTCTrajWriter(const std::string & path, const biospring::spn::SpringNetwork * const spn)
        : TrajectoryWriterBase(path, spn), _xdr(nullptr), _step(0)
    {
        _timestep = 0;
    }
    XTCTrajWriter(const char * const path, const biospring::spn::SpringNetwork * const spn)
        : TrajectoryWriterBase(path, spn), _xdr(nullptr), _step(0)
    {
        _timestep = 0;
    }
    XTCTrajWriter() : XTCTrajWriter("", nullptr) {}

    ~XTCTrajWriter()
    {
        if (_xdr)
            xdrfile_close(_xdr);
    }

    void write() override { writeNextStep(); }

    void update();
    void close();
    void writeNextStep();
    void safeOpen();

  protected:
    XDRFILE * _xdr;
    size_t _step;
};

#endif // __XTCTRAJWRITER_H__
