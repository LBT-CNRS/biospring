#ifndef __XTC_TRAJECTORY_WRITER_HPP__
#define __XTC_TRAJECTORY_WRITER_HPP__

#include "../xdrfile/src/xdrfile.h"
#include "TrajectoryWriterBase.hpp"

#include <string>

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
    //
    // = nullptr MATTERS. This used to be an uninitialised raw pointer, and the
    // constructor below is why nothing ever gave it a value.
    XDRFILE * _xdr = nullptr;

  public:
    // NOT `using TrajectoryWriterBase::TrajectoryWriterBase;`, which is what
    // this class had and what made every XTC run crash.
    //
    // WriterBase::safe_open() is NOT virtual: it opens an std::ofstream, which
    // is exactly right for the PDB and CSV writers and useless here. The base
    // constructor calls it, that call resolves statically to WriterBase's, and
    // this class's own safe_open() -- which opened the xdr handle -- HID the
    // base's rather than overriding it, so nothing on any path ever called it.
    // The file therefore got created and truncated as a text stream nobody
    // wrote to, _xdr kept whatever was on the stack, and write_step() handed
    // that to write_xtc on the first frame.
    //
    // Measured before the fix, on 077.Ribosome (5207 particles) and on 079's
    // 70S (279430), with and without the peptide bond term: SIGSEGV every
    // time, .xtc left at 0 bytes.
    //
    // So the xdr handle is opened HERE, after the base is fully built, and the
    // stream the base opened is closed again: one file, one writer.
    XTCTrajectoryWriter(const std::string & path, const spn::SpringNetwork & topology,
                        size_t write_frequency = 1)
        : TrajectoryWriterBase(path, topology, write_frequency)
    {
        if (_ostream.is_open())
            _ostream.close();
        _open_xdr();
    }

    virtual void write_step();

    ~XTCTrajectoryWriter();

  protected:
    // Named apart from the base's safe_open() on purpose. Sharing that name is
    // what hid this call for as long as it did, and a reader cannot tell a
    // hidden non-virtual from an override by looking at the derived class.
    void _open_xdr();
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __XTC_TRAJECTORY_WRITER_HPP__
