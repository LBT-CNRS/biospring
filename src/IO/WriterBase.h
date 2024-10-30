//
// Defines the WriterBase class.
//
// Every Reader class should inherit from those class.
//

#ifndef __WRITER_BASE__
#define __WRITER_BASE__

#include <fstream>
#include <string>

#include "utils.hpp"

namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class WriterBase
{
  public:
    WriterBase() {}
    WriterBase(const std::string & path) : _filename(path) {}
    WriterBase(const char * const path) { setFileName(path); }
    virtual ~WriterBase()
    {
        if (_ostream.is_open())
            _ostream.close();
    }

    void setFileName(const char * const path) { _filename = path; }
    void setFileName(const std::string & path) { _filename = path; }

    std::string & getFileName(void) { return _filename; }
    const std::string & getFileName(void) const { return _filename; }

    virtual void write() = 0;

    // Open file for writing.
    // Dies if something goes wrong.
    void safeOpen()
    {
        if (not _ostream.is_open())
            biospring::utils::file::openwrite(_filename, _ostream);
    }

    bool isOpen() { return _ostream.is_open(); }

  protected:
    std::string _filename;
    std::ofstream _ostream;
};

class TopologyWriterBase : public WriterBase
{
  public:
    TopologyWriterBase(const std::string & path, const biospring::spn::SpringNetwork * const spn)
        : WriterBase(path), _spn(spn)
    {
    }
    TopologyWriterBase(const char * const path, const biospring::spn::SpringNetwork * const spn)
        : WriterBase(path), _spn(spn)
    {
    }

    TopologyWriterBase() : TopologyWriterBase("", nullptr) {}

    // void setSpringNetwork(const biospring::spn::SpringNetwork * const spn) { _spn = spn; }
    // const biospring::spn::SpringNetwork * const getSpringNetwork(void) const { return _spn; }

  protected:
    const biospring::spn::SpringNetwork * _spn;
};

class TrajectoryWriterBase : public TopologyWriterBase
{
  public:
    TrajectoryWriterBase(const std::string & path, const biospring::spn::SpringNetwork * const spn)
        : TopologyWriterBase(path, spn), _timestep(0)
    {
    }
    TrajectoryWriterBase(const char * const path, const biospring::spn::SpringNetwork * const spn)
        : TopologyWriterBase(path, spn), _timestep(0)
    {
    }

    TrajectoryWriterBase() : TrajectoryWriterBase("", nullptr) {}

    size_t getTimeStep() const { return _timestep; }
    size_t getNumberOfStepsBetweenFrames() { return getTimeStep(); }

    void setSpringNetwork(const biospring::spn::SpringNetwork * const spn) { _spn = spn; }

    void setTimeStep(const size_t n) { _timestep = n; }
    void setNumberOfStepsBetweenFrames(const size_t n) { setTimeStep(n); }

    bool isTimeToWrite(size_t currentTime) const { return currentTime % _timestep == 0; }

  protected:
    size_t _timestep; // number of steps between frames
};

#endif