//
// Defines the ReaderBase and TopologyReaderBase classes.
//
// Every Reader class should inherit from those class.
//

#ifndef __READER_BASE__
#define __READER_BASE__

#include <fstream>
#include <string>

#include "SpringNetwork.h"
#include "topology.hpp"
#include "utils.hpp"

class ReaderBase
{
  public:
    ReaderBase() {}
    ReaderBase(const std::string & path) { setFileName(path); }
    ReaderBase(const char * const path) { setFileName(path); }
    virtual ~ReaderBase() {}
    void setFileName(const char * const path) { _filename = path; }
    void setFileName(const std::string & path) { _filename = path; }
    std::string & getFileName(void) { return _filename; }
    const std::string & getFileName(void) const { return _filename; }

    virtual void read() = 0;

    // Opens file for reading.
    // Dies if something goes wrong.
    void safeOpen() { biospring::utils::file::openread(_filename, _instream, true); }

    // Closes the file.
    void close() { _instream.close(); }

  protected:
    std::string _filename;
    std::ifstream _instream;
};

class SpnReaderBase : public ReaderBase
{
  public:
    SpnReaderBase() : ReaderBase(), _spn(0) {}
    SpnReaderBase(const std::string & path) : ReaderBase(path), _spn(0) {}
    SpnReaderBase(const char * const path) : ReaderBase(path), _spn(0) {}

    void setSpringNetwork(biospring::spn::SpringNetwork * const spn) { _spn = spn; }
    biospring::spn::SpringNetwork * getSpringNetwork(void) { return _spn; }
    const biospring::spn::SpringNetwork * const getSpringNetwork(void) const { return _spn; }

  protected:
    biospring::spn::SpringNetwork * _spn;
};

class TopologyReaderBase : public ReaderBase
{
  public:
    TopologyReaderBase() : ReaderBase() {}
    TopologyReaderBase(const std::string & path) : ReaderBase(path) {}
    TopologyReaderBase(const char * const path) : ReaderBase(path) {}

    biospring::topology::Topology & getTopology(void) { return _topology; }

  protected:
    biospring::topology::Topology _topology;
};

#endif