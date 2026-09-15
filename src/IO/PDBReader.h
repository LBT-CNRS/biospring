#ifndef _PDBREADER_H_
#define _PDBREADER_H_

#include "IO/ReaderBase.h"
#include "topology.hpp"

#include <map>
#include <string>
#include <vector>

#define UNDEFINEDINDEX -1

class PDBReader : public TopologyReaderBase
{
  public:
    PDBReader() : TopologyReaderBase() {}
    PDBReader(const std::string & path) : TopologyReaderBase(path) {}
    PDBReader(const char * const path) : TopologyReaderBase(path) {}

    static biospring::topology::Particle parseAtomLine(const std::string & line);

    // Returns true if a line starts with ATOM or HETATM.
    static bool isAtomLine(const std::string & line)
    {
        std::string record(line.substr(0, 6));
        return (record == "ATOM  " or record == "HETATM");
    }

    static std::vector<std::pair<size_t, size_t>> parseConectLine(const std::string & line);

    // Returns true if a line starts with CONECT.
    static bool isConectLine(const std::string & line)
    {
        std::string record(line.substr(0, 6));
        return (record == "CONECT");
    }

    void read();
    void updatePositions();
    int getIdFromExtid(size_t extid) const;

    void addAtomfilter(const std::string & filter) { _atomfilter.push_back(filter); }
    bool isInFilter(const std::string & filter) const;
    bool isInFilter(const char * const filter) const { return isInFilter(std::string(filter)); }

  protected:
    std::map<int, size_t> _extidtoindex;
    std::vector<std::string> _atomfilter;
};

#endif
