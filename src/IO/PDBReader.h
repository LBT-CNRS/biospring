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

    // SSBOND is the PDB's OWN record for a disulfide bridge and LINK its own
    // record for any other bond between two residues; both were ignored, so a
    // file declaring its bridges the standard way -- GKinase.1S4Q.pdb does --
    // came in with no bond at all, and only a hand-added CONECT worked.
    //
    // Unlike CONECT they name atoms by residue rather than by serial, and they
    // sit in the header, BEFORE the atoms they refer to. So they are collected
    // while reading and resolved once every atom is known.
    struct ResidueAtomRef
    {
        std::string chain;
        int residue_id = 0;
        std::string name;
    };

    static bool isSSBondLine(const std::string & line) { return line.compare(0, 6, "SSBOND") == 0; }
    static bool isLinkLine(const std::string & line) { return line.compare(0, 6, "LINK  ") == 0; }

    // Both return false on a line too short or malformed to read, which is
    // warned about and skipped rather than thrown.
    static bool parseSSBondLine(const std::string & line, ResidueAtomRef & first, ResidueAtomRef & second);
    static bool parseLinkLine(const std::string & line, ResidueAtomRef & first, ResidueAtomRef & second);

    void read();
    void updatePositions();
    int getIdFromExtid(size_t extid) const;

    void addAtomfilter(const std::string & filter) { _atomfilter.push_back(filter); }
    bool isInFilter(const std::string & filter) const;
    bool isInFilter(const char * const filter) const { return isInFilter(std::string(filter)); }

  protected:
    // Adds one spring per collected SSBOND/LINK, once every atom is known.
    void _resolveResidueBonds();

    std::map<int, size_t> _extidtoindex;
    std::vector<std::string> _atomfilter;
    std::vector<std::pair<ResidueAtomRef, ResidueAtomRef>> _residue_bonds;
};

#endif
