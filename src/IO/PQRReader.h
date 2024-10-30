
#ifndef __PQRREADER_H__
#define __PQRREADER_H__

#include "PDBReader.h"
#include "topology.hpp"

class PQRReader : public PDBReader
{
  public:
    void read();
    static biospring::topology::Particle parseAtomLine(const std::string & line);
};

#endif
