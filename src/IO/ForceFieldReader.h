#ifndef _FORCEREADER_H_
#define _FORCEREADER_H_

#include "forcefield/ForceField.h"
#include "IO/ReaderBase.h"

class ForceFieldReader : public ReaderBase
{
  public:
    ForceFieldReader() : ReaderBase() {}
    ForceFieldReader(const std::string & path) : ReaderBase(path) {}
    ForceFieldReader(const char * const path) : ReaderBase(path) {}

    const biospring::forcefield::ForceField & getForceField() const { return _forcefield; }

    void read();

  private:
    biospring::forcefield::ForceField _forcefield;
};

#endif
