#ifndef _PQRWRITER_H_
#define _PQRWRITER_H_

#include "PDBWriter.h"

class PQRWriter : public PDBWriter
{
  public:
    using PDBWriter::PDBWriter;
    void writeModel(size_t modelid) override;
};

#endif
