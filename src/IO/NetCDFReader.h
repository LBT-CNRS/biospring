#ifndef _NETCDFREADER_H_
#define _NETCDFREADER_H_

#include "IO/ReaderBase.h"
#include "IO/SpnBuffer.h"

#include <netcdf>

class NetCDFReader : public TopologyReaderBase
{
  public:
    NetCDFReader() : TopologyReaderBase() {}
    NetCDFReader(const std::string & path) : TopologyReaderBase(path) {}
    NetCDFReader(const char * const path) : TopologyReaderBase(path) {}

    void read();

  protected:
    SpringBuffer _sbuffer;
    ParticleBuffer _pbuffer;

    netCDF::NcFile * _file;
    int _filec;

    void readSprings();
    void readNumberOfSprings();

    void readParticles();
    void readNumberOfParticles();

    void addParticlesToSpn();
    void addSpringsToSpn();

    void checkNDims(const netCDF::NcVar & var, int ref);
    void checkDim(const netCDF::NcVar & var, int dimid, size_t size);
    netCDF::NcVar getNcVar(const char * varname, bool mandatory = true);
};

#endif
