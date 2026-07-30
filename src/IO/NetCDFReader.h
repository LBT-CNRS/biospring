#ifndef _NETCDFREADER_H_
#define _NETCDFREADER_H_

#include "IO/ReaderBase.h"
#include "IO/SpnBuffer.h"

#include <memory>
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
    DihedralSpringBuffer _dihedralbackbonebuffer;
    DihedralSpringBuffer _dihedralsidechainbuffer;
    DihedralSpringBuffer _dihedralplanaritybuffer;

    std::unique_ptr<netCDF::NcFile> _file;
    int _filec = -1;

    void readSprings();
    void readNumberOfSprings();

    // Reads one dihedral ghost-spring family (optional -- old .nc files, or
    // families that were empty at write time, simply have none of these
    // variables; see readDihedralSpringGroup's own dimension check).
    void readDihedralSpringGroup(const char * prefix, DihedralSpringBuffer & buffer);

    void readParticles();
    void readNumberOfParticles();

    void addParticlesToSpn();
    void addSpringsToSpn();
    void addDihedralBackboneSpringsToSpn();
    void addDihedralSidechainSpringsToSpn();
    void addDihedralPlanaritySpringsToSpn();

    void checkNDims(const netCDF::NcVar & var, int ref);
    void checkDim(const netCDF::NcVar & var, int dimid, size_t size);
    netCDF::NcVar getNcVar(const char * varname, bool mandatory = true);
};

#endif
