#ifndef _NETCDFREADER_H_
#define _NETCDFREADER_H_

#include "IO/ReaderBase.h"
#include "IO/SpnBuffer.h"

#include <array>
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
    // Reuses DihedralSpringBuffer's layout (see NetCDFWriter's stretch/bend
    // call sites): structurally identical to a dihedral ghost-spring group,
    // just with springsdcoffsets always zero.
    DihedralSpringBuffer _stretchbuffer;
    DihedralSpringBuffer _bendbuffer;
    // One buffer per dihedral family, indexed by the family itself (see
    // spn::SpringNetwork::DihedralFamilyIndex) -- read back under the .nc
    // variable prefix that same index names, so reader and writer walk the
    // families in one shared order instead of two hand-kept lists.
    std::array<DihedralSpringBuffer, biospring::spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT> _dihedralbuffers;
    GhostParticleBuffer _ghostparticlebuffer;
    DihedralAxisBuffer _dihedralaxisbuffer;

    std::unique_ptr<netCDF::NcFile> _file;
    int _filec = -1;

    void readSprings();
    void readNumberOfSprings();

    // Reads one dihedral ghost-spring family (optional -- old .nc files, or
    // families that were empty at write time, simply have none of these
    // variables; see readDihedralSpringGroup's own dimension check).
    void readDihedralSpringGroup(const char * prefix, DihedralSpringBuffer & buffer);

    // Reads ghost-particle anchor bindings (optional, same convention as
    // readDihedralSpringGroup: an old .nc file, or one with no ghost
    // particles, simply has none of these variables).
    void readGhostParticles();
    void readDihedralAxes();

    void readParticles();
    void readNumberOfParticles();

    void addParticlesToSpn();
    void addSpringsToSpn();
    void addStretchSpringsToSpn();
    void addBendSpringsToSpn();
    void addDihedralSpringsToSpn(unsigned family);
    void addGhostParticlesToSpn();
    void addDihedralAxesToSpn();

    void checkNDims(const netCDF::NcVar & var, int ref);
    void checkDim(const netCDF::NcVar & var, int dimid, size_t size);
    netCDF::NcVar getNcVar(const char * varname, bool mandatory = true);
};

#endif
