#ifndef __NETCDFWRITER_H__
#define __NETCDFWRITER_H__

#include "WriterBase.h"

#include <netcdf>

using netCDF::NcFile;

class NetCDFWriter : public TopologyWriterBase
{
  public:
    using TopologyWriterBase::TopologyWriterBase;

    NcFile * safeOpenBinary();

    void write() override;
    void writeBinary();

  protected:
    void _writeHeaderCDL();
    void _writeParticleDataCDL();
    void _writeSpringDataCDL();

  private:
};

#endif // __NETCDFWRITER_H__
