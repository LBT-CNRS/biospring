#ifndef _CSVSAMPLEWRITER_H_
#define _CSVSAMPLEWRITER_H_

#include "WriterBase.h"

#include <fstream>
#include <string>

class CSVSampleWriter : public TrajectoryWriterBase
{
  public:
    CSVSampleWriter(const std::string & path, const biospring::spn::SpringNetwork * const spn)
        : TrajectoryWriterBase(path, spn), _firstline(true)
    {
        _timestep = 100;
    }
    CSVSampleWriter(const char * const path, const biospring::spn::SpringNetwork * const spn)
        : TrajectoryWriterBase(path, spn), _firstline(true)
    {
        _timestep = 100;
    }
    CSVSampleWriter() : CSVSampleWriter("", nullptr) {}

    void setFilename(const std::string & fn)
    {
        TrajectoryWriterBase::setFileName(fn);
        safeOpen();
    }

    void write();

  protected:
    bool _firstline;
};

#endif // _CSVSAMPLEWRITER_H_
