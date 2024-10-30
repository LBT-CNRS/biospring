#ifndef _OPENDXREADER_H_
#define _OPENDXREADER_H_

#include "IO/ReaderBase.h"
#include "grid/PotentialGrid.hpp"

#include <array>
#include <fstream>
#include <string>

class OpenDXReader : public ReaderBase
{
  public:
    OpenDXReader() : ReaderBase() {}
    OpenDXReader(const std::string & path) : ReaderBase(path) {}
    OpenDXReader(const char * const path) : ReaderBase(path) {}
    ~OpenDXReader() {}

    void read();

    biospring::grid::PotentialGrid & getGrid() { return _grid; }
    const biospring::grid::PotentialGrid & getGrid() const { return _grid; }

  protected:
    std::array<size_t, 3> readSize();
    std::array<double, 3> readOrigin();
    std::array<double, 3> readScalingFactors();
    void readGrid();

  private:
    biospring::grid::PotentialGrid _grid;
};

namespace biospring
{
namespace opendx
{

inline biospring::grid::PotentialGrid readGrid(const std::string & path)
{
    OpenDXReader reader(path);
    reader.read();
    return reader.getGrid();
}

} // namespace opendx
} // namespace biospring

#endif
