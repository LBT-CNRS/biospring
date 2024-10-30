#ifndef __BIOSPRING_IO_H__
#define __BIOSPRING_IO_H__

#include "SpringNetwork.h"
#include "topology.hpp"

#include <string>
#include <vector>

namespace biospring
{

namespace io
{
topology::Topology readTopology(const std::string & path);

void writeTopology(const std::string & path, const topology::Topology & topology);
void writeTopology(const vector<string> & outputfiles, const topology::Topology & topology);

} // namespace io

} // namespace biospring

#endif // __BIOSPRING_IO_H__