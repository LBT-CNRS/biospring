#ifndef __TRAJECTORY_MANAGER_HPP__
#define __TRAJECTORY_MANAGER_HPP__

// Defines a class that is responsible for writing one or more trajectories to
// a file.
//
// In essence, it is a container of `WriterBase` pointers that calls each of
// their `write_step` methods in sequence.

#include <memory>
#include <vector>

#include "TrajectoryWriterBase.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

class TrajectoryManager
{
  protected:
    std::vector<std::unique_ptr<TrajectoryWriterBase>> _writers;

  public:
    // Adds a writer to the list of writers.
    void add_writer(std::unique_ptr<TrajectoryWriterBase> writer) { _writers.push_back(std::move(writer)); }

    // Writes a single step to all writers relative to writing frequency.
    void write_step(size_t frame);

    // Writes a single step to all writers directly when calling this function.
    void write_step();
};

} // namespace modern
} // namespace io
} // namespace biospring

#endif // __TRAJECTORY_MANAGER_HPP__