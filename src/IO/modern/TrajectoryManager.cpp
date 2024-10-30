#include "TrajectoryManager.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

void TrajectoryManager::write_step(size_t frame)
{
    for (const auto & writer : _writers)
    {
        if (frame % writer->write_frequency() == 0)
            writer->write_step();
    }
}

void TrajectoryManager::write_step()
{
    for (const auto & writer : _writers)
    {
        writer->write_step();
    }
}

} // namespace modern
} // namespace io
} // namespace biospring