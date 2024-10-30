#include "PDBTrajectoryWriter.hpp"
#include "pdb_format.hpp"

namespace biospring
{
namespace io
{
namespace modern
{

void PDBTrajectoryWriter::write_step()
{
    // Writes MODEL record.
    _ostream << pdbfmt::model_header(_current_frame) << std::endl;

    // Writes ATOM records.
    for (const auto & particle : _topology.getParticles())
    {
        _ostream << pdbfmt::atom_record(particle) << std::endl;
    }

    // Writes CONECT records for the first step only.
    if (_current_frame == 0)
        for (const auto & spring : _topology.getSprings())
            _ostream << pdbfmt::conect_record(spring) << std::endl;

    // Writes ENDMDL record.
    _ostream << pdbfmt::model_footer() << std::endl;

    // Increments step counter.
    _current_frame++;
}

} // namespace modern
} // namespace io
} // namespace biospring
