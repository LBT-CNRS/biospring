#include "IO/LobeSpec.h"

#include "Particle.h"
#include "logging.h"
#include "utils/string.hpp"

namespace biospring
{
namespace io
{

std::string parse_lobe_spec(const std::string & token, int & mode, float & angle, const char * reader,
                            size_t line_id)
{
    mode = spn::Particle::HBOND_LOBES_NONE;
    angle = 62.0f;
    if (token.empty())
        logging::die("%s: line %d: empty plane or antecedent column", reader, static_cast<int>(line_id));

    std::string atom = token;
    const char marker = atom[0];
    if (marker != '~' && marker != '^')
        return atom; // a plain second antecedent, no lobes

    mode = marker == '~' ? spn::Particle::HBOND_LOBES_IN_PLANE : spn::Particle::HBOND_LOBES_OUT_OF_PLANE;
    atom = atom.substr(1);

    const size_t colon = atom.find(':');
    if (colon != std::string::npos)
    {
        if (!utils::string::from_string(angle, atom.substr(colon + 1)))
            logging::die("%s: line %d: invalid lobe angle '%s'", reader, static_cast<int>(line_id),
                         atom.substr(colon + 1).c_str());
        atom = atom.substr(0, colon);
    }
    else if (marker == '^')
        logging::die("%s: line %d: '^%s' needs an explicit angle, as '^%s:105' -- an out-of-plane "
                     "approach has no default",
                     reader, static_cast<int>(line_id), atom.c_str(), atom.c_str());

    if (atom.empty())
        logging::die("%s: line %d: '%c' names no plane atom", reader, static_cast<int>(line_id), marker);
    if (angle <= 0.0f || angle >= 180.0f)
        logging::die("%s: line %d: lobe angle %g is outside (0, 180)", reader, static_cast<int>(line_id),
                     static_cast<double>(angle));
    return atom;
}

} // namespace io
} // namespace biospring
