
#include "IO/ForceFieldReader.h"
#include "ParticleProperty.h"
#include "logging.h"
#include "utils/string.hpp"

#include <sstream>
#include <string>
#include <vector>

static float getFloat(const std::string & token, const std::string & property, const size_t & lineid)
{
    float f;
    if (not biospring::utils::string::from_string(f, token))
    {
        biospring::logging::die("ForcefieldReader: line %d: cannot convert %s token '%s' to float", lineid,
                                property.c_str(), token.c_str());
    }
    return f;
}

void ForceFieldReader::read()
{
    // Opens file for reading, dies if error occurs.
    safeOpen();

    // Initializes force field, even if was already initialized.
    _forcefield = biospring::forcefield::ForceField();

    std::string name;
    float mass = 0.0;
    float radius = 0.0;
    float charge = 0.0;
    float epsilon = 0.0;
    float transfer = 0.0;
    float hydrophobicity = 0.0;

    std::string buffer;
    size_t lineid = 0;
    while (_instream)
    {
        lineid++;
        std::getline(_instream, buffer);
        biospring::utils::string::trim(buffer);
        if (buffer[0] != '#' && !buffer.empty())
        {
            std::vector<std::string> tokens = biospring::utils::string::split(buffer);
            if (tokens.size() != 6 && tokens.size() != 7)
            {
                biospring::logging::die(
                    "ForcefieldReader: line %d: invalid number of tokens (expected 6 or 7, found %d)", lineid,
                    tokens.size());
            }

            name = tokens[0];
            charge = getFloat(tokens[1], "charge", lineid);
            radius = getFloat(tokens[2], "radius", lineid);
            epsilon = getFloat(tokens[3], "epsilon", lineid);
            mass = getFloat(tokens[4], "mass", lineid);
            transfer = getFloat(tokens[5], "transfert", lineid);

            if (tokens.size() == 7)
                hydrophobicity = getFloat(tokens[6], "hydrophobicity", lineid);

            biospring::spn::ParticleProperty pp;

            pp.setMass(mass);
            pp.setRadius(radius);
            pp.setCharge(charge);
            pp.setEpsilon(epsilon);
            pp.setTransferEnergyByAccessibleSurface(transfer);
            pp.setHydrophobicity(hydrophobicity);
            _forcefield.addPropertiesFromName(name, pp);
        }
    }
    close();
}
