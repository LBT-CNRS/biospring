#include "ForceField.h"
#include "constants.hpp"
#include "energy.hpp"

#include <iostream>
#include <math.h>

namespace biospring
{
namespace forcefield
{

// ======================================================================================
// Steric energy and force.
float ForceField::computeStericEnergy(float radius1, float radius2, float epsilon1, float epsilon2,
                                      float distance) const
{
    return 0.0;
}

float ForceField::computeStericForceModule(float radius1, float radius2, float epsilon1, float epsilon2,
                                           float distance) const
{
    return 0.0;
}

// ======================================================================================
// Electrostatic energy and force.
float ForceField::computeElectrostaticEnergy(float charge1, float charge2, float distance) const
{
    return _coulombscale * electrostatic_energy(charge1, charge2, distance, _dielectric);
}

float ForceField::computeElectrostaticForceModule(float charge1, float charge2, float distance) const
{
    return _coulombscale * electrostatic_force_module(charge1, charge2, distance, _dielectric);
}

// ======================================================================================
// Spring energy and force.
float ForceField::computeSpringEnergy(float distance, float stiffness, float equilibrium) const
{
    return _springscale * spring_energy(distance, stiffness, equilibrium);
}

float ForceField::computeSpringForceModule(float distance, float stiffness, float equilibrium) const
{
    return _springscale * spring_force_module(distance, stiffness, equilibrium);
}

// ======================================================================================
// IMP energy and force.

/// @brief Compute IMPALA energy (double membrane version).
/// @link https://doi.org/10.3390/membranes13030362
/// @details See @ref imp_energy function for more informations.
/// @callergraph @callgraph
/// @param x x coordinate
/// @param y y coordinate
/// @param z z coordinate
/// @param surface Solvent accessible surface of the particle
/// @param transfer Transfer energy of the particle
/// @return Return IMPALA energy of the particle in kJ.mol-1
float ForceField::computeIMPEnergy(float x, float y, float z, float surface, float transfer) const
{
    return _impscale * imp_energy(x, y, z, 
                                  surface, transfer, 
                                  _impuppermebraneoffset, _implowermembraneoffset, 
                                  _uppermembtubecurv, _lowermembtubecurv);
}

/// @brief Compute IMPALA force module (double membrane version)
/// @link https://doi.org/10.3390/ membranes13030362
/// @details See @ref imp_force_vector function for more informations.
/// @callergraph @callgraph
/// @param x x coordinate
/// @param y y coordinate
/// @param z z coordinate
/// @param surface Solvent accessible surface of the particle
/// @param transfer Transfer energy of the particle
/// @return IMPALA force vector of the particle in Da.A.fs-2
Vector3f ForceField::computeIMPForceVector(float x, float y, float z, float surface, float transfer) const
{
    return -imp_force_vector(x, y, z,
                             surface, transfer,
                             _impuppermebraneoffset, _implowermembraneoffset, 
                             _uppermembtubecurv, _lowermembtubecurv) * _impscale;
}

// ======================================================================================
// Hydrophobicity energy and force.
float ForceField::computeHydrophobicityEnergy(float hydrophobicity1, float hydrophobicity2, float distance) const
{
    return _hydrophobicityscale * hydrophobic_energy(hydrophobicity1, hydrophobicity2, distance);
}

float ForceField::computeHydrophobicityForceModule(float hydrophobicity1, float hydrophobicity2, float distance) const
{
    return _hydrophobicityscale * hydrophobic_force_module(hydrophobicity1, hydrophobicity2, distance);
}

// ======================================================================================
// Electrostatic field energy.

float ForceField::computeElectrostaticFieldEnergy(float potential, float charge) const
{
    return _forcefieldscale * electrostatic_field_energy(potential, charge);
}

// ======================================================================================
// Other methods.
void ForceField::addPropertiesFromName(const std::string name, const biospring::spn::ParticleProperty & pp)
{
    _propertiesfromname.insert(propertiesmap::value_type(name, pp));
}

biospring::spn::ParticleProperty ForceField::getPropertiesFromName(const std::string & name) const
{
    const auto it = _propertiesfromname.find(name);
    if (it != _propertiesfromname.end())
        return biospring::spn::ParticleProperty(it->second);
    return biospring::spn::ParticleProperty();
}

void ForceField::print()
{
    for (propertiesmap::const_iterator it = _propertiesfromname.begin(); it != _propertiesfromname.end(); ++it)
    {
        std::cout << it->first << " " << it->second.getCharge() << " " << it->second.getRadius() << " "
                  << it->second.getEpsilon() << " " << it->second.getMass() << " "
                  << it->second.getTransferEnergyByAccessibleSurface() << " " << std::endl;
    }
    std::cout << std::endl;
}

} // namespace forcefield
} // namespace biospring