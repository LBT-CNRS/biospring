#ifndef __FORCEFIELD_CONSTANTS_HPP__
#define __FORCEFIELD_CONSTANTS_HPP__

// ======================================================================================
// Unit convention used throughout the forcefield module.
//
// Two different conventions coexist by design and must not be mixed:
//
// - Energies (spring, electrostatic, steric, IMP, hydrophobic, kinetic) are
//   reported in kJ.mol-1: the real, single-particle-pair interaction (computed
//   in SI units where relevant, e.g. real Coulombs and meters for
//   electrostatics) is scaled by AVOGADRO_NUMBER and converted from J to kJ,
//   purely so logged/compared energies read as familiar molar quantities.
//   User-facing "already molar" inputs (spring stiffness, steric stiffness,
//   hydrophobicity, IMP transfer energy: all in kJ.mol-1[.A-2]) skip that
//   Avogadro step since they are molar to begin with.
//
// - Forces are integrated in real, per-particle mechanical units: mass in
//   Dalton (Da), distance in Angstrom (A), time in femtosecond (fs), so
//   force = mass * acceleration must come out in Da.A.fs-2 (see
//   Particle::IntegrateVelocityVerlet / IntegrateEuler). Every force_module
//   function that starts from a molar (kJ.mol-1[.A-1|.A-2]) raw expression
//   must multiply the result by GLOBAL_SPRING_FORCE_CONVERT (or the
//   equivalent GLOBAL_IMP_FORCE_CONVERT / GLOBAL_ELECTROSTATIC_FORCE_CONVERT,
//   tailored to their own raw input units) to land in Da.A.fs-2. Forgetting
//   this conversion is what caused the hydrophobic force bug (~1e21 too
//   large): see hydrophobic_force_module in energy/hydrophobic.hpp.

namespace biospring
{
namespace forcefield
{

// ======================================================================================
// Energy conversion factors.
const double KJOULE_TO_JOULE = 1.0E3;
const double JOULE_TO_KJOULE = 1.0E-3;

const double ELECTRONVOLT_TO_JOULE = 1.6022E-19;
const double JOULE_TO_ELECTRONVOLT = 1.0 / ELECTRONVOLT_TO_JOULE;

const double ELECTRONCHARGE_TO_COULOMB = 1.60217653E-19;
const double COULOMB_TO_ELECTRONCHARGE = 1.0 / ELECTRONCHARGE_TO_COULOMB;

const double KJOULE_TO_KCAL = 0.23901;
const double KCAL_TO_KJOULE = 4.1840;

const double KCALPERMOL_TO_ELECTRONVOLPERTATOM = 4.3364E-2;

// ======================================================================================
// Distance conversion factors.
const double ANGSTROM_TO_METER = 1.0E-10;
const double METER_TO_ANGSTROM = 1.0E10;

// ======================================================================================
// Mass conversion factors.
const double DALTON_TO_UMA = 1.00794;
const double DALTON_TO_KILOGRAM = 1.66054E-27;
const double UMA_TO_KILOGRAM = 1.660538782E-27;
const double KILOGRAM_TO_DALTON = 1.0 / DALTON_TO_KILOGRAM;

// ======================================================================================
// Time conversion factors.
const double SECOND_TO_FEMTOSECOND = 1.0E15;
const double FEMTOSECOND_TO_SECOND = 1.0E-15;

// ======================================================================================
// Other constants.
const double PI = 3.1415927;
const double BOLTZMANJPERK = 1.3806504E-23;       // in J / Kelvin
const double BOLTZMANEVPERK = 8.617343E-5;        // in eV / Kelvin
const double VACUUMPERMITIVITY = 8.854187817E-12; // in F/m
const double AVOGADRO_NUMBER = 6.02214169E23;

// Coulomb's constant.
const double COULOMB_CONSTANT = 1.0 / (4.0 * PI * VACUUMPERMITIVITY);

// Converts a force in Newton (kg.m.s-2) to the simulation's mechanical unit
// system (Da.A.fs-2), matching mass in Da, distance in A and time in fs.
const double NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2 =
    METER_TO_ANGSTROM * KILOGRAM_TO_DALTON * FEMTOSECOND_TO_SECOND * FEMTOSECOND_TO_SECOND;

// ======================================================================================
// Global cutoffs.
const double MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF = 0.001;
const double MINIMAL_DISTANCE_VDW_CUTOFF = 0.1;

// ======================================================================================
// Global energy/force conversion factors.

// Converts 0.5*mass[Da]*velocity[A.fs-1]^2 (real, per-particle) into the
// molar kJ.mol-1 convention used for reported/logged kinetic energy.
const double GLOBAL_KINETIC_ENERGY_CONVERT = (ANGSTROM_TO_METER * SECOND_TO_FEMTOSECOND) *
                                             (ANGSTROM_TO_METER * SECOND_TO_FEMTOSECOND) * AVOGADRO_NUMBER *
                                             DALTON_TO_KILOGRAM * JOULE_TO_KJOULE;

// Converts a raw "stiffness[kJ.mol-1.A-2] * distance[A]" (or equivalent
// molar, kJ.mol-1.A-1-per-unit-distance) force expression into Da.A.fs-2.
// Used by spring and every steric potential (linear, Lennard-Jones
// variants): all of them take user/forcefield-supplied stiffness/epsilon
// values that are already molar (kJ.mol-1), so this factor only needs to
// strip the "per mole" convention (/AVOGADRO_NUMBER) and convert N -> Da.A.fs-2.
const double GLOBAL_SPRING_FORCE_CONVERT =
    KJOULE_TO_JOULE / AVOGADRO_NUMBER * METER_TO_ANGSTROM * NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2;

// Same role as GLOBAL_SPRING_FORCE_CONVERT, for the IMPALA membrane
// insertion force (raw expression built from ALIP/transfer energy, both in
// kJ.mol-1).
const double GLOBAL_IMP_FORCE_CONVERT =
    KJOULE_TO_JOULE / AVOGADRO_NUMBER * METER_TO_ANGSTROM * NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2;

// Converts the raw Coulomb force expression -(q1*q2)/(4.pi.dielectric.d^2),
// with q1/q2 in elementary charge (e) and d in Angstrom, directly into
// Da.A.fs-2. Unlike GLOBAL_SPRING_FORCE_CONVERT, no /AVOGADRO_NUMBER is
// needed here: charges are real per-particle quantities, not molar.
const double GLOBAL_ELECTROSTATIC_FORCE_CONVERT = NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2 *
                                                  (ELECTRONCHARGE_TO_COULOMB * ELECTRONCHARGE_TO_COULOMB) /
                                                  (VACUUMPERMITIVITY * ANGSTROM_TO_METER * ANGSTROM_TO_METER);

} // namespace forcefield
} // namespace biospring

#endif // __FORCEFIELD_CONSTANTS_HPP__
