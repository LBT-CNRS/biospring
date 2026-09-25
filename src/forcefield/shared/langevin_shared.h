#ifndef __BIOSPRING_LANGEVIN_SHARED_H__
#define __BIOSPRING_LANGEVIN_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU -- see spring_shared.h for the rules this
// file obeys and why they exist.
//
// The solvent this model does not simulate did two things to a bead: it slowed
// down whatever moved, and it kicked whatever stood still. BioSpring kept the
// first -- viscosity.value -- and dropped the second, so energy could only ever
// leave and the structure cooled towards a standstill. Measured on 072: the
// kinetic energy falls from 15.93 to 0.63 kJ.mol-1 over 4000 steps, where 300 K
// on that system is around 10500.
//
// The two are not independent knobs. They are one physical effect, and the
// fluctuation-dissipation relation fixes the size of the kick from the friction
// and the temperature. So the same gamma now does both.
//
// The velocity update is the EXACT solution of dv = -(gamma/m) v dt + noise
// over one step, not a first-order approximation of it:
//
//     v <- c v + sigma xi,   c = exp(-gamma dt / m),   sigma^2 = kB T (1-c^2)/m
//
// Exact matters because the approximation's stationary temperature drifts with
// the timestep: the same .msp at 2 fs and at 3 fs would settle at two different
// temperatures, and nothing would say which one was asked for.

/// How much of the velocity survives one step of friction.
/// @param gamma    viscosity.value, in Da.fs-1.
/// @param mass     Particle mass in Da. A massless particle is not integrated
///                 and is not thermostatted either, hence the 1.
/// @param timestep In fs.
/// @return The dimensionless factor c in (0, 1].
inline float biospring_langevin_decay(float gamma, float mass, float timestep)
{
    return mass > 0.0f ? exp(-gamma * timestep / mass) : 1.0f;
}

/// The matching kick, so that the two together hold the bath temperature.
/// @param decay The c above -- passed rather than recomputed, because the two
///     must be the same number to the bit or the temperature is not the one
///     that was asked for.
/// @param mass  Particle mass in Da.
/// @param boltzmanntemperature kB*T in the integrator's own units
///     (Da.A^2.fs-2), which the host computes once. Zero is a legitimate
///     value: it leaves the decay alone and the thermostat is then the pure
///     brake this code had before.
/// @return The standard deviation of one velocity component's kick, in A.fs-1.
inline float biospring_langevin_kick(float decay, float mass, float boltzmanntemperature)
{
    if (mass <= 0.0f || boltzmanntemperature <= 0.0f)
        return 0.0f;
    return sqrt(boltzmanntemperature / mass * (1.0f - decay * decay));
}

#endif // __BIOSPRING_LANGEVIN_SHARED_H__
