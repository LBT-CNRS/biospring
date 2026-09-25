#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Particle.h"
#include "SpringNetwork.h"
#include "Vector3f.h"
#include "configuration/Configuration.hpp"
#include "forcefield/constants.hpp"
#include "forcefield/shared/random_shared.h"

using namespace biospring;

namespace
{

// Free particles: no spring, no pairwise term, nothing but the integrator and
// the bath. Whatever temperature comes out is the thermostat's doing and
// nobody else's, which is the only way to test it without arguing about a
// force field at the same time.
//
// 2000 of them because the test is statistical: the relative spread of the
// kinetic energy of N particles is sqrt(2/3N), which is 1.8 % here and would
// be 8 % at N = 100 -- a bar loose enough to pass a thermostat that is 5 %
// wrong.
void buildFreeGas(spn::SpringNetwork & network, configuration::Configuration & config, unsigned n, double temperature,
                  double timestep)
{
    for (unsigned i = 0; i < n; ++i)
    {
        spn::Particle p;
        // Spread out so nothing would interact even if a term were on.
        p.setPosition(Vector3f(10.0f * (i % 100), 10.0f * ((i / 100) % 100), 10.0f * (i / 10000)));
        p.setMass(12.0f);
        p.setDynamic(true);
        network.addParticle(p);
    }

    config = configuration::defaultConfiguration();
    config.sim.nbsteps = 4000;
    config.sim.timestep = timestep;
    config.sim.samplerate = 4000;
    config.spring.enable = false;
    config.steric.enable = false;
    config.electrostatic.enable = false;
    config.viscosity.enable = true;
    config.viscosity.value = 0.1;
    config.thermostat.enable = temperature > 0.0;
    config.thermostat.temperature = temperature;

    network.setup(config);
}

} // namespace

// The one thing a thermostat is for. Without it BioSpring's viscosity only
// ever removes energy -- measured on 072.GK, the kinetic energy falls by a
// factor of 25 over 4000 steps and keeps going -- so a run is a viscous
// relaxation and not a dynamics at any temperature.
TEST(Thermostat, HoldsTheTemperatureItIsAskedFor)
{
    spn::SpringNetwork network;
    configuration::Configuration config;
    buildFreeGas(network, config, 2000, 300.0, 1.0);
    network.run();

    // sqrt(2/3N) = 1.8 % at N = 2000, so 5 % is about three standard
    // deviations: tight enough to catch a thermostat that aims at the wrong
    // temperature, loose enough not to fail on the fluctuation it is supposed
    // to have.
    EXPECT_NEAR(network.getInstantaneousTemperature(), 300.0f, 15.0f);
}

// The same system with the thermostat off must still cool, because that is
// what the viscosity alone does and turning the thermostat off has to give
// exactly the behaviour this code had before it existed.
TEST(Thermostat, WithoutItTheSystemCools)
{
    spn::SpringNetwork network;
    configuration::Configuration config;
    buildFreeGas(network, config, 2000, 0.0, 1.0);
    // Give them something to lose: the bath is off, so only the initial
    // velocities can feed the kinetic energy.
    for (unsigned i = 0; i < network.getNumberOfParticles(); ++i)
        network.getParticle(i).setVelocity(Vector3f(0.05f, -0.03f, 0.02f));
    network.run();

    EXPECT_LT(network.getInstantaneousTemperature(), 1.0f)
        << "viscosity with no thermostat must drain the kinetic energy, not hold it";
}

// The temperature must not depend on the timestep. It does for the obvious
// implementation -- damping folded into the force, noise added on top -- which
// integrates the Langevin equation to first order only; the exact update this
// uses does not. Getting this wrong is quiet: the same .msp settles at two
// different temperatures at 2 fs and at 3 fs, and nothing says which one was
// asked for.
//
// Note what this does and does not cover. These particles feel no force, so
// the splitting is exact and what is tested is the BATH alone. With a force
// present, BAOAB's kinetic temperature carries an O(dt^2) error -- measured on
// 042.FepA it reads 268 K at 2 fs against 300 at 0.5 -- while its
// CONFIGURATIONAL averages do not, which is the whole reason to use it. Judging
// that integrator by its kinetic temperature measures the one thing it does not
// set out to improve.
TEST(Thermostat, TheTemperatureDoesNotDependOnTheTimestep)
{
    float measured[3];
    const double timesteps[3] = {2.0, 1.0, 0.5};
    for (int k = 0; k < 3; ++k)
    {
        spn::SpringNetwork network;
        configuration::Configuration config;
        buildFreeGas(network, config, 2000, 300.0, timesteps[k]);
        network.run();
        measured[k] = network.getInstantaneousTemperature();
    }

    for (int k = 0; k < 3; ++k)
        EXPECT_NEAR(measured[k], 300.0f, 15.0f) << "at timestep " << timesteps[k] << " fs";
}

// Every shipped example enables viscosity and none enables the thermostat, so
// the switch being off has to mean EXACTLY what it meant before the thermostat
// existed -- not "almost". Measured over the 46 CPU modes that use viscosity,
// 199 reported energies are identical to the bit against a run predating the
// feature; this pins the property so a later change cannot quietly break it.
//
// The trap it guards against is the obvious refactor: making the TEMPERATURE
// alone decide, so that a .msp carrying a leftover thermostat.temperature
// silently starts thermostatting. The switch decides, and nothing else.
TEST(Thermostat, AnUnaskedTemperatureChangesNothing)
{
    Vector3f cold, warm;
    for (int pass = 0; pass < 2; ++pass)
    {
        spn::SpringNetwork network;
        configuration::Configuration config;
        buildFreeGas(network, config, 200, 0.0, 1.0);
        // buildFreeGas leaves the thermostat off at temperature 0; the second
        // pass leaves it off with a temperature set, which must not matter.
        config.thermostat.enable = false;
        config.thermostat.temperature = (pass == 0) ? 0.0 : 300.0;
        network.setup(config);
        for (unsigned i = 0; i < network.getNumberOfParticles(); ++i)
            network.getParticle(i).setVelocity(Vector3f(0.05f, -0.03f, 0.02f));
        network.run();
        (pass == 0 ? cold : warm) = network.getParticle(7).getPosition();
    }

    EXPECT_FLOAT_EQ(cold.getX(), warm.getX());
    EXPECT_FLOAT_EQ(cold.getY(), warm.getY());
    EXPECT_FLOAT_EQ(cold.getZ(), warm.getZ());
}

// The generator has no state on purpose: it is called from an OpenMP loop and
// from a kernel, and the two have to draw the SAME numbers or the backends
// cannot be compared. Two properties are worth pinning: the same key gives the
// same number twice, and neighbouring keys do not.
TEST(Thermostat, TheNoiseIsReproducibleAndDecorrelated)
{
    EXPECT_FLOAT_EQ(biospring_random_normal(7, 42, 1, 0x5eed1234U),
                    biospring_random_normal(7, 42, 1, 0x5eed1234U));
    EXPECT_NE(biospring_random_normal(7, 42, 1, 0x5eed1234U), biospring_random_normal(7, 43, 1, 0x5eed1234U));
    EXPECT_NE(biospring_random_normal(7, 42, 1, 0x5eed1234U), biospring_random_normal(8, 42, 1, 0x5eed1234U));

    // Mean 0 and variance 1, over enough draws for the bar to mean something.
    double sum = 0.0, sumsq = 0.0;
    const unsigned n = 200000;
    for (unsigned i = 0; i < n; ++i)
    {
        const double x = biospring_random_normal(i / 3, i % 3, i % 3, 0x5eed1234U);
        sum += x;
        sumsq += x * x;
    }
    const double mean = sum / n;
    const double variance = sumsq / n - mean * mean;
    EXPECT_NEAR(mean, 0.0, 0.02);
    EXPECT_NEAR(variance, 1.0, 0.02);
}
