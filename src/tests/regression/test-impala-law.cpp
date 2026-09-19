#include <gtest/gtest.h>

#include <cmath>

#include "Vector3f.h"
#include "forcefield/constants.hpp"
#include "forcefield/energy/imp.hpp"
#include "forcefield/shared/imp_shared.h"

using namespace biospring::forcefield;

// The OpenCL backend evaluates IMPALA from shared/imp_shared.h and the CPU from
// energy/imp.hpp. The two must agree on the flat single membrane, which is the
// only case the device implements and the only one any batch run takes: the
// four geometry parameters default to 0 and no .msp key reaches them, only an
// MDDriver client can set them.
//
// This is what holds the two backends together. Nothing else compares them
// outside the OpenCL parity tests, which need a device.
TEST(ImpalaLaw, SharedHeaderMatchesTheCPUOnAFlatMembrane)
{
    double worstenergy = 0.0, worstforce = 0.0;
    unsigned points = 0;

    // Out to 400 A, far past the 60 A where exp(alpha*(|z| - z0)) overflows a
    // float. That is not a hypothetical: example 051 reaches |z| = 76 A and the
    // first version of the shared header returned a NaN there, which killed the
    // run with "Found non-finite position for particle 1".
    for (float z = -400.0f; z <= 400.0f; z += 0.37f)
        for (float surface : {0.0f, 12.5f, 91.2f, 200.0f})
            for (float transfer : {-0.4393f, -0.0561f, 0.0f, 0.1686f, 0.4686f})
            {
                // x and y deliberately non-zero: a flat membrane has no lateral
                // gradient and the force must not acquire one.
                const float e_cpu = imp_energy(1.3f, -2.7f, z, surface, transfer, 0.0f, 0.0f, 0.0f, 0.0f);
                const float e_shared = biospring_imp_energy(z, surface, transfer, ALIP, ALPHA, Z0);

                const Vector3f f_cpu = imp_force_vector(1.3f, -2.7f, z, surface, transfer, 0.0f, 0.0f, 0.0f, 0.0f);
                // ForceField::computeIMPForceVector returns the NEGATIVE of
                // imp_force_vector, and the shared helper carries that sign.
                const float f_shared = biospring_imp_force_z(z, surface, transfer, ALIP, ALPHA, Z0,
                                                             static_cast<float>(GLOBAL_IMP_FORCE_CONVERT));

                EXPECT_TRUE(std::isfinite(e_shared)) << "energy is not finite at z = " << z;
                EXPECT_TRUE(std::isfinite(f_shared)) << "force is not finite at z = " << z;
                EXPECT_FLOAT_EQ(f_cpu.getX(), 0.0f) << "a flat membrane pushed sideways at z = " << z;
                EXPECT_FLOAT_EQ(f_cpu.getY(), 0.0f) << "a flat membrane pushed sideways at z = " << z;

                worstenergy = std::max<double>(worstenergy, std::fabs(e_cpu - e_shared));
                worstforce = std::max<double>(worstforce, std::fabs(-f_cpu.getZ() - f_shared));
                points++;
            }

    // Measured 7.6e-06 and 9.3e-10 over these 43260 points. The energy is of
    // order 20 kJ.mol-1 there, so that is single-precision rounding; the bars
    // are two orders above and still far below anything a wrong sign, a missing
    // factor or a different profile would cause.
    EXPECT_GT(points, 40000u) << "the sweep did not cover what it claims to";
    EXPECT_LT(worstenergy, 1.0e-3) << "energy disagrees by " << worstenergy << " kJ.mol-1";
    EXPECT_LT(worstforce, 1.0e-7) << "force disagrees by " << worstforce << " Da.A.fs-2";
}

// A membrane whose offset and curvature are both zero is an absent membrane,
// not one sitting at z = 0. It used to be the latter: the branch test asked
// whether ALL FOUR parameters were zero, so nudging one of them switched on a
// second membrane that coincided with the first, and the energy and the force
// both jumped by a factor of exactly two -- measured at an offset of 1e-30 A,
// long before the two are physically distinct.
TEST(ImpalaLaw, NudgingOneOffsetDoesNotDoubleTheEnergy)
{
    const float surface = 91.2f, transfer = -0.4393f;
    const float epsilon = 1.0e-30f; // physically zero, but the condition sees it

    for (float z : {2.0f, 8.0f, 14.0f, 20.0f, 30.0f})
    {
        const float flat = imp_energy(0.0f, 0.0f, z, surface, transfer, 0.0f, 0.0f, 0.0f, 0.0f);
        const float nudged = imp_energy(0.0f, 0.0f, z, surface, transfer, epsilon, 0.0f, 0.0f, 0.0f);
        EXPECT_NEAR(nudged, flat, 1.0e-4f) << "energy jumps at z = " << z;

        const float fflat = imp_force_vector(0.0f, 0.0f, z, surface, transfer, 0.0f, 0.0f, 0.0f, 0.0f).getZ();
        const float fnudged = imp_force_vector(0.0f, 0.0f, z, surface, transfer, epsilon, 0.0f, 0.0f, 0.0f).getZ();
        EXPECT_NEAR(fnudged, fflat, std::fabs(fflat) * 1.0e-4f + 1.0e-20f) << "force jumps at z = " << z;
    }

    // And two membranes really set apart still give two contributions.
    const float one = imp_energy(0.0f, 0.0f, 4.0f, surface, transfer, 20.0f, 0.0f, 0.0f, 0.0f);
    const float two = imp_energy(0.0f, 0.0f, 4.0f, surface, transfer, 20.0f, 20.0f, 0.0f, 0.0f);
    EXPECT_NE(one, two) << "the lower membrane stopped contributing when it was actually requested";
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
