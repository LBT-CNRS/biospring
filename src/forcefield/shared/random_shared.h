#ifndef __BIOSPRING_RANDOM_SHARED_H__
#define __BIOSPRING_RANDOM_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU -- see spring_shared.h for the rules this
// file obeys and why they exist.
//
// A COUNTER-BASED generator: there is no state and no stream to advance. The
// n-th number is computed directly from (step, particle, component), so
//
//   - the CPU and the device draw exactly the SAME numbers, which is what
//     keeps the two backends comparable to the digit. A thermostat seeded from
//     two different generators would make every parity test meaningless, and
//     those tests are the only thing that has caught the real defects here;
//   - an OpenMP loop needs no per-thread state and no lock. rand() is neither
//     thread-safe nor free of correlation between threads;
//   - a run is reproducible from its seed alone, and a single particle's noise
//     at a single step can be recomputed in isolation when something has to be
//     debugged.
//
// The mixer is the "lowbias32" finalizer, applied three times over the key.
// This is not a cryptographic generator and does not need to be: what a
// thermostat asks of it is a flat spectrum and no correlation between
// neighbouring keys, which an avalanche mixer of this quality gives.

/// One round of avalanche. Every input bit reaches every output bit.
inline unsigned biospring_hash32(unsigned x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

/// Uniform bits for one (step, particle, component) triple.
/// @param seed Run seed, so two runs of the same input can differ on purpose.
inline unsigned biospring_random_bits(unsigned step, unsigned index, unsigned component, unsigned seed)
{
    // The odd multipliers are golden-ratio and Murmur constants: their only
    // job is to keep the three coordinates from colliding before the mixer
    // sees them.
    unsigned k = biospring_hash32(index + 0x9E3779B9U * component);
    k = biospring_hash32(k ^ (step + 0x85EBCA6BU));
    return biospring_hash32(k ^ seed);
}

/// Uniform in (0, 1] -- never exactly zero, because the logarithm below
/// cannot take it. 24 bits, which is every bit a float has.
inline float biospring_random_uniform(unsigned bits)
{
    return ((float)(bits >> 8) + 1.0f) * (1.0f / 16777216.0f);
}

/// A standard normal, mean 0 and variance 1, by Box-Muller.
///
/// Box-Muller produces numbers in pairs and this returns one, so half of each
/// pair is discarded. That is deliberate: keeping the other half would mean
/// carrying state between calls, which is exactly what makes a generator
/// unusable from a kernel and from an OpenMP loop. Two hashes are cheaper than
/// the state would be.
inline float biospring_random_normal(unsigned step, unsigned index, unsigned component, unsigned seed)
{
    unsigned b0 = biospring_random_bits(step, index, 2U * component, seed);
    unsigned b1 = biospring_random_bits(step, index, 2U * component + 1U, seed);
    float u1 = biospring_random_uniform(b0);
    float u2 = biospring_random_uniform(b1);
    return sqrt(-2.0f * log(u1)) * cos(6.28318530718f * u2);
}

#endif // __BIOSPRING_RANDOM_SHARED_H__
