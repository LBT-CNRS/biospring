#ifndef __BIOSPRING_CELLGRID_SHARED_H__
#define __BIOSPRING_CELLGRID_SHARED_H__

// SHARED BETWEEN THE CPU AND THE GPU, on the same terms as spring_shared.h:
// compiled twice, once as C++ into biospring-core and once as OpenCL C
// prepended to biospring.cl. The rules that keep it compilable by both
// toolchains are spelled out there. This is not a force law -- it is the
// geometry of the neighbour search -- but it has the same reason to live here:
// the CPU and the GPU must walk the SAME cells, or their forces differ and
// nothing downstream can tell a porting bug from a physics one.
//
// WHY A CELL IS NOT A CUTOFF.
//
// Both backends used to size a cell at the cutoff and walk the 3x3x3 block
// around a particle. That block is a cube of side 3.rc; the part of it that can
// actually hold a neighbour is a sphere of radius rc. The cube is 27.rc^3 and
// the sphere is 4.19.rc^3, so 84% of the distances computed were thrown away --
// measured on 023: 2275 distances per bead per step, 353 of them inside a
// cutoff.
//
// Narrower cells fix that, because the block of cells hugs the sphere more
// closely: at w = rc/2 the stencil is 5x5x5 but its volume is 15.6 rc^3, at
// w = rc/3 it is 7x7x7 and 12.7 rc^3. More cells, fewer candidates -- and
// visiting a cell is one load, testing a candidate is a subtraction, a dot
// product and a comparison.
//
// It also decouples the cell width from the cutoff, which is what lets ONE grid
// serve several terms at once: each term takes its own stencil radius out of
// the same bins. That matters here because the cutoffs differ by a factor of
// two (steric 9 A, coulomb 16 A on 023) and a grid sized for either one is
// wrong for the other.
//
// See SpringNetwork::getCellWidth for how the width is picked.

/// How many cells of `width` a `cutoff` reaches across.
///
/// The stencil is then the cube of cells from -k to +k on each axis, which is
/// what both backends loop over, and `biospring_cell_in_range` trims its
/// corners.
///
/// Rounded UP, so the stencil always covers the cutoff sphere. The one pair it
/// can miss is at a distance exactly equal to the cutoff, when the width
/// divides the cutoff exactly; that is the same half-open convention the 3x3x3
/// stencil has always had, and it is not reachable in floating point.
///
/// Not written as ceil(cutoff / width), because that quotient is a float and
/// rounds the wrong way exactly where it matters: a width chosen to divide the
/// cutoff -- which is the interesting case, since it is the one with no slack --
/// comes out at 3.0000002 rather than 3, and the stencil gains a whole dead
/// shell (13^3 cells instead of 11^3, for nothing). So the last step is taken
/// only when it buys more than a part in a million of the cutoff, which is four
/// orders of magnitude below the 0.001 A a coordinate is written with.
inline int biospring_stencil_radius(float cutoff, float width)
{
    if (!(width > 0.0f) || !(cutoff > 0.0f))
        return 1;

    int k = (int)floor(cutoff / width);
    if ((float)k * width < cutoff * (1.0f - 1.0e-6f))
        k += 1;
    return k < 1 ? 1 : k;
}

/// The squared distance between the closest two points of a particle's own cell
/// and of the cell offset from it by (dx, dy, dz).
///
/// A particle sits somewhere in its own cell, so what bounds the pair is cell to
/// cell, not centre to centre: along one axis, a cell |d| steps away starts
/// (|d| - 1) widths from the far edge of the origin cell. Adjacent and
/// overlapping cells (|d| <= 1) touch, so they contribute nothing.
inline float biospring_cell_gap_squared(int dx, int dy, int dz, float width)
{
    const int ax = dx < 0 ? -dx : dx;
    const int ay = dy < 0 ? -dy : dy;
    const int az = dz < 0 ? -dz : dz;
    const float gx = ax > 1 ? (float)(ax - 1) * width : 0.0f;
    const float gy = ay > 1 ? (float)(ay - 1) * width : 0.0f;
    const float gz = az > 1 ? (float)(az - 1) * width : 0.0f;
    return gx * gx + gy * gy + gz * gz;
}

/// Whether the cell at that offset can hold anything within `cutoff`.
///
/// This is what pays for the finer grid: the stencil is a cube and the cutoff is
/// a sphere, so the wider the stencil the more of its corners are dead. At
/// w = 3.2 A on 023 it drops 23% of the cells before a single position is read
/// -- 1674 cells down to 1290 -- for three multiplies and a comparison each.
///
/// Conservative at the boundary (a cell exactly `cutoff` away is kept), so it
/// never removes a pair either backend's distance test would have accepted.
inline int biospring_cell_in_range(int dx, int dy, int dz, float width, float cutoffsquared)
{
    return biospring_cell_gap_squared(dx, dy, dz, width) <= cutoffsquared;
}

#endif // __BIOSPRING_CELLGRID_SHARED_H__
