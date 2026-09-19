#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

typedef struct{
	unsigned id1,id2;
	float equilibrium;
	float stiffness;
	} Springocl;


// ============================================================================
// Cell list
// ============================================================================
//
// Every non-bonded term needs the same thing: given a particle, the handful of
// particles close enough to matter. The CPU answers it with an infinite grid of
// cells (nsearch.hpp); this is the device's version of the same structure.
//
// ONE GRID PER TERM, not one shared. The cell width IS the cutoff, and the
// three cutoffs differ -- steric 8 A, electrostatic 16, hydrophobicity 15 by
// default -- so a shared grid at the longest of them would make the shortest
// term walk the longest one's volume, eight times the candidates for the same
// answer. The kernels below take the frame as arguments for exactly that
// reason: the same code bins into whichever grid it is handed.
//
// Built as a linked list per cell rather than a sorted array, after the method
// in Marcus Bannerman's OpenCL course (exercise 3, "sorting particles"): each
// particle atomically exchanges itself into its cell's head and keeps whatever
// was there as its own successor. One atomic per particle, no counting pass, no
// sort, and the build is O(N) with nothing to size in advance.
//
//   cellhead[c] -- the last particle binned into cell c, or EMPTY
//   nextincell[p] -- the particle binned into p's cell just before it, or EMPTY
//
// Two things the course's version does not need and this one does. Its
// particles live in [0,1) with a periodic box; a protein sits wherever the PDB
// put it, so cells are indexed from an origin the host recomputes, and a
// stencil that falls outside the grid is skipped instead of wrapped.

#define BIOSPRING_EMPTY_CELL ((uint)-1)

// Which cell a position falls in. Returns the flat index, or BIOSPRING_EMPTY_CELL
// when the position is outside the grid -- which happens between two rebuilds,
// since a particle keeps moving after the box was measured.
inline uint biospring_cell_of(const float4 position, const float4 origin,
                              const float cellwidth, const int4 ncells)
	{
	int x = (int)floor((position.x - origin.x) / cellwidth);
	int y = (int)floor((position.y - origin.y) / cellwidth);
	int z = (int)floor((position.z - origin.z) / cellwidth);

	if (x < 0 || y < 0 || z < 0 || x >= ncells.x || y >= ncells.y || z >= ncells.z)
		return BIOSPRING_EMPTY_CELL;

	return (uint)((z * ncells.y + y) * ncells.x + x);
	}


__kernel void blankCells(const uint ncellstotal, __global uint * cellhead)
	{
	const uint c = get_global_id(0);
	if (c < ncellstotal)
		cellhead[c] = BIOSPRING_EMPTY_CELL;
	}


__kernel void binParticles(const __global float4 * positions,
                           const float4 origin, const float cellwidth, const int4 ncells,
                           __global uint * cellhead, __global uint * nextincell,
                           const uint N)
	{
	const uint p = get_global_id(0);
	if (p >= N) return;

	const uint c = biospring_cell_of(positions[p], origin, cellwidth, ncells);
	if (c == BIOSPRING_EMPTY_CELL)
		{
		// Outside the frame: linked to nothing, and no cell points at it, so it
		// is invisible to every neighbour walk. The host checks for this before
		// binning -- see _frameStillHolds -- and remeasures rather than letting
		// it happen; reaching here at all means the frame was already stale.
		nextincell[p] = BIOSPRING_EMPTY_CELL;
		return;
		}

	// Atomic: several particles land in the same cell in the same instant, and
	// a plain read-modify-write loses all but one of them.
	nextincell[p] = atom_xchg(cellhead + c, p);
	}

// The interactive probe: one extra particle that every other one feels, and
// that feels every other one back.
//
// It has no cutoff and no cell list -- it is an all-pairs term with a single
// partner, so one work item per particle is the whole parallelisation. The
// laws are the ordinary steric and Coulomb ones, from the same shared headers
// the pairwise kernels use; what differs is only that one side of every pair
// is the same particle.
//
// THE FORCE GOES BOTH WAYS. Particle::addStericProbeForce adds +f to the
// particle and -f to the probe, which is why the CPU runs this serially: a
// parallel loop would race on the probe. Here each work item writes its own -f
// into `probeforces[tid]` and the host sums that, which is a reduction it was
// going to pay for anyway -- the probe is integrated on the host, as one
// particle among N it is not worth a kernel for.
//
// Only DYNAMIC particles take part, matching the CPU's loop over
// _dynamicparticules, and the probe is skipped against itself.
__kernel void probe(const __global float4 * positions,
                    const __global float * radii,
                    const __global float * epsilons,
                    const __global float * charges,
                    const __global int * dynamicstate,
                    __global float4 * forces,
                    __global float4 * probeforces,
                    const uint probeid,
                    const int stericenabled, const int coulombenabled,
                    const int mode, const float proberadius, const float probeepsilon,
                    const float probecharge, const float dielectric,
                    const float linearstiffness, const float stericmindistance,
                    const float coulombmindistance, const float fourpi,
                    const float stericconvert, const float coulombconvert,
                    const float stericscale, const float coulombscale,
                    const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	probeforces[tid] = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
	if (tid == probeid || !dynamicstate[tid]) return;

	const float4 axis = positions[probeid] - positions[tid];
	const float distsq = axis.x * axis.x + axis.y * axis.y + axis.z * axis.z;
	if (distsq == 0.0f) return;
	const float dist = sqrt(distsq);

	float3 f = (float3)(0.0f, 0.0f, 0.0f);

	if (stericenabled)
		{
		const float module = biospring_steric_force_module(
		    mode, proberadius, radii[tid], probeepsilon, epsilons[tid], dist,
		    linearstiffness, stericmindistance, stericconvert);
		f += (axis.xyz / dist) * (stericscale * module);
		}

	if (coulombenabled)
		{
		const float module = biospring_electrostatic_force_module(
		    probecharge, charges[tid], dist, dielectric, coulombmindistance, fourpi,
		    coulombconvert);
		f += (axis.xyz / dist) * (coulombscale * module);
		}

	forces[tid].xyz += f;
	probeforces[tid].xyz = -f;
	}

// Sums what the probe kernel scattered and gives the probe its own force.
//
// A reduction rather than an atomic, because OpenCL 1.2 has no atomic add on
// floats. ONE work group: the result is a single float4, so a second pass would
// cost more than it saves. Each lane walks the array with a stride of the group
// size -- coalesced -- then the group folds its partial sums in local memory.
//
// Launched with exactly WORK_GROUP_SIZE work items, and `partial` sized to
// match.
__kernel void probegather(const __global float4 * probeforces,
                          __global float4 * forces,
                          __local float4 * partial,
                          const uint probeid,
                          const uint N)
	{
	const uint lid = get_local_id(0);
	const uint groupsize = get_local_size(0);

	float4 sum = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
	for (uint i = lid; i < N; i += groupsize)
		sum += probeforces[i];
	partial[lid] = sum;

	barrier(CLK_LOCAL_MEM_FENCE);
	for (uint stride = groupsize / 2; stride > 0; stride >>= 1)
		{
		if (lid < stride)
			partial[lid] += partial[lid + stride];
		barrier(CLK_LOCAL_MEM_FENCE);
		}

	if (lid == 0)
		forces[probeid].xyz += partial[0].xyz;
	}

__kernel void spring(const __global float4 * positions,
                     const __global Springocl * springs,
                     const __global int * springoffsets,
                     __global float4 * forces,
                     __global float * energy,
                     const uint N,
                     const float scale)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) return;

	// CSR offsets: where this particle's springs start, and where they stop.
	// A particle with no spring has begin == end and the loop does not run.
	int begin = springoffsets[tid];
	int end = springoffsets[tid+1];

	float3 here = positions[tid].xyz;
	float3 sum = (float3)(0.0f, 0.0f, 0.0f);
	float e = 0.0f;

	for(int i=begin;i<end;i++)
		{
		float3 axis = positions[springs[i].id2].xyz - here;
		float dist = length(axis);
		sum += normalize(axis) * biospring_spring_force_module(dist,
		                                                      springs[i].stiffness,
		                                                      springs[i].equilibrium,
		                                                      scale);
		// Half, because the CSR holds each spring from both of its ends and
		// this kernel runs once per end. Summed over every particle, the halves
		// make each spring's energy exactly once.
		e += 0.5f * biospring_spring_energy(dist, springs[i].stiffness, springs[i].equilibrium);
		}

	forces[tid].xyz += sum;
	energy[tid] = e;
	}


// Is `other` joined to `self` by a spring?
//
// The exclusion the CPU makes in Particle::addElectrostaticForce, guarded by
// spring.enable there and by `springsenabled` here. It is not optional: the
// mesh puts bonded atoms 1 to 1.5 A apart with opposite partial charges, and
// a Coulomb term that sees them pulls with hundreds of kJ.mol-1.A-1 on pairs
// a spring is already holding.
//
// The spring CSR is already on the device for the spring kernel, so this costs
// a walk over the four to ten springs a particle has.
inline bool biospring_sprung_together(const __global Springocl * springs,
                                      const __global int * springoffsets,
                                      const uint self, const uint other)
	{
	const int begin = springoffsets[self];
	const int end = springoffsets[self + 1];
	for (int i = begin; i < end; i++)
		if (springs[i].id2 == other)
			return true;
	return false;
	}


// Coulomb, gathered over the cell list.
//
// Work item tid owns particle tid and is the only writer of forces[tid], so no
// atomics and no deferred second pass. Each pair is therefore evaluated twice,
// once from each end -- the CPU evaluates it once and hands the other half over
// (see the deferred scratch in Particle::addElectrostaticForce), which is the
// right trade on a core and the wrong one on a device, where the arithmetic is
// cheaper than the bookkeeping.
//
// `cellwidth` is the grid's, which may be WIDER than `cutoff` (see
// _measureCellGrid): the stencil then covers more than asked and the distance
// test below drops the surplus. Reading the cutoff from the grid instead would
// silently extend the term's range.
__kernel void electrostatic(const __global float4 * positions,
                            const __global float * charges,
                            __global float4 * forces,
                            const __global uint * cellhead,
                            const __global uint * nextincell,
                            const float4 origin, const float cellwidth, const int4 ncells,
                            const __global Springocl * springs,
                            const __global int * springoffsets,
                            const int springsenabled,
                            const float cutoff, const float dielectric, const float mindistance,
                            const float fourpi, const float convert, const float coulombscale,
                            const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 here = positions[tid];
	const float q = charges[tid];
	const float cutoffsq = cutoff * cutoff;

	const int cx = (int)floor((here.x - origin.x) / cellwidth);
	const int cy = (int)floor((here.y - origin.y) / cellwidth);
	const int cz = (int)floor((here.z - origin.z) / cellwidth);

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	for (int dz = -1; dz <= 1; dz++)
		for (int dy = -1; dy <= 1; dy++)
			for (int dx = -1; dx <= 1; dx++)
				{
				const int x = cx + dx, y = cy + dy, z = cz + dz;
				// Nothing is periodic here: a stencil cell outside the grid is
				// absent, never the cell on the opposite face.
				if (x < 0 || y < 0 || z < 0 || x >= ncells.x || y >= ncells.y || z >= ncells.z)
					continue;

				const uint cell = (uint)((z * ncells.y + y) * ncells.x + x);
				for (uint p = cellhead[cell]; p != BIOSPRING_EMPTY_CELL; p = nextincell[p])
					{
					if (p == tid)
						continue;

					float3 axis = positions[p].xyz - here.xyz;
					const float distsq = axis.x * axis.x + axis.y * axis.y + axis.z * axis.z;
					if (distsq > cutoffsq || distsq == 0.0f)
						continue;

					if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
						continue;

					const float dist = sqrt(distsq);
					const float module = biospring_electrostatic_force_module(
					    charges[p], q, dist, dielectric, mindistance, fourpi, convert);
					sum += (axis / dist) * (coulombscale * module);
					}
				}

	forces[tid].xyz += sum;
	}


// Steric, gathered over its own cell list.
//
// Same shape as the electrostatic kernel above, and the same reasons: one
// writer per particle so no atomics, each pair evaluated from both ends, the
// spring exclusion walked over the CSR already on the device, and `cellwidth`
// taken from the grid while `cutoff` comes from the term.
//
// `mode` selects among the four laws the .msp can ask for. One kernel rather
// than four builds: the branch is uniform across the whole grid -- every work
// item takes the same one -- so it costs nothing a separate build would save.
__kernel void steric(const __global float4 * positions,
                     const __global float * radii,
                     const __global float * epsilons,
                     __global float4 * forces,
                     const __global uint * cellhead,
                     const __global uint * nextincell,
                     const float4 origin, const float cellwidth, const int4 ncells,
                     const __global Springocl * springs,
                     const __global int * springoffsets,
                     const int springsenabled,
                     const int mode,
                     const float cutoff, const float linearstiffness, const float mindistance,
                     const float convert, const float stericscale,
                     const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 here = positions[tid];
	const float radius = radii[tid];
	const float epsilon = epsilons[tid];
	const float cutoffsq = cutoff * cutoff;

	const int cx = (int)floor((here.x - origin.x) / cellwidth);
	const int cy = (int)floor((here.y - origin.y) / cellwidth);
	const int cz = (int)floor((here.z - origin.z) / cellwidth);

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	for (int dz = -1; dz <= 1; dz++)
		for (int dy = -1; dy <= 1; dy++)
			for (int dx = -1; dx <= 1; dx++)
				{
				const int x = cx + dx, y = cy + dy, z = cz + dz;
				if (x < 0 || y < 0 || z < 0 || x >= ncells.x || y >= ncells.y || z >= ncells.z)
					continue;

				const uint cell = (uint)((z * ncells.y + y) * ncells.x + x);
				for (uint p = cellhead[cell]; p != BIOSPRING_EMPTY_CELL; p = nextincell[p])
					{
					if (p == tid)
						continue;

					float3 axis = positions[p].xyz - here.xyz;
					const float distsq = axis.x * axis.x + axis.y * axis.y + axis.z * axis.z;
					if (distsq > cutoffsq || distsq == 0.0f)
						continue;

					if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
						continue;

					const float dist = sqrt(distsq);
					// Neighbour first, self second, as Particle::addStericForce
					// calls it. Both combination rules are symmetric, so this is
					// for the reader rather than for the arithmetic.
					const float module = biospring_steric_force_module(
					    mode, radii[p], radius, epsilons[p], epsilon, dist,
					    linearstiffness, mindistance, convert);
					sum += (axis / dist) * (stericscale * module);
					}
				}

	forces[tid].xyz += sum;
	}


// Hydrophobic attraction, gathered over its own cell list.
//
// The third of the pairwise terms, and the same shape as the two above: one
// writer per particle, each pair from both ends, the spring exclusion over the
// CSR already on the device, `cellwidth` from the grid and `cutoff` from the
// term.
//
// No separate test for "is this particle hydrophobic": the law is a product of
// the two hydrophobicities, so a particle with none contributes nothing on its
// own. The CPU's isHydrophobic() guard skips the neighbour walk entirely for
// those, which is worth it on a core walking a map and not on a device where
// every work item runs anyway.
__kernel void hydrophobic(const __global float4 * positions,
                          const __global float * hydrophobicities,
                          __global float4 * forces,
                          const __global uint * cellhead,
                          const __global uint * nextincell,
                          const float4 origin, const float cellwidth, const int4 ncells,
                          const __global Springocl * springs,
                          const __global int * springoffsets,
                          const int springsenabled,
                          const float cutoff, const float convert, const float decaylength,
                          const float hydrophobicityscale,
                          const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 here = positions[tid];
	const float h = hydrophobicities[tid];
	const float cutoffsq = cutoff * cutoff;

	const int cx = (int)floor((here.x - origin.x) / cellwidth);
	const int cy = (int)floor((here.y - origin.y) / cellwidth);
	const int cz = (int)floor((here.z - origin.z) / cellwidth);

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	for (int dz = -1; dz <= 1; dz++)
		for (int dy = -1; dy <= 1; dy++)
			for (int dx = -1; dx <= 1; dx++)
				{
				const int x = cx + dx, y = cy + dy, z = cz + dz;
				if (x < 0 || y < 0 || z < 0 || x >= ncells.x || y >= ncells.y || z >= ncells.z)
					continue;

				const uint cell = (uint)((z * ncells.y + y) * ncells.x + x);
				for (uint p = cellhead[cell]; p != BIOSPRING_EMPTY_CELL; p = nextincell[p])
					{
					if (p == tid)
						continue;

					float3 axis = positions[p].xyz - here.xyz;
					const float distsq = axis.x * axis.x + axis.y * axis.y + axis.z * axis.z;
					if (distsq > cutoffsq || distsq == 0.0f)
						continue;

					if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
						continue;

					const float dist = sqrt(distsq);
					const float module = biospring_hydrophobic_force_module(
					    hydrophobicities[p], h, dist, decaylength, convert);
					sum += (axis / dist) * (hydrophobicityscale * module);
					}
				}

	forces[tid].xyz += sum;
	}


// Tabulated torsions, gathered per particle.
// IMPALA: the implicit membrane, as a force on each particle's accessible
// surface.
//
// ONE BODY, no neighbour walk: the membrane is a profile in z and each particle
// answers to it alone. That makes this the cheapest kernel here and the one
// whose cost is exactly proportional to the particle count.
//
// WHAT IT NEEDS THAT THE OTHER TERMS DO NOT: a solvent-accessible surface per
// particle. That is why impala.enable implies --sasa. This kernel is written
// for the STATIC surface -- computed once at setup and constant afterwards, as
// examples 052 and 054 arrange by injecting it into the CDL. When FreeSASA
// runs in its dynamic mode the surface changes during the run and the host has
// to re-upload it; see SpringNetworkOpenCL::computeOpenCLSurfaces.
//
// FLAT SINGLE MEMBRANE ONLY. imp.hpp also carries a double membrane with a tube
// curvature on each, reachable only from an MDDriver client (the "dmou",
// "dmol" and "dmtc" custom data). That is not a parameterisation of this one:
// measured, the general branch returns exactly TWICE this one when its
// parameters are set to zero, because the lower membrane then coincides with
// the upper and is counted again. So the host does not run this kernel at all
// once any of the four is non-zero -- it cannot approximate a model it is not.
__kernel void impala(const __global float4 * positions,
                     const __global float * surfaces,
                     const __global float * transfers,
                     __global float4 * forces,
                     const float alip, const float alpha, const float z0,
                     const float convert, const float impscale,
                     const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float surface = surfaces[tid];
	if (surface == 0.0f) return;        // no surface, no term -- and most beads of a
	                                    // network built without --sasa are in that case

	const float fz = biospring_imp_force_z(positions[tid].z, surface, transfers[tid],
	                                       alip, alpha, z0, convert);
	forces[tid].z += fz * impscale;
	}

// The precomputed electrostatic potential map (APBS/OpenDX), as a force on each
// charge.
//
// WHY IT IS WORTH HAVING HERE. The map is CONSTANT for the whole run: it is
// read from the .dx once in SpringNetwork::_setupElectrostatic and nothing
// changes it afterwards, not even MDDriver. So it is uploaded once and never
// transferred again, and the eleven examples that carry a .dx stop needing the
// host in their step at all.
//
// Each cell carries the potential in .x and the field in .yzw: the gradient
// PRECOMPUTED on the host by PotentialGrid::compute_gradient, by central
// differences, already negated and already scaled, and zero on the boundary
// cells where a central difference has no neighbour. Nothing is differentiated
// here.
//
// NEAREST cell, no interpolation, because that is what the CPU does:
// DenseGrid::at goes through cell_coordinates, which truncates towards zero.
//
// THE MAP IS ANISOTROPIC. Eleven of the twelve in the examples have a different
// step on each axis -- 011 is 0.5895 x 0.5325 x 0.5717 A, 022 is 1.3782 x
// 0.8902 x 0.9798 -- and the counts differ too (161 x 225 x 385 for 032). Hence
// one inverse step per axis rather than a single cell width, and the row-major
// index below. The delta matrix is diagonal in every one of them, so no shear
// has to be handled; a sheared map would need more than this.
//
// The density map next door does the same lookup with its own weight and its
// own scale. Kept as two kernels rather than one parameterised on the weight:
// the two terms are independent, they are read separately in the timings, and
// nothing says the laws stay identical.
__kernel void electrostaticfield(const __global float4 * positions,
                                 const __global float * charges,
                                 __global float4 * forces,
                                 const __global float4 * cells,
                                 const float4 origin,
                                 const float4 invstep,
                                 const int4 shape,
                                 const float4 boxmin,
                                 const float4 boxmax,
                                 const float gridscale,
                                 const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 p = positions[tid];

	// Mirrors GridCoordinatesSystem::is_out_of_grid, whose upper bound is the
	// box minus 1e-6 (applied on the host, so boxmax already carries it). A
	// particle outside contributes no force and no energy, rather than throwing
	// as DenseGrid::at would.
	if (p.x < boxmin.x || p.x > boxmax.x ||
	    p.y < boxmin.y || p.y > boxmax.y ||
	    p.z < boxmin.z || p.z > boxmax.z)
		return;

	const int i = (int)((p.x - origin.x) * invstep.x);
	const int j = (int)((p.y - origin.y) * invstep.y);
	const int k = (int)((p.z - origin.z) * invstep.z);

	// The bounds test above is on the box and this one is on the indices: they
	// are not the same test, and rounding can put a particle a hair inside the
	// box and a hair past the last cell.
	if (i < 0 || j < 0 || k < 0 || i >= shape.x || j >= shape.y || k >= shape.z)
		return;

	const float4 cell = cells[((size_t)i * shape.y + j) * shape.z + k];
	const float q = charges[tid] * gridscale;

	forces[tid].xyz += cell.yzw * q;
	}

// The precomputed density map, as a force on each particle's burying factor.
//
// The same lookup as the electrostatic map above -- same .dx format, same
// PotentialGrid, same precomputed field, same nearest-cell rule, same
// anisotropy, same off-grid behaviour -- with its own weight and its own scale
// (densitygrid.scale, which is NOT the steric grid scale; see
// SpringNetwork::getDensityGridScale).
//
// On the weight: ParticleProperty initialises the burying factor to 1.0 and
// setBurying() has no caller, so today it is 1 everywhere and this kernel
// reduces to the map times the scale. It is read per particle anyway rather
// than assumed, so that wiring setBurying() up does not leave the device
// behind.
__kernel void densityfield(const __global float4 * positions,
                           const __global float * buryings,
                           __global float4 * forces,
                           const __global float4 * cells,
                           const float4 origin,
                           const float4 invstep,
                           const int4 shape,
                           const float4 boxmin,
                           const float4 boxmax,
                           const float gridscale,
                           const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 p = positions[tid];

	// Mirrors DenseGrid::is_out_of_grid, as above: a particle that has left the
	// map contributes nothing rather than throwing.
	if (p.x < boxmin.x || p.x > boxmax.x ||
	    p.y < boxmin.y || p.y > boxmax.y ||
	    p.z < boxmin.z || p.z > boxmax.z)
		return;

	const int i = (int)((p.x - origin.x) * invstep.x);
	const int j = (int)((p.y - origin.y) * invstep.y);
	const int k = (int)((p.z - origin.z) * invstep.z);

	if (i < 0 || j < 0 || k < 0 || i >= shape.x || j >= shape.y || k >= shape.z)
		return;

	const float4 cell = cells[((size_t)i * shape.y + j) * shape.z + k];
	const float w = buryings[tid] * gridscale;

	forces[tid].xyz += cell.yzw * w;
	}

//
// A torsion is a FOUR-atom term, so unlike every other kernel here it cannot
// simply own its output: one work item per torsion would have four particles to
// write and two torsions sharing an atom would race. Gathering instead -- one
// work item per particle, walking the torsions that particle takes part in and
// keeping only its own share -- costs computing each torsion four times and
// needs no atomic at all. The same trade, and the same CSR shape, as the spring
// kernel.
//
//   torsionoffsets[p] .. torsionoffsets[p+1]   this particle's entries
//   each entry packs (torsion index << 2 | slot), slot being which of the four
//   atoms this particle is.
//
// The scalars come from torsion_shared.h, which the CPU compiles too. The
// four-line vector assembly below is the one thing that cannot be shared -- a
// Vector3f is not a float3 -- so it is spelled exactly as SpringNetwork::
// computeTorsionForces spells it, signs included. Those signs are not
// guessable: the other convention passes the sum-to-zero test while being wrong
// by 5 rad/A.
__kernel void torsion(const __global float4 * positions,
                      __global float4 * forces,
                      __global float * energy,
                      const __global uint4 * torsionatoms,
                      const __global uint * torsiontable,
                      const __global uint * torsionfamily,
                      const __global float * tableenergy,   // unused here, kept for symmetry
                      const __global float * tabletorque,
                      const uint bins,
                      const __global int * torsionoffsets,
                      const __global uint * torsionentries,
                      const int familymask,
                      const float pi,
                      const float unit,
                      const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const int begin = torsionoffsets[tid];
	const int end = torsionoffsets[tid + 1];
	float3 sum = (float3)(0.0f, 0.0f, 0.0f);
	float esum = 0.0f;

	for (int k = begin; k < end; k++)
		{
		const uint packed = torsionentries[k];
		const uint ti = packed >> 2;
		const uint slot = packed & 3u;

		if (!(familymask & (1 << torsionfamily[ti])))
			continue;

		const uint4 q = torsionatoms[ti];
		const float3 p1 = positions[q.x].xyz;
		const float3 p2 = positions[q.y].xyz;
		const float3 p3 = positions[q.z].xyz;
		const float3 p4 = positions[q.w].xyz;

		const float3 b1 = p2 - p1;
		const float3 b2 = p3 - p2;
		const float3 b3 = p4 - p3;

		const float3 n1 = cross(b1, b2);
		const float3 n2 = cross(b2, b3);
		const float n1sq = dot(n1, n1);
		const float n2sq = dot(n2, n2);
		const float b2len = length(b2);
		// Three atoms in line: the dihedral is not defined.
		if (n1sq < 1e-12f || n2sq < 1e-12f || b2len < 1e-6f)
			continue;

		const float phi = atan2(b2len * dot(b1, n2), dot(n1, n2));

		const int b = biospring_torsion_bin(phi, pi, (int)bins);
		const float f = biospring_torsion_fraction(phi, pi, (int)bins);
		const uint base = torsiontable[ti] * (bins + 1);
		// A quarter, because the CSR holds each torsion from each of its four
		// atoms and this kernel runs once per atom.
		esum += 0.25f * (tableenergy[base + b] + f * (tableenergy[base + b + 1] - tableenergy[base + b]));
		const float torque = tabletorque[base + b] + f * (tabletorque[base + b + 1] - tabletorque[base + b]);
		if (torque == 0.0f)
			continue;

		const float scale = unit * torque;
		const float3 F1 = n1 * biospring_torsion_k1(scale, b2len, n1sq);
		const float3 F4 = n2 * biospring_torsion_k4(scale, b2len, n2sq);
		const float b2lensq = b2len * b2len;
		const float c1 = biospring_torsion_c(dot(b1, b2), b2lensq);
		const float c3 = biospring_torsion_c(dot(b3, b2), b2lensq);
		const float3 F2 = F1 * (-(c1 + 1.0f)) + F4 * c3;
		const float3 F3 = F1 * c1 - F4 * (c3 + 1.0f);

		// Only this particle's share. A particle appearing twice in the same
		// quadruplet would be a malformed torsion; the generator cannot emit
		// one and the reader rejects it.
		if (slot == 0u)      sum += F1;
		else if (slot == 1u) sum += F2;
		else if (slot == 2u) sum += F3;
		else                 sum += F4;
		}

	forces[tid].xyz += sum;
	energy[tid] = esum;
	}


__kernel void damping(__global float4 * forces,   const __global float4 * velocities, float damping, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) return;
	forces[tid]+=-damping*velocities[tid];
	}


__kernel void external(__global float4 * forces,   const __global float4 * externalforces, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) 
		return;
	forces[tid]+=externalforces[tid];
	}



// Same two steps, in the same order, as Particle::IntegrateEuler: the velocity
// takes the acceleration, then the position takes the UPDATED velocity.
//
// The division by the mass is what this kernel used to be missing. It read
// velocities += forces*timestep, which is only the CPU's answer when every
// particle weighs 1 Da. The zero guard matches the CPU's, which exists because
// a topology may declare a mass of 0.
__kernel void integration(__global float4 * positions, __global float4 * velocities,
                          __global float4 * forces, const __global float * masses,
                          const __global int * isdynamic, const float timestep, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) 
		return;

	// A static particle is left exactly as the CPU leaves it: not integrated,
	// and its force not reset either, because SpringNetwork's resetForce() sits
	// inside the same loop over the dynamic list. Filtering on the mass instead,
	// as this kernel used to, only agrees with that while every static particle
	// happens to be a massless ghost -- true of every example shipped today,
	// and untrue for any network built with pdb2spn --static, which freezes
	// particles without touching their masses.
	if(!isdynamic[tid])
		return;

	float mass = masses[tid];
	if(mass>0.0f)
		velocities[tid]+=(forces[tid]/mass)*timestep;
	positions[tid]+=velocities[tid]*timestep;	
	forces[tid]=(float4)0;
	}
