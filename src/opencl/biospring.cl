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
// ONE GRID, SHARED, and cells narrower than any cutoff. What separates the
// terms is the stencil radius each kernel is handed -- how many cells out of
// its own it walks -- and not the width of a cell. See
// shared/cellgrid_shared.h, which is where that arithmetic lives and which the
// host compiles too: the two backends have to walk the same cells.
//
// The kernels still take the frame as arguments rather than reading a global,
// so the same code bins into and walks whichever grid it is handed.
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
                           const __global uchar * included,
                           const uint N)
	{
	const uint p = get_global_id(0);
	if (p >= N) return;

	// A grid may hold a SUBSET. Coulomb only ever asks about charged particles,
	// so binning the others means walking past them at every query for nothing:
	// 85% of the beads on 034, 79% on 024. A null pointer means everyone, which
	// is what the steric term wants.
	if (included != 0 && included[p] == 0)
		{
		nextincell[p] = BIOSPRING_EMPTY_CELL;
		return;
		}

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

// ============================================================================
// Neighbour list
// ============================================================================
//
// The cell walk above answers "who is near me", and the force kernels used to
// ask it again at every step. Measured with the force law compiled out of the
// kernel, that walk is 87% of what a non-bonded kernel does and the physics is
// 13%. So the answer is kept instead: these two kernels write it once, and the
// force kernels read it until something has moved far enough to invalidate it.
//
// Built at `radius` = cutoff + skin, so a pair can cross into the cutoff
// without the list having to be rebuilt. nsearch.hpp carries the same
// construction, and the same argument, for the CPU.
//
// Two passes, because the slices have to be packed: `countneighbours` says how
// long each particle's slice is, the host turns that into offsets, and
// `fillneighbours` writes the indices. A fixed capacity per particle would be
// one pass, but it would have to be sized for the worst particle and paid for
// by every one of them.
//
// TWO filters, because the two ends of a pair are not the same question.
//
// `targets` is who gets a list at all: a particle whose force nobody will ever
// read does not need its neighbours found. That is every STATIC particle -- the
// CPU only ever loops over _dynamicparticules -- and, for Coulomb, every
// uncharged one. 013.GLIC and 041 are entirely static, so today the device
// computes 25000 particles' worth of non-bonded force and throws all of it
// away; 022.RecA is 89% static and 11% charged, of which 1% is both.
//
// `candidates` is who may APPEAR in someone's list. Not the same set: a static
// charged particle exerts a perfectly good force on a dynamic one, so it has to
// stay a candidate even though it needs no list of its own.
//
// A null pointer means "everyone", which is what the steric term's candidates
// are. The CPU has had the candidate half through its per-term searchers since
// before the device existed; it had the target half through looping over the
// dynamic particles only.

#define BIOSPRING_WALK_AT_RADIUS(BODY)                                                      \
	const int cx = (int)floor((here.x - origin.x) / cellwidth);                             \
	const int cy = (int)floor((here.y - origin.y) / cellwidth);                             \
	const int cz = (int)floor((here.z - origin.z) / cellwidth);                             \
	const float radiussq = radius * radius;                                                 \
	for (int dz = -stencilradius; dz <= stencilradius; dz++)                                \
		for (int dy = -stencilradius; dy <= stencilradius; dy++)                            \
			for (int dx = -stencilradius; dx <= stencilradius; dx++)                        \
				{                                                                           \
				if (!biospring_cell_in_range(dx, dy, dz, cellwidth, radiussq))              \
					continue;                                                               \
				const int x = cx + dx, y = cy + dy, z = cz + dz;                            \
				if (x < 0 || y < 0 || z < 0 ||                                              \
				    x >= ncells.x || y >= ncells.y || z >= ncells.z)                         \
					continue;                                                               \
				const uint cell = (uint)((z * ncells.y + y) * ncells.x + x);                \
				for (uint p = cellhead[cell]; p != BIOSPRING_EMPTY_CELL; p = nextincell[p]) \
					{                                                                       \
					if (p == tid)                                                           \
						continue;                                                           \
					if (candidates != 0 && candidates[p] == 0)                              \
						continue;                                                           \
					const float3 axis = positions[p].xyz - here.xyz;                        \
					const float distsq = axis.x*axis.x + axis.y*axis.y + axis.z*axis.z;     \
					if (distsq > radiussq || distsq == 0.0f)                                \
						continue;                                                           \
					BODY                                                                    \
					}                                                                       \
				}

__kernel void countneighbours(const __global float4 * positions,
                              const __global uint * cellhead,
                              const __global uint * nextincell,
                              const float4 origin, const float cellwidth, const int4 ncells,
                              const int stencilradius,
                              const __global uchar * targets,
                              const __global uchar * candidates,
                              const float radius,
                              __global uint * counts,
                              const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	counts[tid] = 0u;
	if (targets != 0 && targets[tid] == 0) return;

	const float4 here = positions[tid];
	uint n = 0u;
	BIOSPRING_WALK_AT_RADIUS(n++;)
	counts[tid] = n;
	}

/// Whether the list just built is complete: the total the scan produced against
/// the room the host had allocated. One word, read by every kernel that would
/// otherwise walk it, so the decision never leaves the device.
__kernel void markListUsable(const __global uint * total, const uint capacity,
                             __global uint * guard)
	{
	if (get_global_id(0) != 0u)
		return;
	guard[0] = (total[0] <= capacity && capacity > 0u) ? 1u : 0u;
	}


__kernel void fillneighbours(const __global float4 * positions,
                             const __global uint * cellhead,
                             const __global uint * nextincell,
                             const float4 origin, const float cellwidth, const int4 ncells,
                             const int stencilradius,
                             const __global uchar * targets,
                             const __global uchar * candidates,
                             const float radius,
                             const __global uint * offsets,
                             __global uint * items,
                             const __global uint * listguard,
                             const uint capacity,
                             const uint N)
	{
	// Nothing to write if it would not fit; the walk falls back to the cells.
	if (listguard[0] == 0u)
		return;
	const uint tid = get_global_id(0);
	if (tid >= N) return;
	if (targets != 0 && targets[tid] == 0) return;

	const float4 here = positions[tid];
	uint at = offsets[tid];
	// Clamped as well as gated: the guard above already turned this launch off
	// when the total did not fit, and this makes a write past the end
	// impossible even if it ever did not.
	BIOSPRING_WALK_AT_RADIUS(if (at < capacity) items[at] = p; at++;)
	}

#undef BIOSPRING_WALK_AT_RADIUS


// Enumerating the candidates of particle `tid`, either out of the stored list
// or by walking the cells. One or the other per launch, chosen by
// `listoffsets != 0`: a list exists only when a skin was configured, and
// without one it would be rebuilt every step to serve a single step.
//
// The physics goes in BODY and is written once for both. The two enumerations
// have to agree about which pairs exist, and the surest way to keep them
// agreeing is for the force to be the same text.
#define BIOSPRING_FOR_EACH_CANDIDATE(BODY)                                                    \
	if (listoffsets != 0 && listguard[0] != 0u)                                               \
		{                                                                                     \
		const uint last = listoffsets[tid + 1];                                               \
		for (uint slot = listoffsets[tid]; slot < last; slot++)                               \
			{                                                                                 \
			const uint p = listitems[slot];                                                   \
			const float3 axis = positions[p].xyz - here.xyz;                                  \
			const float distsq = axis.x*axis.x + axis.y*axis.y + axis.z*axis.z;               \
			if (distsq > cutoffsq || distsq == 0.0f)                                          \
				continue;                                                                     \
			BODY                                                                              \
			}                                                                                 \
		}                                                                                     \
	else                                                                                      \
		{                                                                                     \
		const int cx = (int)floor((here.x - origin.x) / cellwidth);                           \
		const int cy = (int)floor((here.y - origin.y) / cellwidth);                           \
		const int cz = (int)floor((here.z - origin.z) / cellwidth);                           \
		for (int dz = -stencilradius; dz <= stencilradius; dz++)                              \
		for (int dy = -stencilradius; dy <= stencilradius; dy++)                              \
		for (int dx = -stencilradius; dx <= stencilradius; dx++)                              \
			{                                                                                 \
			if (!biospring_cell_in_range(dx, dy, dz, cellwidth, cutoffsq))                    \
				continue;                                                                     \
			const int x = cx + dx, y = cy + dy, z = cz + dz;                                  \
			if (x < 0 || y < 0 || z < 0 || x >= ncells.x || y >= ncells.y || z >= ncells.z)   \
				continue;                                                                     \
			const uint cell = (uint)((z * ncells.y + y) * ncells.x + x);                      \
			for (uint p = cellhead[cell]; p != BIOSPRING_EMPTY_CELL; p = nextincell[p])       \
				{                                                                             \
				if (p == tid)                                                                 \
					continue;                                                                 \
				const float3 axis = positions[p].xyz - here.xyz;                              \
				const float distsq = axis.x*axis.x + axis.y*axis.y + axis.z*axis.z;           \
				if (distsq > cutoffsq || distsq == 0.0f)                                      \
					continue;                                                                 \
				BODY                                                                          \
				}                                                                             \
			}                                                                                 \
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
                     const float scale,
                     // spring.scale ALONE. `scale` above carries it multiplied
                     // by GLOBAL_SPRING_FORCE_CONVERT, which turns kJ/mol/A
                     // into the integrator's units -- right for a force, wrong
                     // for an energy by that same 1e-4.
                     const float energyscale)
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
		// spring.scale belongs here, as ForceField::computeSpringEnergy
		// applies it on the CPU. Without it this kernel reported an energy
		// spring.scale times too small -- 65.23 kJ/mol against the CPU's 6522
		// on 034.VirusCA, whose spring.scale is 100. The FORCES were always
		// right, which is why the two backends' trajectories agree: only the
		// reported number was wrong.
		e += 0.5f * energyscale * biospring_spring_energy(dist, springs[i].stiffness, springs[i].equilibrium);
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
                            const int stencilradius,
                            const __global uint * listoffsets,
                            const __global uint * listitems,
                            const __global uint * listguard,
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

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	BIOSPRING_FOR_EACH_CANDIDATE(
		if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
			continue;
		const float dist = sqrt(distsq);
		const float module = biospring_electrostatic_force_module(
		    charges[p], q, dist, dielectric, mindistance, fourpi, convert);
		sum += (axis / dist) * (coulombscale * module);
	)
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
                     const int stencilradius,
                     const __global uint * listoffsets,
                     const __global uint * listitems,
                     const __global uint * listguard,
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

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	BIOSPRING_FOR_EACH_CANDIDATE(
		if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
			continue;
		/* Neighbour first, self second, as Particle::addStericForce calls it. Both */
		/* combination rules are symmetric, so this is for the reader. */
		const float dist = sqrt(distsq);
		const float module = biospring_steric_force_module(
		    mode, radii[p], radius, epsilons[p], epsilon, dist,
		    linearstiffness, mindistance, convert);
		sum += (axis / dist) * (stericscale * module);
	)
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
                          const int stencilradius,
                          const __global uint * listoffsets,
                          const __global uint * listitems,
                          const __global uint * listguard,
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

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);

	BIOSPRING_FOR_EACH_CANDIDATE(
		if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
			continue;
		const float dist = sqrt(distsq);
		const float module = biospring_hydrophobic_force_module(
		    hydrophobicities[p], h, dist, decaylength, convert);
		sum += (axis / dist) * (hydrophobicityscale * module);
	)
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
                      // One weight per dihedral family, in the host's family
                      // order. The torsion term has to be weighable against the
                      // rest of the model on its own: spring.scale cannot do it,
                      // because the rigid-body mesh reads the same number.
                      const __global float * familyscale,
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

		const uint fam = torsionfamily[ti];
		if (!(familymask & (1 << fam)))
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
		const float fs = familyscale[fam];
		esum += 0.25f * fs * (tableenergy[base + b] + f * (tableenergy[base + b + 1] - tableenergy[base + b]));
		const float torque = fs * (tabletorque[base + b] + f * (tabletorque[base + b + 1] - tabletorque[base + b]));
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



// ======================================================================
// HYDROGEN BONDS
//
// Unlike every other pairwise term, this one CHOOSES its pairs: a donor slot
// and an acceptor slot are consumed by a bond and stay consumed until the bond
// breaks, so the assignment persists between steps and is not a neighbour list.
// The CPU does it in four rounds of "propose, then confirm the reciprocal
// bests", and the same four rounds run here, as three small kernels, entirely
// on the device -- the host never sees the assignment, so nothing is paid to
// bring it back and send it down again.
//
// Why that is safe to parallelise: each round writes ONE best partner per
// particle, so "i's best is j and j's best is i" designates a couple without
// ambiguity, and a particle belongs to at most one such couple per round.
// Confirmations are therefore independent and need no lock -- only the i < j
// test, so that exactly one of the two acts.

// A free slot, or -1 when the particle has none left.
inline int biospring_hb_free_slot(const __global uint * offsets, const __global int * slots, uint i)
	{
	for (uint s = offsets[i]; s < offsets[i + 1]; s++)
		if (slots[s] < 0)
			return (int)s;
	return -1;
	}

inline bool biospring_hb_already_bonded(const __global uint * donoroffsets, const __global int * donorslots,
                                        const __global uint * acceptoroffsets, const __global int * acceptorslots,
                                        uint a, uint b)
	{
	for (uint s = donoroffsets[a]; s < donoroffsets[a + 1]; s++)
		if (donorslots[s] == (int)b)
			return true;
	for (uint s = acceptoroffsets[a]; s < acceptoroffsets[a + 1]; s++)
		if (acceptorslots[s] == (int)b)
			return true;
	return false;
	}

// "Where the hydrogen points", or the lone pair: away from the antecedent, or
// along the bisector when there are two (planar sp2, exact). Returns false
// when the particle has no antecedent at all, in which case the side imposes
// no direction and its weight is 1.
inline bool biospring_hb_direction(const __global float4 * positions, const __global int2 * antecedents,
                                   uint i, float3 * outdir)
	{
	int a1 = antecedents[i].x;
	if (a1 < 0)
		return false;
	float3 here = positions[i].xyz;
	float3 u1 = here - positions[a1].xyz;
	float l1 = length(u1);
	if (l1 <= 1e-6f)
		return false;
	float3 h = u1 / l1;
	int a2 = antecedents[i].y;
	if (a2 >= 0)
		{
		float3 u2 = here - positions[a2].xyz;
		float l2 = length(u2);
		if (l2 > 1e-6f)
			h += u2 / l2;
		}
	float hlen = length(h);
	if (hlen <= 1e-6f)
		return false;
	*outdir = h / hlen;
	return true;
	}

// Step 1 of a step: release any bond whose two ends have drifted past the
// cutoff. Only the donor side is walked; the acceptor side is cleared with it,
// so the two can never disagree.
__kernel void hbondBreak(const __global float4 * positions,
                         const __global uint * donoroffsets, __global int * donorslots,
                         const __global uint * acceptoroffsets, __global int * acceptorslots,
                         const float cutoff, const uint N)
	{
	uint tid = get_global_id(0);
	if (tid >= N) return;
	for (uint s = donoroffsets[tid]; s < donoroffsets[tid + 1]; s++)
		{
		int j = donorslots[s];
		if (j < 0)
			continue;
		if (distance(positions[tid].xyz, positions[j].xyz) <= cutoff)
			continue;
		donorslots[s] = -1;
		for (uint sa = acceptoroffsets[j]; sa < acceptoroffsets[j + 1]; sa++)
			if (acceptorslots[sa] == (int)tid)
				{
				acceptorslots[sa] = -1;
				break;
				}
		}
	}

// Step 2, once per round: every particle with a free slot proposes its best
// partner. Ranked by what the bond is WORTH -- the Morse well times both
// angular factors -- and not by distance, because the nearest candidate is
// very often one the angle forbids, and it would still occupy the slot.
__kernel void hbondScore(const __global float4 * positions,
                         const __global uint * donoroffsets, const __global int * donorslots,
                         const __global uint * acceptoroffsets, const __global int * acceptorslots,
                         const __global int2 * antecedents,
                         const __global int * resids, const __global int * chains,
                         const __global uint * listoffsets, const __global uint * listitems,
                         const __global uint * listguard,
                         const __global Springocl * springs, const __global int * springoffsets,
                         const int springsenabled, const int probeid,
                         const float cutoff, const float welldepth, const float equilibrium,
                         const float width, const float hbondscale,
                         __global int * nearest, __global float * strength, const uint N)
	{
	uint tid = get_global_id(0);
	if (tid >= N) return;
	nearest[tid] = -1;
	// Zero, not -infinity: a candidate whose angular factor kills it is worth
	// exactly nothing and must not take a slot from a real partner.
	strength[tid] = 0.0f;
	if (tid == (uint)probeid) return;

	bool i_donor = biospring_hb_free_slot(donoroffsets, donorslots, tid) >= 0;
	bool i_acceptor = biospring_hb_free_slot(acceptoroffsets, acceptorslots, tid) >= 0;
	if (!i_donor && !i_acceptor) return;

	float3 here = positions[tid].xyz;
	float3 dirhere;
	bool hashere = biospring_hb_direction(positions, antecedents, tid, &dirhere);

	int best = -1;
	float beststrength = 0.0f;
	uint last = listoffsets[tid + 1];
	for (uint slot = listoffsets[tid]; slot < last; slot++)
		{
		uint j = listitems[slot];
		if (j == (uint)probeid)
			continue;

		// A free donor slot needs a free acceptor slot facing it. A particle
		// that is both -- a hydroxyl -- may pair either way, but never twice
		// with the SAME partner: that would be one bond counted as two.
		bool roles_match = (i_donor && biospring_hb_free_slot(acceptoroffsets, acceptorslots, j) >= 0)
		                || (i_acceptor && biospring_hb_free_slot(donoroffsets, donorslots, j) >= 0);
		if (!roles_match)
			continue;
		if (biospring_hb_already_bonded(donoroffsets, donorslots, acceptoroffsets, acceptorslots, tid, j))
			continue;

		// A residue's own backbone N and O sit at a fixed covalent distance.
		// That is not a hydrogen bond, and without this it is invariably the
		// closest candidate and starves the real inter-residue one.
		if (resids[tid] == resids[j] && chains[tid] == chains[j])
			continue;
		if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, j))
			continue;

		float3 axis = positions[j].xyz - here;
		float dist = length(axis);
		if (dist >= cutoff || dist <= 1e-6f)
			continue;
		float3 vhat = axis / dist;

		float weight = 1.0f;
		if (hashere)
			weight *= biospring_hbond_angular_factor(dot(dirhere, vhat));
		float3 dirthere;
		if (biospring_hb_direction(positions, antecedents, j, &dirthere))
			weight *= biospring_hbond_angular_factor(dot(dirthere, -vhat));

		float s = hbondscale * biospring_hbond_energy(dist, welldepth, equilibrium, width) * weight;
		if (s < beststrength)
			{
			beststrength = s;
			best = (int)j;
			}
		}
	nearest[tid] = best;
	strength[tid] = beststrength;
	}

// Step 3, once per round: confirm the reciprocal bests. `nearest` is fixed by
// now, so the outcome does not depend on the order pairs are confirmed in --
// which is what lets this run in parallel at all. i < j so exactly one of the
// two work items acts on the couple.
__kernel void hbondConfirm(const __global int * nearest,
                           const __global uint * donoroffsets, __global int * donorslots,
                           const __global uint * acceptoroffsets, __global int * acceptorslots,
                           const uint N)
	{
	uint tid = get_global_id(0);
	if (tid >= N) return;
	int j = nearest[tid];
	if (j < 0 || (uint)j <= tid) return;
	if (nearest[j] != (int)tid) return;

	// Which way round: whoever has a free donor slot facing the other's free
	// acceptor slot. Tried donor-side-first, exactly as the CPU does.
	int sd = biospring_hb_free_slot(donoroffsets, donorslots, tid);
	int sa = biospring_hb_free_slot(acceptoroffsets, acceptorslots, (uint)j);
	if (sd >= 0 && sa >= 0)
		{
		donorslots[sd] = j;
		acceptorslots[sa] = (int)tid;
		return;
		}
	sd = biospring_hb_free_slot(donoroffsets, donorslots, (uint)j);
	sa = biospring_hb_free_slot(acceptoroffsets, acceptorslots, tid);
	if (sd >= 0 && sa >= 0)
		{
		donorslots[sd] = (int)tid;
		acceptorslots[sa] = j;
		}
	}


// The core repulsion: the SAME Morse well, on every donor/acceptor pair that
// is NOT engaged, and only below the equilibrium distance.
//
// Cutting it there introduces no discontinuity, because the Morse force is
// exactly zero at equilibrium; what it does is leave the attractive range to
// whichever pair actually won each other's slot, so a particle cannot be
// pulled by every candidate around it at once. It is a separate term from the
// bonds above and the CPU reports the two added together -- on 072 it carries
// -898.37 of the -1064.32, so leaving it out is not a detail.
//
// Evaluated from BOTH ends like every other pairwise kernel here: one writer
// per particle, no atomics, and each side takes half the pair energy so the
// total is counted once.
__kernel void hbondCoreRepulsion(const __global float4 * positions, __global float4 * forces,
                                 const __global uint * cellhead, const __global uint * nextincell,
                                 const float4 origin, const float cellwidth, const int4 ncells,
                                 const int stencilradius,
                                 const __global uint * listoffsets, const __global uint * listitems,
                         const __global uint * listguard,
                                 const __global Springocl * springs, const __global int * springoffsets,
                                 const int springsenabled,
                                 const __global uint * donoroffsets, const __global int * donorslots,
                                 const __global uint * acceptoroffsets, const __global int * acceptorslots,
                                 const float cutoff, const float welldepth, const float equilibrium,
                                 const float width, const float hbondscale, const float convert,
                                 __global float * energyper, const uint N)
	{
	const uint tid = get_global_id(0);
	if (tid >= N) return;

	const float4 here = positions[tid];
	const float cutoffsq = cutoff * cutoff;
	// Capacity, not occupancy: whether the particle CAN donate or accept.
	const bool self_donor = donoroffsets[tid + 1] > donoroffsets[tid];
	const bool self_acceptor = acceptoroffsets[tid + 1] > acceptoroffsets[tid];

	float3 sum = (float3)(0.0f, 0.0f, 0.0f);
	float e = 0.0f;

	BIOSPRING_FOR_EACH_CANDIDATE(
		if (!((self_donor && acceptoroffsets[p + 1] > acceptoroffsets[p])
		   || (self_acceptor && donoroffsets[p + 1] > donoroffsets[p])))
			continue;
		// An engaged pair is handled in full -- attraction and repulsion --
		// by hbondForce, so counting it here too would double it.
		if (biospring_hb_already_bonded(donoroffsets, donorslots, acceptoroffsets, acceptorslots, tid, p))
			continue;
		if (springsenabled && biospring_sprung_together(springs, springoffsets, tid, p))
			continue;
		const float dist = sqrt(distsq);
		if (dist >= equilibrium)
			continue;
		const float module = hbondscale * biospring_hbond_force_module(dist, welldepth, equilibrium, width, convert);
		sum += (axis / dist) * module;
		e += 0.5f * hbondscale * biospring_hbond_energy(dist, welldepth, equilibrium, width);
	)

	forces[tid].xyz += sum;
	energyper[tid] += e;
	}

// A float add that several work items may aim at the same address. OpenCL 1.2
// has no atomic float, so this is the usual compare-and-swap loop over the
// bit pattern; cl_khr_global_int32_base_atomics is what makes it legal, and
// the device advertises it.
//
// It makes the summation order between bonds unspecified, so two runs can
// differ in the last bits. That is already true of this backend (see the
// non-determinism the parity tests are written against) and the error is far
// below what the comparison against the CPU asserts.
inline void biospring_atomic_add_float(volatile __global float * address, float value)
	{
	volatile __global int * as_int = (volatile __global int *)address;
	int expected, wanted;
	do
		{
		expected = *as_int;
		wanted = as_int(as_float(expected) + value);
		}
	while (atomic_cmpxchg(as_int, expected, wanted) != expected);
	}

inline void biospring_atomic_add_float3(volatile __global float4 * forces, uint i, float3 value)
	{
	volatile __global float * base = (volatile __global float *)(forces + i);
	biospring_atomic_add_float(base + 0, value.x);
	biospring_atomic_add_float(base + 1, value.y);
	biospring_atomic_add_float(base + 2, value.z);
	}

// The force of every engaged bond. One work item per particle, walking its own
// DONOR slots, so each bond is evaluated exactly once and by the side that
// donated -- which is what the angular weight needs to know.
//
// Six atoms receive a contribution: the donor, the acceptor, and up to two
// antecedents on each side, because the angular weight makes them third
// bodies. The three sub-terms are each balanced on their OWN atoms rather than
// letting one global "donor takes the rest" absorb everything: the acceptor
// term acts on the DONOR through the same axis the donor term acts on the
// acceptor through, and folding both into one balance puts that reaction on
// the wrong atom with the wrong sign.
__kernel void hbondForce(const __global float4 * positions, volatile __global float4 * forces,
                         const __global uint * donoroffsets, const __global int * donorslots,
                         const __global int2 * antecedents,
                         const float welldepth, const float equilibrium, const float width,
                         const float hbondscale, const float convert,
                         __global float * energyper, const uint N)
	{
	uint tid = get_global_id(0);
	if (tid >= N) return;
	energyper[tid] = 0.0f;

	float3 pd = positions[tid].xyz;
	float e = 0.0f;

	for (uint s = donoroffsets[tid]; s < donoroffsets[tid + 1]; s++)
		{
		int acceptor = donorslots[s];
		if (acceptor < 0)
			continue;

		float3 pa = positions[acceptor].xyz;
		float3 v = pa - pd;
		float dist = length(v);
		if (dist < 1e-6f)
			continue;
		float3 vhat = v / dist;

		float morse  = hbondscale * biospring_hbond_energy(dist, welldepth, equilibrium, width);
		float dmorse = hbondscale * biospring_hbond_force_module(dist, welldepth, equilibrium, width, convert);

		// Each side's direction, and the lengths the gradient needs.
		int d1i = antecedents[tid].x, d2i = antecedents[tid].y;
		int a1i = antecedents[acceptor].x, a2i = antecedents[acceptor].y;
		float3 dhat = (float3)(0.0f), ahat = (float3)(0.0f);
		float3 d1 = (float3)(0.0f), d2 = (float3)(0.0f), a1 = (float3)(0.0f), a2 = (float3)(0.0f);
		float dl1 = 0.0f, dl2 = 0.0f, dhlen = 0.0f, al1 = 0.0f, al2 = 0.0f, ahlen = 0.0f;
		bool hasD = false, hasA = false;

		if (d1i >= 0)
			{
			float3 u1 = pd - positions[d1i].xyz;
			dl1 = length(u1);
			if (dl1 > 1e-6f)
				{
				d1 = u1 / dl1;
				float3 h = d1;
				if (d2i >= 0)
					{
					float3 u2 = pd - positions[d2i].xyz;
					dl2 = length(u2);
					if (dl2 > 1e-6f) { d2 = u2 / dl2; h += d2; }
					}
				dhlen = length(h);
				if (dhlen > 1e-6f) { dhat = h / dhlen; hasD = true; }
				}
			}
		if (a1i >= 0)
			{
			float3 u1 = pa - positions[a1i].xyz;
			al1 = length(u1);
			if (al1 > 1e-6f)
				{
				a1 = u1 / al1;
				float3 h = a1;
				if (a2i >= 0)
					{
					float3 u2 = pa - positions[a2i].xyz;
					al2 = length(u2);
					if (al2 > 1e-6f) { a2 = u2 / al2; h += a2; }
					}
				ahlen = length(h);
				if (ahlen > 1e-6f) { ahat = h / ahlen; hasA = true; }
				}
			}

		float cd = hasD ? dot(dhat, vhat) : 1.0f;
		float ca = hasA ? dot(ahat, -vhat) : 1.0f;
		float wd = hasD ? biospring_hbond_angular_factor(cd) : 1.0f;
		float wa = hasA ? biospring_hbond_angular_factor(ca) : 1.0f;
		float dwd = hasD ? biospring_hbond_angular_derivative(cd) : 0.0f;
		float dwa = hasA ? biospring_hbond_angular_derivative(ca) : 0.0f;
		float w = wd * wa;

		float3 f_donor = (float3)(0.0f), f_acceptor = (float3)(0.0f);
		float3 f_d1 = (float3)(0.0f), f_d2 = (float3)(0.0f);
		float3 f_a1 = (float3)(0.0f), f_a2 = (float3)(0.0f);

		// (A) radial, on the pair.
		float3 radial = vhat * (-dmorse * w);
		f_acceptor += radial;
		f_donor -= radial;

		// (B) the donor's angular factor: explicit on its antecedents, its
		// axis dependence on the acceptor, the donor balancing the three.
		if (hasD && dwd != 0.0f)
			{
			float3 t = vhat - dhat * cd;
			float g = morse * wa * dwd * convert / dhlen;
			f_d1 = (t - d1 * dot(d1, t)) * (g / dl1);
			if (d2i >= 0 && dl2 > 1e-6f)
				f_d2 = (t - d2 * dot(d2, t)) * (g / dl2);
			float3 on_acceptor = (dhat - vhat * cd) * (-morse * wa * dwd * convert / dist);
			f_acceptor += on_acceptor;
			f_donor -= (f_d1 + f_d2 + on_acceptor);
			}

		// (C) the acceptor's, the same with the roles swapped: its partner
		// direction is -vhat and its axis dependence lands on the DONOR.
		if (hasA && dwa != 0.0f)
			{
			float3 t = -vhat - ahat * ca;
			float g = morse * wd * dwa * convert / ahlen;
			f_a1 = (t - a1 * dot(a1, t)) * (g / al1);
			if (a2i >= 0 && al2 > 1e-6f)
				f_a2 = (t - a2 * dot(a2, t)) * (g / al2);
			float3 on_donor = (ahat + vhat * ca) * (-morse * wd * dwa * convert / dist);
			f_donor += on_donor;
			f_acceptor -= (f_a1 + f_a2 + on_donor);
			}

		// The donor's own contribution is the only one nobody else can be
		// writing at the same time, but it goes through the same path for the
		// sake of one expression rather than two.
		if (d1i >= 0)            biospring_atomic_add_float3(forces, (uint)d1i, f_d1);
		if (d2i >= 0)            biospring_atomic_add_float3(forces, (uint)d2i, f_d2);
		if (a1i >= 0)            biospring_atomic_add_float3(forces, (uint)a1i, f_a1);
		if (a2i >= 0)            biospring_atomic_add_float3(forces, (uint)a2i, f_a2);
		biospring_atomic_add_float3(forces, tid, f_donor);
		biospring_atomic_add_float3(forces, (uint)acceptor, f_acceptor);

		e += morse * w;
		}

	// Owned by the donor, so no atomic: a bond is counted once, by the side
	// whose slot holds it.
	energyper[tid] = e;
	}

__kernel void external(__global float4 * forces,   const __global float4 * externalforces, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) 
		return;
	forces[tid]+=externalforces[tid];
	}



// ======================================================================
// BOUNDING BOX
//
// The box the structure occupies, which is what a cell grid is measured
// against. The host used to walk every position for it, which is why the
// positions had to come down every step -- along with the frame check and the
// list rebuild criterion, all three host loops over the same array.
//
// Two passes: each work group reduces its own block, then one work item
// reduces the blocks. With 256 to a block there are 146 of them for the
// capsid, so the second pass is a hundred-odd iterations and a parallel
// version of it would cost more in launch than it saves.
//
// A non-finite coordinate is reported rather than folded in. A diverging
// structure has no box to measure, and a min/max that quietly swallowed a NaN
// would hand back a frame describing nothing.

#define BIOSPRING_BOUNDS_FLOATS 6u   // minx miny minz maxx maxy maxz

__kernel void measureBoundsBlocks(const __global float4 * positions,
                                  __global float * blockbounds, __global int * blockfinite,
                                  __local float * lmin, __local float * lmax, __local int * lfinite,
                                  const uint N)
	{
	uint gid = get_global_id(0);
	uint lid = get_local_id(0);
	uint wg = get_local_size(0);

	// A lane past the end contributes nothing: +inf to a minimum and -inf to a
	// maximum are the identities of those operations.
	float3 p = (float3)(INFINITY, INFINITY, INFINITY);
	float3 q = (float3)(-INFINITY, -INFINITY, -INFINITY);
	int finite = 1;
	if (gid < N)
		{
		float4 here = positions[gid];
		finite = (isfinite(here.x) && isfinite(here.y) && isfinite(here.z)) ? 1 : 0;
		if (finite)
			{
			p = here.xyz;
			q = here.xyz;
			}
		}
	lmin[3u*lid+0u] = p.x; lmin[3u*lid+1u] = p.y; lmin[3u*lid+2u] = p.z;
	lmax[3u*lid+0u] = q.x; lmax[3u*lid+1u] = q.y; lmax[3u*lid+2u] = q.z;
	lfinite[lid] = finite;
	barrier(CLK_LOCAL_MEM_FENCE);

	for (uint stride = wg >> 1; stride > 0u; stride >>= 1)
		{
		if (lid < stride)
			{
			for (uint d = 0u; d < 3u; d++)
				{
				lmin[3u*lid+d] = fmin(lmin[3u*lid+d], lmin[3u*(lid+stride)+d]);
				lmax[3u*lid+d] = fmax(lmax[3u*lid+d], lmax[3u*(lid+stride)+d]);
				}
			lfinite[lid] = lfinite[lid] & lfinite[lid + stride];
			}
		barrier(CLK_LOCAL_MEM_FENCE);
		}

	if (lid == 0u)
		{
		uint g = get_group_id(0);
		for (uint d = 0u; d < 3u; d++)
			{
			blockbounds[BIOSPRING_BOUNDS_FLOATS*g + d]      = lmin[d];
			blockbounds[BIOSPRING_BOUNDS_FLOATS*g + 3u + d] = lmax[d];
			}
		blockfinite[g] = lfinite[0];
		}
	}

// The blocks reduced in turn, into bounds[0..5] and finite[0].
__kernel void measureBoundsFinal(const __global float * blockbounds, const __global int * blockfinite,
                                 __global float * bounds, __global int * finite, const uint nblocks)
	{
	if (get_global_id(0) != 0u)
		return;

	// Scalars rather than private float lo[3]. A private array indexed by a
	// loop variable is legal and was silently miscompiled here: every block
	// was read, nblocks arrived correct, the per-block minima were correct --
	// and the reduction still returned the cloud's extent instead of the two
	// corners deliberately placed in the last block. Written out, it is also
	// three registers instead of an array the compiler has to place.
	float lo0 = INFINITY, lo1 = INFINITY, lo2 = INFINITY;
	float hi0 = -INFINITY, hi1 = -INFINITY, hi2 = -INFINITY;
	int allfinite = 1;
	for (uint b = 0u; b < nblocks; b++)
		{
		uint o = BIOSPRING_BOUNDS_FLOATS * b;
		lo0 = fmin(lo0, blockbounds[o + 0u]);
		lo1 = fmin(lo1, blockbounds[o + 1u]);
		lo2 = fmin(lo2, blockbounds[o + 2u]);
		hi0 = fmax(hi0, blockbounds[o + 3u]);
		hi1 = fmax(hi1, blockbounds[o + 4u]);
		hi2 = fmax(hi2, blockbounds[o + 5u]);
		allfinite = allfinite & blockfinite[b];
		}
	bounds[0] = lo0; bounds[1] = lo1; bounds[2] = lo2;
	bounds[3] = hi0; bounds[4] = hi1; bounds[5] = hi2;
	finite[0] = allfinite;
	}


// ======================================================================
// FRAME CHECK
//
// Does every particle still fall inside the grid's frame? A particle outside
// is binned nowhere and is invisible to every neighbour walk, so the frame has
// to be re-measured as soon as the first one leaves.
//
// The host pass this replaces walks every position, which is one of the three
// reasons they have to come down at every step. The test is written exactly as
// that pass writes it, negation included: `!(local >= 0)` rather than
// `local < 0`, so that a coordinate which is not a number fails it instead of
// slipping through -- every comparison with a NaN is false, and the naive form
// would call it inside.

// One buffer for every question the host used to answer by walking the
// positions, so it asks once a step instead of once per grid. Slot names are
// on the host side, BIOSPRING_FLAG_*.
// The positions the current lists were built from. Copied on the device, so
// the reference never travels either.
__kernel void snapshotPositions(const __global float4 * positions, __global float4 * reference,
                                const uint N)
	{
	uint gid = get_global_id(0);
	if (gid >= N)
		return;
	reference[gid] = positions[gid];
	}


// ---------------------------------------------------------------------------
// Every question the host used to answer by walking the positions, in ONE pass.
//
// There were four kernels here -- finite, drift, and one frame test per grid --
// and up to seven launches a step. They all read the same float4 and they all
// reduce to a yes/no, so they are one kernel over one read, and their answers
// are BITS OF ONE INT. That makes the per-work-group reduction a single AND and
// the global combine a single atomic, instead of one of each per question.
//
// The host reads the two ints ONCE, at the end of the step, in the same batch
// as the positions and behind the same finish(). Asking at the top of the next
// step instead costs a second synchronisation per step, which is the whole of
// the difference -- measured, see _enqueueDeviceFlags.
//
// Bit set = that question's answer is still yes.
#define BIOSPRING_BIT_FINITE 1
#define BIOSPRING_BIT_DRIFT  2
#define BIOSPRING_BIT_FRAME  4    // and the next grid's bit is the one above

__kernel void resetFlags(__global int * flags, const int mask, const int high)
	{
	if (get_global_id(0) != 0u)
		return;
	flags[0] = mask;   // every question asked starts at yes
	flags[1] = high;   // the bad index is a minimum: start above every index
	}


/// @param reference   Where the particles were when the lists were built.
/// @param frames      One per grid; .w carries the cell width, so the frame is
///                    one vector rather than three arguments per grid.
/// @param flags       [0] the answer bits, [1] the first non-finite index.
__kernel void checkPositions(const __global float4 * positions,
                             const __global float4 * reference,
                             const float halfskinsquared, const int probeid,
                             const __global float4 * frames,
                             const __global int4 * framencells,
                             const int questions, const uint nframes,
                             __local int * scratch, __global int * flags, const uint N)
	{
	uint gid = get_global_id(0);
	uint lid = get_local_id(0);
	uint wg = get_local_size(0);

	int ok = questions;   // a particle past the end answers yes to everything
	if (gid < N)
		{
		float4 p = positions[gid];

		if (!(isfinite(p.x) && isfinite(p.y) && isfinite(p.z)))
			{
			ok &= ~BIOSPRING_BIT_FINITE;
			// Whichever wins the race names the particle; any of them is a
			// true answer to "which one went".
			atomic_min(flags + 1, (int)gid);
			}

		// The probe moves under someone's hand, by as much as they like, and
		// it is not what the lists are about: its interactions are the probe
		// kernel's, pair by pair against everything, with no list at all.
		// Counting its drift would rebuild every list at every step of an
		// interactive session for nothing. The CPU's own searcher skips it too.
		if ((questions & BIOSPRING_BIT_DRIFT) && gid != (uint)probeid)
			{
			float3 d = p.xyz - reference[gid].xyz;
			float moved = dot(d, d);
			// Negated, as the host writes it: a NaN fails every comparison, so
			// `moved <= limit` is false for one and the list is rebuilt. The
			// naive `moved > limit` would call it unmoved.
			if (!(moved <= halfskinsquared))
				ok &= ~BIOSPRING_BIT_DRIFT;
			}

		for (uint g = 0u; g < nframes; g++)
			{
			int bit = BIOSPRING_BIT_FRAME << g;
			if (!(questions & bit))
				continue;
			float4 frame = frames[g];
			int4 ncells = framencells[g];
			float lx = (p.x - frame.x) / frame.w;
			float ly = (p.y - frame.y) / frame.w;
			float lz = (p.z - frame.z) / frame.w;
			if (!(lx >= 0.0f) || lx >= (float)ncells.x ||
			    !(ly >= 0.0f) || ly >= (float)ncells.y ||
			    !(lz >= 0.0f) || lz >= (float)ncells.z)
				ok &= ~bit;
			}
		}

	scratch[lid] = ok;
	barrier(CLK_LOCAL_MEM_FENCE);

	for (uint stride = wg >> 1; stride > 0u; stride >>= 1)
		{
		if (lid < stride)
			scratch[lid] = scratch[lid] & scratch[lid + stride];
		barrier(CLK_LOCAL_MEM_FENCE);
		}

	// One atomic per work group rather than one per particle, and only ever
	// clearing bits, so the order the groups arrive in cannot matter.
	if (lid == 0u && scratch[0] != questions)
		atomic_and(flags, scratch[0]);
	}

// ======================================================================
// EXCLUSIVE PREFIX SCAN
//
// The per-particle neighbour COUNTS become the OFFSETS where each particle's
// slice of the item array begins. The host used to do it: read the counts
// down, add them up serially, send the offsets back.
//
// On a machine with unified memory that costs almost nothing in bytes -- the
// two arrays never travel -- so this buys no speed HERE. It buys it on a
// discrete card, where 100 kB each way is a real transfer, and that is the
// machine most OpenCL deployments are.
//
// Hillis-Steele inside a block rather than Blelloch: more work in principle
// (O(n log n) against O(n)) and half the code, with no bank conflict to reason
// about, on blocks of 256.
__kernel void scanCounts(const __global uint * counts, __global uint * offsets,
                         __global uint * blocksums, __local uint * scratch, const uint N)
	{
	uint gid = get_global_id(0);
	uint lid = get_local_id(0);
	uint wg = get_local_size(0);

	uint mine = (gid < N) ? counts[gid] : 0u;
	scratch[lid] = mine;
	barrier(CLK_LOCAL_MEM_FENCE);

	for (uint stride = 1; stride < wg; stride <<= 1)
		{
		uint add = (lid >= stride) ? scratch[lid - stride] : 0u;
		barrier(CLK_LOCAL_MEM_FENCE);
		scratch[lid] += add;
		barrier(CLK_LOCAL_MEM_FENCE);
		}

	uint inclusive = scratch[lid];
	if (gid < N)
		offsets[gid] = inclusive - mine;      // exclusive is inclusive shifted
	if (lid == wg - 1u)
		blocksums[get_group_id(0)] = inclusive;
	}

// The block totals, scanned in turn. One work item walking them: there are
// N/256 -- 98 for the nucleosome, 146 for the capsid -- and a parallel scan of
// a hundred entries costs more in launch than it saves.
__kernel void scanBlockSums(__global uint * blocksums, __global uint * total, const uint nblocks)
	{
	if (get_global_id(0) != 0)
		return;
	uint running = 0;
	for (uint b = 0; b < nblocks; b++)
		{
		uint v = blocksums[b];
		blocksums[b] = running;
		running += v;
		}
	total[0] = running;
	}

// Each block's offsets shifted by everything before it, and the grand total
// parked at offsets[N] where the fill kernel looks for the last slice's end.
//
// blocksize is PASSED, not read from get_local_size(0): this kernel is
// launched over N+1 items, so its groups need not line up with the ones
// scanCounts used, and which block an offset belongs to was decided by THAT
// launch.
__kernel void addBlockSums(__global uint * offsets, const __global uint * blocksums,
                           const __global uint * total, const uint blocksize, const uint N)
	{
	uint gid = get_global_id(0);
	if (gid < N)
		offsets[gid] += blocksums[gid / blocksize];
	else if (gid == N)
		offsets[N] = total[0];
	}

// BAOAB, first half: half kick on the force standing from the previous step,
// half drift, the bath over the whole step, half drift. Mirrors
// Particle::baoabKickDriftBathDrift line for line, including the order of the
// operations -- floating point is not associative and the two backends are
// compared to the digit.
//
// The forces are zeroed here because what the closing kick needs is the force
// at the position this leaves behind, which nobody has computed yet.
//
// Static particles return before anything: they are not integrated, so they
// have no temperature to hold either.
__kernel void baoabDriftBath(__global float4 * positions, __global float4 * velocities,
                             __global float4 * forces, __global float4 * bondedforces,
                             const __global float * masses,
                             const __global int * isdynamic, const float timestep,
                             const float gamma, const float boltzmanntemperature,
                             const uint step, const uint seed, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N)
		return;
	if(!isdynamic[tid])
		return;

	// Not `half`: that is a reserved type name in OpenCL C, and naming a
	// variable after it makes the kernel fail to compile at runtime with a
	// message about a function-style cast.
	float halfstep = 0.5f * timestep;
	float mass = masses[tid];
	if(mass>0.0f)
		velocities[tid]+=((forces[tid]+bondedforces[tid])/mass)*halfstep;   // B

	positions[tid]+=velocities[tid]*halfstep;              // A

	float decay = biospring_langevin_decay(gamma, mass, timestep);
	float kick = biospring_langevin_kick(decay, mass, boltzmanntemperature);
	velocities[tid] *= decay;                              // O, over the whole step
	if(kick > 0.0f)
		{
		float4 xi = (float4)(biospring_random_normal(step, (unsigned)tid, 0, seed),
		                     biospring_random_normal(step, (unsigned)tid, 1, seed),
		                     biospring_random_normal(step, (unsigned)tid, 2, seed),
		                     0.0f);
		velocities[tid] += kick * xi;
		}

	positions[tid]+=velocities[tid]*halfstep;              // A

	forces[tid]=(float4)0;
	bondedforces[tid]=(float4)0;
	}

// BAOAB, closing half kick, with the force at the position the half above
// left. The velocity this produces is the one the step reports, which is why
// the two half kicks are not merged into one: the merged form reports a
// velocity half a step out of date.
__kernel void baoabFinalKick(__global float4 * velocities, __global float4 * forces,
                             const __global float4 * bondedforces,
                             const __global float * masses, const __global int * isdynamic,
                             const float timestep, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N)
		return;
	if(!isdynamic[tid])
		return;
	float mass = masses[tid];
	if(mass>0.0f)
		velocities[tid]+=((forces[tid]+bondedforces[tid])/mass)*(0.5f*timestep);
	}

// Same two steps, in the same order, as Particle::IntegrateEuler: the velocity
// takes the acceleration, then the position takes the UPDATED velocity.
//
// The division by the mass is what this kernel used to be missing. It read
// velocities += forces*timestep, which is only the CPU's answer when every
// particle weighs 1 Da. The zero guard matches the CPU's, which exists because
// a topology may declare a mass of 0.
// The force terms no longer all land in one array. The bonded ones accumulate
// into their own, so that they can run at the same time as the non-bonded ones
// without two kernels doing a read-modify-write on the same forces[tid]. The
// two are regrouped HERE, where an integrator already reads the force and
// already zeroes it: the regrouping costs one more load and one more store per
// particle, and no kernel of its own.
__kernel void integration(__global float4 * positions, __global float4 * velocities,
                          __global float4 * forces, __global float4 * bondedforces,
                          const __global float * masses,
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
		velocities[tid]+=((forces[tid]+bondedforces[tid])/mass)*timestep;
	positions[tid]+=velocities[tid]*timestep;	
	forces[tid]=(float4)0;
	bondedforces[tid]=(float4)0;
	}
