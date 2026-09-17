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
// cells (nsearch.hpp); this is the device's version of the same structure, and
// it is shared -- the steric and the electrostatic terms differ in their force
// law and in their cutoff, not in who is near whom.
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
		// Outside the box: linked to nothing, and no cell points at it. It is
		// then invisible to every neighbour walk, which is why the host rebuilds
		// the box rather than letting particles drift out of it.
		nextincell[p] = BIOSPRING_EMPTY_CELL;
		return;
		}

	// Atomic: several particles land in the same cell in the same instant, and
	// a plain read-modify-write loses all but one of them.
	nextincell[p] = atom_xchg(cellhead + c, p);
	}

__kernel void linearstericprobeonparticle(const __global float4 * positions,const __global float * radii,  __global float4 * forces, const uint probeid, const float proberadius,  const uint N, const float unitscale)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) return;	
	if(tid==probeid) return;
		
	float dist=0.0f, diff=0.0f, intersectdist=proberadius+radii[tid]; 
	dist=distance(positions[probeid], positions[tid]);	
	diff=dist-intersectdist;
	
	if(diff<0.0)
		{
		float4 unit=(float4)0;
		float4 dir=(float4)0;
		float4 force=(float4)0;
		dir=positions[probeid]-positions[tid];
		unit=normalize(dir);
		force=unit*diff*unitscale;
		forces[tid]+=force;
		}
	}

// Springs, gathered: work item tid owns particle tid and is the only writer of
// forces[tid], so no atomics and no serial pass.
//
// The force module comes from biospring_spring_force_module(), which is the
// SAME TEXT the CPU compiles -- see spring_shared.h, prepended to this file
// when the kernel source is embedded. The scale argument carries the force
// field's spring scale times the unit conversion; this kernel used to receive
// only the conversion and so pulled spring.scale times too weakly.
__kernel void spring(const __global float4 * positions,
                     const __global Springocl * springs,
                     const __global int * springoffsets,
                     __global float4 * forces,
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

	for(int i=begin;i<end;i++)
		{
		float3 axis = positions[springs[i].id2].xyz - here;
		float dist = length(axis);
		sum += normalize(axis) * biospring_spring_force_module(dist,
		                                                      springs[i].stiffness,
		                                                      springs[i].equilibrium,
		                                                      scale);
		}

	forces[tid].xyz += sum;
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
