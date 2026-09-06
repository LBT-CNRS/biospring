#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

typedef struct{
	unsigned id1,id2;
	float equilibrium;
	float stiffness;
	} Springocl;

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
                          const float timestep, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) 
		return;
	float mass = masses[tid];
	if(mass>0.0f)
		velocities[tid]+=(forces[tid]/mass)*timestep;
	positions[tid]+=velocities[tid]*timestep;	
	forces[tid]=(float4)0;
	}
