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

__kernel void spring(const __global float4 * positions, const  __global Springocl * springs, const __global int * startspringindexes, __global float4 * forces, const uint N, const float unitscale)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) return;
	float dist=0.0f, diff=0.0f; 
	float4 unit=(float4)0;
	float4 dir=(float4)0;
	int springstartindex= startspringindexes[tid];
	if(springstartindex<0)
		return;
	int springendindex= N;
	if((tid+1)<N)
		springendindex= startspringindexes[tid+1];
	float equilibrium=0.0, stiffness=0.0;	
	int particleid1= springs[springstartindex].id1;
	int particleid2 = 0;
	float4 positionid1=positions[particleid1];
	int i=0;
	
	for(i=springstartindex;i<springendindex;i++)
		{
		particleid2=springs[i].id2;
		dir=positions[particleid2]-positionid1;
		dist=length(dir);
		equilibrium=springs[i].equilibrium;
		stiffness=springs[i].stiffness;	
		unit=normalize(dir);
		diff=dist-equilibrium;
		forces[tid]+=unit*diff*stiffness*unitscale;
		}
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



__kernel void integration(__global float4 * positions,  __global float4 * velocities, __global float4 * forces, const float timestep, const uint N)
	{
	size_t tid = get_global_id(0);
	if(tid>=N) 
		return;
	velocities[tid]+=forces[tid]*timestep;
	positions[tid]+=velocities[tid]*timestep;	
	forces[tid]=(float4)0;
	}
