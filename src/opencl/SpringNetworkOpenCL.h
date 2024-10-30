#ifdef OPENCL_SUPPORT

#ifndef _SPRINGNETWORKOPENCL_H_
#define _SPRINGNETWORKOPENCL_H_

#include "SpringNetwork.h"
#ifdef OPENGL_SUPPORT
	#include "viewer/SpringNetworkViewer.h"
#endif
#define __CL_ENABLE_EXCEPTIONS
#include "cl.hpp"
#include <utility>
//#define __NO_STD_VECTOR // Use cl::vector instead of STL version


#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <vector>

using namespace std;

typedef struct
	{
	float x, y, z;
	} Vec3ocl;

typedef struct
{
	float x, y, z, w;
} float4;



typedef struct{
	Vec3ocl pos, accel, vel;
	float mass;
	float radius;
} Particleocl;

typedef struct{
	unsigned id1,id2;
	float equilibrium;
	float stiffness;
	} Springocl;


class SpringNetworkOpenCL : public SpringNetwork
	{
	public :
		SpringNetworkOpenCL();
		virtual ~SpringNetworkOpenCL();
		virtual void run();
		virtual void idleRun();
		virtual void initRun();



		cl::Context * getContext() ;
		static const char* oclErrorString(cl_int error);
		static GLuint createVBO(const void* data, int dataSize, GLenum target, GLenum usage);

		/*int _particlevelocitiesvbo;
		int _particleforcesvbo;
		int _springsvbo;
		int _particletospringindexesvbo;
		int _particleexternalforcesvbo;*/
		unsigned * _springparticlesindexes;
		void getOpenCLRessources();

	private :

		unsigned _nbparticlesocl;
		unsigned _nbspringsocl;


		float4 * _particlepositions;
		float4 * _particlevelocities;
		float4 * _particleforces;
		float4 * _particleexternalforces;
		int * _particletospringindexes;



		Springocl * _springsocl;
		Vec3ocl * _forcesocl;
		float _cutoff;


		cl_int _err;
		std::vector<cl::Memory> _allvbos;


		vector< cl::Platform > _platforms;
		vector<cl::Device> _devices;
		cl::Buffer _inoutPositionBuffer;
		cl::Buffer _inoutVelocityBuffer;
		cl::Buffer _inoutForceBuffer;

		cl::Buffer _inSpringBuffer;
		cl::Buffer _inSpringIndexesBuffer;

		cl::Buffer _inExternalForceBuffer;

		cl_context_properties * _contextproperties;
		cl::Context _context;
		cl::Program::Sources _source;
		cl::Program _program;
		cl::Kernel _kernelspring;
		cl::Kernel _kernelintegration;
		cl::Kernel _kerneldamping;
		cl::Kernel _kernelexternal;

		cl::KernelFunctor _kernelfunctorspring;
		cl::KernelFunctor _kernelfunctordamping;
		cl::KernelFunctor _kernelfunctorintegration;
		cl::KernelFunctor _kernelfunctorexternal;


		cl::CommandQueue _queue;
		cl::Event _event;

		void checkErr(const char * name);
		static float distance (const float4 p1,const float4 p2) ;

		void computeOpenCLSprings();
		void computeOpenCLPositions();
		void computeOpenCLVelocities();
		void computeOpenCLForces();

		void computeParticleToSpringIndexes();
		void wrappingOcl();
		void InitOcl();
		void createBuffer();
		static void convertSpringtoSpringocl(const Spring & spin, Springocl & spout) ;
		void computeRandomSet();
		virtual void getParticlePosition(unsigned i, float position[3]) const;
		virtual unsigned getNumberOfParticles() const;
		virtual unsigned getNumberOfSprings() const;

		virtual void setForce(unsigned i, float force[3]);




		};
#endif
#endif
