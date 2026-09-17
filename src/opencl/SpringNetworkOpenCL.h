#ifdef OPENCL_SUPPORT

#ifndef _SPRINGNETWORKOPENCL_H_
#define _SPRINGNETWORKOPENCL_H_

#include "SpringNetwork.h"
#ifdef OPENGL_SUPPORT
	#include "viewer/SpringNetworkViewer.h"
#endif
#define __CL_ENABLE_EXCEPTIONS

// cl.hpp is a vendored, unmodified copy of the legacy Khronos OpenCL C++
// bindings. It intentionally targets deprecated OpenCL 1.x APIs and predates
// modern C++ warning conventions; silence warnings from it here rather than
// editing third-party code to satisfy this project's own compiler flags.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma clang diagnostic ignored "-Wignored-qualifiers"
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996) // deprecated declarations
#endif

#include "cl.hpp"

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning(pop)
#endif

#include <utility>

using biospring::spn::Spring;
using biospring::spn::SpringNetwork;
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
		virtual void endRun();

		// Walks the device's cell list the way a force kernel has to, and
		// returns the particles within `cutoff` of particle `i`.
		//
		// Public because it is what makes the grid testable. A neighbour
		// structure that quietly misses pairs does not crash and does not look
		// wrong: it just makes every term built on it too weak, by an amount
		// nothing reports. The parity test compares this against the O(N^2)
		// answer, which is the only way to see it.
		//
		// Empty when no grid has been built -- no non-bonded term is enabled,
		// or the box needed more cells than this build allocates.
		std::vector<unsigned> neighborsFromCellList(unsigned i, float cutoff);

		// The cell width of the grid currently built, or 0 if there is none.
		float cellListWidth() const { return _ncellstotal == 0 ? 0.0f : _cellwidth; }



		cl::Context * getContext() ;
		static const char* oclErrorString(cl_int error);
#ifdef OPENGL_SUPPORT
		static GLuint createVBO(const void* data, int dataSize, GLenum target, GLenum usage);
#endif

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
		float * _particlemasses;
		int * _particledynamic;
		int * _particletospringindexes;



		Springocl * _springsocl;


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
		cl::Buffer _inMassBuffer;
		cl::Buffer _inDynamicBuffer;

		cl_context_properties * _contextproperties;
		cl::Context _context;
		cl::Program::Sources _source;
		cl::Program _program;
		cl::Kernel _kernelspring;
		cl::Kernel _kernelintegration;
		cl::Kernel _kerneldamping;
		cl::Kernel _kernelexternal;
		cl::Kernel _kernelblankcells;
		cl::Kernel _kernelbinparticles;

		cl::KernelFunctor _kernelfunctorspring;
		cl::KernelFunctor _kernelfunctordamping;
		cl::KernelFunctor _kernelfunctorintegration;
		cl::KernelFunctor _kernelfunctorexternal;

		// ------------------------------------------------------------------
		// Cell list
		// ------------------------------------------------------------------
		//
		// The device's answer to "which particles are near this one", shared by
		// every non-bonded term: they differ in their force law and their
		// cutoff, not in who is near whom. See biospring.cl for the structure
		// (a linked list per cell, after Bannerman's exercise 3).
		//
		// The grid is rebuilt from the box the particles currently occupy, so
		// the origin and the cell counts are recomputed rather than fixed.
		// Must match BIOSPRING_EMPTY_CELL in biospring.cl.
		static const unsigned EMPTY_CELL = static_cast<unsigned>(-1);

		cl::Buffer _cellHeadBuffer;      // one head per cell
		cl::Buffer _nextInCellBuffer;    // one successor per particle
		unsigned * _cellhead = nullptr;
		unsigned * _nextincell = nullptr;
		unsigned _ncellstotal = 0;       // 0 = no grid built yet
		cl_int4 _ncells = {{0, 0, 0, 0}};
		cl_float4 _cellorigin = {{0.0f, 0.0f, 0.0f, 0.0f}};
		float _cellwidth = 0.0f;

		// Measures the box, sizes the grid to it and fills the cell list.
		// Returns false when there is nothing to bin, or when the cutoff makes
		// no grid possible.
		bool _buildCellList(float cutoff);

		// The cutoff the grid must resolve: the largest one among the enabled
		// non-bonded terms, since one grid serves them all.
		float _cellListCutoff() const;


		cl::CommandQueue _queue;
		cl::Event _event;

		void checkErr(const char * name);
		static float distance (const float4 p1,const float4 p2) ;

		void computeOpenCLSprings();
		void computeOpenCLPositions();
		void computeOpenCLVelocities();
		void computeOpenCLMasses();
		void computeOpenCLDynamicState();
		void computeOpenCLForces();

		void computeParticleToSpringIndexes();
		void _syncParticlesFromDevice();

		// The kernels replace computeForces(), which is where the CPU fills
		// _energies as a by-product of evaluating each term. Nothing did so
		// here, so a GPU run reported no energy at all -- not zero, absent --
		// and there was no way to compare the two backends on anything but
		// coordinates. These recompute what the device actually evaluated,
		// each from the same side of the integration as the CPU.
		float _springEnergyOfCurrentState();
		void _computeEnergiesFromDeviceState();

		// Measured before the kernels run, reported after (_resetEnergies sits
		// between the two).
		float _springenergybeforestep = 0.0f;

		// Says once, at startup, which terms the .msp turns on that the device
		// does not evaluate. Silence there means a --opencl run quietly
		// computes different physics from the same .msp.
		void _warnAboutTermsTheDeviceIgnores() const;

		// Prints only what this backend measured (see SpringNetwork's).
		virtual void _displayFrameData();
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
