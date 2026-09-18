#ifdef OPENCL_SUPPORT

#ifndef _SPRINGNETWORKOPENCL_H_
#define _SPRINGNETWORKOPENCL_H_

#include "SpringNetwork.h"
#include "forcefield/shared/torsion_shared.h"
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

		// A cell list: the linked-list grid of biospring.cl, plus the frame it
		// is counted in. One per non-bonded term -- see the members below for
		// why the cell width has to be that term's own cutoff.
		struct CellGrid
		{
			cl::Buffer headbuffer;       // one head per cell
			cl::Buffer nextbuffer;       // one successor per particle
			unsigned * head = nullptr;
			unsigned * next = nullptr;
			unsigned ncellstotal = 0;    // 0 = never measured
			cl_int4 ncells = {{0, 0, 0, 0}};
			cl_float4 origin = {{0.0f, 0.0f, 0.0f, 0.0f}};
			float cutoff = 0.0f;         // what this term asked for
			float width = 0.0f;          // the cells actually built; >= cutoff
		};

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
		std::vector<unsigned> neighborsFromCellList(const CellGrid & grid, unsigned i,
		                                            float cutoff);

		// The grid of a given term, for the parity test to walk. Each term has
		// its own; see the CellGrid declaration for why.
		const CellGrid & stericCells() const { return _stericcells; }
		const CellGrid & electrostaticCells() const { return _electrostaticcells; }
		const CellGrid & hydrophobicCells() const { return _hydrophobiccells; }



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
		float * _particlecharges = nullptr;
		float * _particleradii = nullptr;
		float * _particleepsilons = nullptr;
		float * _particlehydrophobicities = nullptr;
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
		cl::Buffer _inChargeBuffer;
		cl::Buffer _inRadiusBuffer;
		cl::Buffer _inEpsilonBuffer;
		cl::Buffer _inHydrophobicityBuffer;

		// The precomputed electrostatic potential map, flattened row-major, one
		// float4 per cell: potential in .x, the field PotentialGrid::
		// compute_gradient already derived in .yzw.
		//
		// Uploaded once. The map is read from the .dx at setup and nothing
		// changes it afterwards, so it never crosses the bus again -- which is
		// the whole reason it is worth holding here. Held as CL_MEM_USE_HOST_PTR
		// like every other input, so on unified memory it is not copied at all.
		//
		// The frame is kept beside it because the grid is ANISOTROPIC: almost
		// every map in the examples has a different step on each axis, so what
		// the kernel needs is one inverse step per axis, not a cell width.
		// _gridboxmax already carries GridCoordinatesSystem's -1e-6.
		cl::Buffer _inElectrostaticGridBuffer;
		cl_float4 * _electrostaticgridcells = nullptr;
		size_t _electrostaticgridcellcount = 0;
		cl_float4 _gridorigin = {{0.0f, 0.0f, 0.0f, 0.0f}};
		cl_float4 _gridinvstep = {{0.0f, 0.0f, 0.0f, 0.0f}};
		cl_float4 _gridboxmin = {{0.0f, 0.0f, 0.0f, 0.0f}};
		cl_float4 _gridboxmax = {{0.0f, 0.0f, 0.0f, 0.0f}};
		cl_int4 _gridshape = {{0, 0, 0, 0}};

		// Torsions: the quadruplets and their tables, plus a CSR from each
		// particle to the torsions it takes part in. See the torsion kernel for
		// why it gathers rather than scatters.
		cl::Buffer _inTorsionAtomsBuffer;
		cl::Buffer _inTorsionTableBuffer;
		cl::Buffer _inTorsionFamilyBuffer;
		cl::Buffer _inTorsionEnergyBuffer;
		cl::Buffer _inTorsionTorqueBuffer;
		cl::Buffer _inTorsionOffsetsBuffer;
		cl::Buffer _inTorsionEntriesBuffer;
		cl_uint4 * _torsionatoms = nullptr;
		unsigned * _torsiontable = nullptr;
		unsigned * _torsionfamily = nullptr;
		float * _torsiontableenergy = nullptr;
		float * _torsiontabletorque = nullptr;
		int * _torsionoffsets = nullptr;
		unsigned * _torsionentries = nullptr;
		unsigned _nbtorsionsocl = 0;
		unsigned _nbtorsionentries = 0;
		unsigned _torsionbins = 0;
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
		cl::Kernel _kernelelectrostatic;
		cl::Kernel _kernelsteric;
		cl::Kernel _kernelhydrophobic;
		cl::Kernel _kernelelectrostaticfield;
		cl::Kernel _kerneltorsion;

		cl::KernelFunctor _kernelfunctorspring;
		cl::KernelFunctor _kernelfunctordamping;
		cl::KernelFunctor _kernelfunctorintegration;
		cl::KernelFunctor _kernelfunctorexternal;

		// ------------------------------------------------------------------
		// Cell lists
		// ------------------------------------------------------------------
		//
		// The device's answer to "which particles are near this one". See
		// biospring.cl for the structure: a linked list per cell, after
		// Bannerman's exercise 3.
		//
		// ONE PER TERM, because the cell width IS the cutoff and the three
		// cutoffs differ -- steric 8 A, electrostatic 16, hydrophobicity 15 by
		// default. A single grid at the longest of them would make the shortest
		// term walk the longest one's volume: 27 cells of 16 A is 110592 A^3
		// against 13824, so at protein density the steric term would sift some
		// 7400 candidates to find the ~144 within its own 8 A. Eight times the
		// work for the same answer. Binning is O(N) with one atomic, walking is
		// O(N x candidates), so a second bin pass is the cheap side of that
		// trade.
		//
		// Must match BIOSPRING_EMPTY_CELL in biospring.cl.
		static const unsigned EMPTY_CELL = static_cast<unsigned>(-1);

		CellGrid _stericcells;
		CellGrid _electrostaticcells;
		CellGrid _hydrophobiccells;

		// Two different things, deliberately not one function.
		//
		// The GRID -- the cells, their width, the origin they are counted from
		// -- is a frame, and it does not move every step. What moves is which
		// cell each particle is in. Measuring the frame costs a pass over every
		// position; placing the particles in it costs two kernel launches. Only
		// the second is per-step work.
		//
		// The frame is remeasured when the device reports that a particle fell
		// outside it, which is a four-byte read rather than that pass.
		bool _measureCellGrid(CellGrid & grid, float cutoff);
		void _binParticlesIntoCells(CellGrid & grid);
		bool _frameStillHolds(const CellGrid & grid) const;
		bool _buildCellList(CellGrid & grid, float cutoff);

		// Refreshes the grid of every enabled non-bonded term.
		void _updateCellLists();



		cl::CommandQueue _queue;
		cl::Event _event;

		void checkErr(const char * name);
		static float distance (const float4 p1,const float4 p2) ;

		void computeOpenCLSprings();
		void computeOpenCLPositions();
		void computeOpenCLVelocities();
		void computeOpenCLMasses();
		void computeOpenCLCharges();
		void computeOpenCLStericParameters();
		void computeOpenCLHydrophobicity();
		// Flattens the .dx potential map into _electrostaticgridcells and fills
		// the frame beside it. Called once: the map never changes.
		void computeOpenCLElectrostaticGrid();
		void computeOpenCLTorsions();

		// Which families the .msp turns on, as the bitmask the kernel takes.
		int _torsionFamilyMask() const;



		// The .msp's steric.mode as the kernel's integer. Resolved once rather
		// than compared as a string per step; see steric_shared.h.
		int _stericMode() const;
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
		void _computeEnergiesFromDeviceState();

		// Measured before the kernels run, reported after (_resetEnergies sits
		// between the two).

		// Per-particle energy, written by the spring and torsion kernels as a
		// by-product of the force they were already computing: an energy and a
		// force are the same function, one the derivative of the other, so once
		// the distance or the angle is in hand the second costs three flops or
		// one more lerp of the SAME table at the SAME index.
		//
		// Read back only on the steps that report it -- that part is not free,
		// being a transfer and a synchronisation.
		cl::Buffer _springEnergyBuffer;
		cl::Buffer _torsionEnergyBuffer;
		float * _springenergyper = nullptr;
		float * _torsionenergyper = nullptr;

		// Whether the step now running will have its energies reported.
		bool _measuringthisstep = false;
		bool _willLogAfterThisStep() const
			{
			const unsigned rate = static_cast<unsigned>(getSampleRate());
			return rate == 0 || (static_cast<unsigned>(_nbiter) + 1u) % rate == 0u;
			}

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
