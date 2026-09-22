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
		// is counted in.
		//
		// ONE grid for every non-bonded term, not one each. A cell is not a
		// cutoff -- see forcefield/shared/cellgrid_shared.h -- so a term is not
		// tied to a grid of its own size: it takes its own stencil radius out of
		// the same bins, `biospring_stencil_radius(its cutoff, this width)`.
		// That is what removes two of the three rangings per step, and it is
		// what lets the width be narrower than any cutoff.
		struct CellGrid
		{
			cl::Buffer headbuffer;       // one head per cell
			cl::Buffer nextbuffer;       // one successor per particle
			unsigned * head = nullptr;
			unsigned * next = nullptr;
			unsigned ncellstotal = 0;    // 0 = never measured
			cl_int4 ncells = {{0, 0, 0, 0}};
			cl_float4 origin = {{0.0f, 0.0f, 0.0f, 0.0f}};
			cl::Buffer includedbuffer;   // N uchars: the subset binned here, if any
			bool restricted = false;
			float requestedwidth = 0.0f; // what getCellWidth() asked for
			float width = 0.0f;          // the cells actually built; >= requested
			int maxstencil = 1;          // the widest stencil any term needs here
		};

		// A term's stored neighbours, in the same packed layout as the cells:
		// particle i owns [offsets[i], offsets[i + 1]) of `items`.
		//
		// Built at cutoff + skin so that a pair can cross into the cutoff
		// without a rebuild, and rebuilt when something has drifted half the
		// skin -- half, because both ends of a pair are stale. nsearch.hpp
		// carries the same construction for the CPU, and the same reasoning.
		//
		// `included` is what the device never had: the particles a term acts
		// on. Coulomb only concerns the charged ones, and on 034 that is 15% of
		// the beads at both ends of every pair.
	public:
		struct NeighbourList
		{
			cl::Buffer offsetsbuffer;   // N + 1 uints
			cl::Buffer itemsbuffer;     // `total` uints
			cl::Buffer countsbuffer;    // N uints, the first pass's answer
			cl::Buffer targetsbuffer;    // N uchars: who gets a list at all
			cl::Buffer candidatesbuffer; // N uchars: who may appear in one
			std::vector<unsigned> offsets;
			std::vector<unsigned> counts;
			unsigned total = 0;
			unsigned capacity = 0;      // what itemsbuffer currently holds
			float radius = 0.0f;
			bool valid = false;
			bool hastargets = false;
			// The target mask is invariant, so it is uploaded on the first
			// build and kept. It used to be rebuilt as a fresh COPY_HOST_PTR
			// buffer on every step of every term.
			bool targetsuploaded = false;
			bool hascandidates = false;
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

		// The grid, for the parity test to walk. One for every term; see the
		// CellGrid declaration for why.
		const CellGrid & cells() const { return _cells; }

		// How many times the stored neighbours were rebuilt. Zero without a
		// skin, since there is then no list; otherwise far below the step count
		// or the list is not paying for itself.
		unsigned neighbourListRebuilds() const { return _listrebuilds; }


		// The stored neighbours of particle `i`, read back off the device.
		//
		// Public for the same reason neighborsFromCellList is: a neighbour
		// structure that quietly drops pairs does not crash and does not look
		// wrong, it just makes every term built on it too weak. Comparing the
		// forces it produces is far too blunt to see that -- a dropped pair
		// near the cutoff is worth almost nothing -- so the structure itself
		// has to be held against the O(N^2) answer.
		std::vector<unsigned> neighboursFromList(const NeighbourList & list, unsigned i);

		const NeighbourList & stericList() const { return _stericlist; }
		const NeighbourList & electrostaticList() const { return _electrostaticlist; }
		const NeighbourList & hydrophobicList() const { return _hydrophobiclist; }

		// Three cells, where the CPU wants two. Measured the same way, on the
		// same two examples, against the kernels' own timers: 5.33 A is the best
		// width on both -- the optimum outright on the capsid and within 1% of it
		// on the nucleosome -- and that is the 16 A coulomb term at three cells.
		//
		// Finer than the CPU's because the device blanks and bins in parallel, so
		// what is left of visiting a cell is one scattered read. Still not free,
		// which is why it stops at three and not at the eight the geometry alone
		// would ask for: at 2.67 A the capsid's kernels take 2.405 s against
		// 1.341 s at 5.33 A.
		// No skin on the device. A tight list rebuilt every step beats any
		// skin here -- measured monotone on both 023 and 034 -- because the
		// build is parallel and costs little, while a skin only makes the
		// walk longer for every work item on every step.
		float defaultNeighborSkin() const override { return 0.0f; }

		int targetStencilRadius() const override { return 3; }



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

		// IMPALA's two per-particle inputs: the solvent-accessible surface in A2
		// and the transfer energy PER UNIT of that surface in kJ.mol-1.A-2 (the
		// .ff's sixth column -- not the seventh, which is the pairwise
		// hydrophobic term). Uploaded once, since this port covers the STATIC
		// surface; a dynamic FreeSASA would have to re-upload _particlesurfaces.
		cl::Buffer _inSurfaceBuffer;
		cl::Buffer _inTransferBuffer;
		float * _particlesurfaces = nullptr;
		float * _particletransfers = nullptr;

		// A .dx map on the device: the cells, flattened row-major into one
		// float4 each -- the scalar in .x and the field PotentialGrid::
		// compute_gradient already derived in .yzw -- plus the frame needed to
		// find a cell from a position.
		//
		// Uploaded once. Both maps are read from their .dx at setup and nothing
		// changes them afterwards, not even MDDriver, so they never cross the
		// bus again. Held as CL_MEM_USE_HOST_PTR like every other input, so on
		// unified memory they are not copied at all.
		//
		// The frame carries one INVERSE STEP PER AXIS rather than a cell width,
		// because the maps are anisotropic: eleven of the twelve in the examples
		// have a different step on each axis. boxmax already carries
		// GridCoordinatesSystem's -1e-6.
		struct MapOnDevice
		{
			cl::Buffer buffer;
			cl_float4 * cells = nullptr;
			size_t cellcount = 0;
			cl_float4 origin = {{0.0f, 0.0f, 0.0f, 0.0f}};
			cl_float4 invstep = {{0.0f, 0.0f, 0.0f, 0.0f}};
			cl_float4 boxmin = {{0.0f, 0.0f, 0.0f, 0.0f}};
			cl_float4 boxmax = {{0.0f, 0.0f, 0.0f, 0.0f}};
			cl_int4 shape = {{0, 0, 0, 0}};
		};
		MapOnDevice _electrostaticmap;
		MapOnDevice _densitymap;

		// The per-particle weight each map is multiplied by: the charge for the
		// electrostatic one, the burying factor for the density one. Same
		// kernel, different weights.
		cl::Buffer _inBuryingBuffer;

		// Where the probe kernel scatters each particle's reaction on the probe,
		// for probegather to fold into one float4. Device-only: nothing reads it
		// on the host.
		cl::Buffer _probeForceBuffer;
		float * _particleburyings = nullptr;

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
		cl::Kernel _kerneldensityfield;
		cl::Kernel _kernelprobe;
		cl::Kernel _kernelprobegather;
		cl::Kernel _kernelimpala;

		// What one unit of a profiling timestamp is worth, in nanoseconds.
		//
		// OpenCL says those timestamps ARE nanoseconds. Apple's implementation
		// returns mach ticks instead, and on Apple Silicon a tick is
		// mach_timebase_info's 125/3 = 41.6667 ns, the 24 MHz timebase. Taking
		// them for nanoseconds under-reports every kernel by that factor --
		// measured 41.8 to 42.1 against a host clock on a kernel whose duration
		// was varied over three orders of magnitude, and it is why the kernels
		// looked like 0.9% of a step when they are most of it.
		//
		// CL_DEVICE_PROFILING_TIMER_RESOLUTION does not help: this device
		// answers 1000 ns, which is neither the resolution nor the unit.
		//
		// 1.0 anywhere else, where the spec is followed.
		double _profilingtickns = 1.0;

		// The kernels of one step, with the timer each one feeds, queried after
		// the step's single synchronisation instead of one at a time.
		//
		// The queue is in-order (created with CL_QUEUE_PROFILING_ENABLE alone),
		// so a kernel cannot start before its predecessor has finished and the
		// waits bought nothing: nine cl::Event::wait() and eighteen
		// getProfilingInfo() per step, blocking the host each time, for 0.27 ms
		// of actual kernel work on a 25069-particle system. Profiling info is
		// readable once an event has completed, which the blocking read at the
		// end of the step guarantees for all of them.
		std::vector<std::pair<cl::Event, double *>> _pendingevents;
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
		// ONE, for every non-bonded term at once. There used to be one each,
		// because the cell width was the cutoff and the cutoffs differ --
		// steric 8 A, electrostatic 16, hydrophobicity 15 by default -- so a
		// shared grid meant a shared cutoff, and the shortest term would have
		// walked the longest one's volume. Once the width is its own number
		// (forcefield/shared/cellgrid_shared.h) the terms differ by their
		// stencil radius instead, which costs one integer rather than a grid,
		// a ranging and a binning pass each.
		//
		// Must match BIOSPRING_EMPTY_CELL in biospring.cl.
		static const unsigned EMPTY_CELL = static_cast<unsigned>(-1);

		CellGrid _cells;
		// The same frame binned over a subset: Coulomb walks only the charged
		// particles, the pairwise hydrophobic term only the hydrophobic ones.
		// One more cellhead array each, which is the memory this trades for not
		// walking past four beads in five.
		CellGrid _chargedcells;
		CellGrid _hydrophobiccells;


		NeighbourList _stericlist;
		NeighbourList _electrostaticlist;
		NeighbourList _hydrophobiclist;

		// Where every particle was when the lists were last built, and whether
		// building them is worth it at all -- see _updateNeighbourLists.
		std::vector<float4> _listreference;
		bool _listsarebuilt = false;
		// How many times the lists were rebuilt over the run. A list that is
		// rebuilt every step is a list that never served, which is the failure
		// mode a skin has to be checked against.
		unsigned _listrebuilds = 0;

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
		bool _measureCellGrid(CellGrid & grid, float width);
		void _binParticlesIntoCells(CellGrid & grid);
		// Gives `subset` the frame `_cells` already measured and bins the masked
		// particles into it. No remeasuring: one frame serves every subset, so
		// only the binning differs.
		void _binSubsetIntoCells(CellGrid & subset, const std::vector<unsigned char> & mask);

		// Who each term computes a force FOR (targets) and who may appear as
		// someone's neighbour (candidates). isDynamic, isCharged, isHydrophobic
		// and the probe's identity are all fixed for the run, so these are
		// built once rather than once per step -- and the same target mask
		// gates the force kernel, which is what keeps a static or uncharged
		// particle from walking a neighbourhood whose result nobody reads.
		struct TermMasks
			{
			std::vector<unsigned char> dynamic;            // steric targets
			std::vector<unsigned char> dynamiccharged;     // Coulomb targets
			std::vector<unsigned char> charged;            // Coulomb candidates
			std::vector<unsigned char> dynamichydrophobic; // hydrophobic targets
			std::vector<unsigned char> hydrophobic;        // hydrophobic candidates
			unsigned builtfor = 0;                         // particle count they were built for
			};
		TermMasks _masks;
		void _buildTermMasks();
		bool _frameStillHolds(const CellGrid & grid) const;
		bool _buildCellList(CellGrid & grid, float width);

		// The stored neighbours of every enabled pairwise term, rebuilt only
		// when the positions have moved past half the skin.
		void _updateNeighbourLists();
		void _buildNeighbourList(NeighbourList & list, float cutoff,
		                         const std::vector<unsigned char> & targets,
		                         const std::vector<unsigned char> & candidates,
		                         const CellGrid & grid);
		bool _listsNeedRebuilding();
		// The kernels that fill a list, held beside the force kernels.
		cl::Kernel _kernelcountneighbours;
		cl::Kernel _kernelfillneighbours;
		// The widest stencil any enabled term asks of cells of `width`.
		int _maxStencilRadius(float width) const;
		// The longest reach of any enabled pairwise term, in A: what the grid's
		// margin is measured against.
		float _largestPairwiseCutoff() const;

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
		// Flattens a .dx map into a MapOnDevice and fills its frame. Called once
		// per map: neither changes during a run.
		void computeOpenCLMap(MapOnDevice & map, const biospring::grid::PotentialGrid & grid,
		                      const char * what);
		void computeOpenCLBuryings();
		// Fills _particlesurfaces and _particletransfers, IMPALA's per-particle
		// inputs. Called once: this port covers the static surface.
		void computeOpenCLSurfaces();
		// Re-reads the accessible surfaces from the Particle objects and, if any
		// changed, sends them to the device. For --sasa-dynamic, where a worker
		// thread recomputes them during the run.
		void _refreshSurfacesIfChanged();

		// True when the membrane is the flat single one the kernel implements.
		// The four geometry parameters default to 0 and only an MDDriver client
		// can change them, mid-run, so this is re-read every step.
		bool _membraneIsFlat() const;
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
