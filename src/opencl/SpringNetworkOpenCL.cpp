#ifdef OPENCL_SUPPORT



#if defined(__APPLE__)
#include <mach/mach_time.h>
#endif
#include "forcefield/energy/imp.hpp"
#include "SpringNetworkOpenCL.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/CSVSampleWriter.h"
#include "KernelSource.h"
#include <cstring>
#include <cmath>
#include "logging.h"

#include <fstream>


#include "Spring.h"
#include "Particle.h"
#include "Vector3f.h"
#include "forcefield/constants.hpp"
#include "forcefield/energy/steric.hpp"   // STERIC_LINEAR_STIFFNESS and the mode constants

using biospring::spn::Particle;
using biospring::spn::Spring;
using biospring::spn::SpringNetwork;




#ifdef OPENGL_SUPPORT
    #if defined(__APPLE__) || defined(MACOSX)
        #include <OpenGL/OpenGL.h>
    #elif defined(_WIN32)
        #include <windows.h>
    #else
        // Needed only for OpenCL/OpenGL context sharing on X11.
        #include <GL/glx.h>
    #endif
#endif



#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <stdexcept>


#define WORK_GROUP_SIZE 256


SpringNetworkOpenCL::SpringNetworkOpenCL()
    : SpringNetwork(), _springparticlesindexes(nullptr), _nbparticlesocl(0), _nbspringsocl(0),
      _particlepositions(nullptr), _particlevelocities(nullptr), _particleforces(nullptr),
      _particleexternalforces(nullptr), _particlemasses(nullptr), _particledynamic(nullptr), _particletospringindexes(nullptr), _springsocl(nullptr),
      _err(CL_SUCCESS), _contextproperties(nullptr)
    {
    getOpenCLRessources();
    }





void SpringNetworkOpenCL::getOpenCLRessources()
	{
	 cl::Platform::get(&_platforms);
	_err=_platforms.size()!=0 ? CL_SUCCESS : -1;
   	 checkErr( "cl::Platform::get");
	for(unsigned i=0;i<_platforms.size();i++)
		{
		std::cerr << "Available Platforms name are: "<<_platforms.front().getInfo<CL_PLATFORM_NAME>() <<std::endl;
		_platforms[i].getDevices( CL_DEVICE_TYPE_ALL, &_devices) ;
		_err=_devices.size()!=0 ? CL_SUCCESS : -1;
		checkErr( "cl::Device::get");
		for(unsigned i=0;i<_devices.size();i++)
			{
			std::cerr << "	Available Devices name are: "<<_devices[i].getInfo<CL_DEVICE_NAME>() <<std::endl;
			std::cerr << "	Available Devices Extensions are : "<<_devices[i].getInfo<CL_DEVICE_EXTENSIONS>() <<std::endl;
			}
		}
	}

SpringNetworkOpenCL::~SpringNetworkOpenCL()
    {
    delete[] _springparticlesindexes;
    delete[] _particlepositions;
    delete[] _particlevelocities;
    delete[] _particleforces;
    delete[] _particleexternalforces;
    delete[] _particletospringindexes;
    delete[] _particlemasses;
    delete[] _particlecharges;
    delete[] _particleradii;
    delete[] _particleepsilons;
    delete[] _particlehydrophobicities;
    delete[] _torsionatoms;
    delete[] _torsiontable;
    delete[] _torsionfamily;
    delete[] _torsiontableenergy;
    delete[] _torsiontabletorque;
    delete[] _torsionoffsets;
    delete[] _torsionentries;
    delete[] _springenergyper;
    delete[] _torsionenergyper;
    delete[] _particledynamic;
    delete[] _springsocl;
    for (CellGrid * g : {&_cells, &_chargedcells, &_hydrophobiccells, &_hydrogenbondcells})
        {
        delete[] g->head;
        delete[] g->next;
        }
    delete[] _contextproperties;
    }

cl::Context * SpringNetworkOpenCL::getContext()
	{
	return &_context;
	}



void SpringNetworkOpenCL::checkErr(const char * name)
	{
    if (_err != CL_SUCCESS)
		{
        std::cerr << "ERROR: " << name
		<< " (" << oclErrorString(_err) << ")" << std::endl;
        exit(EXIT_FAILURE);
		}
	}

void SpringNetworkOpenCL::createBuffer()
	{
	#ifdef OPENGL_SUPPORT
	if (_sharingwithgl)
		{
		_inoutPositionBuffer=cl::BufferGL(_context, CL_MEM_READ_WRITE,_viewer->getParticlesPositionVBO() , &_err);
		_allvbos.push_back(_inoutPositionBuffer);
		checkErr( "Buffer::Buffer() 1");
		glFinish();
		}
	else
		{
		// No viewer: an ordinary buffer over our own array, exactly as a build
		// without the viewer uses.
		_inoutPositionBuffer=cl::Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		                                sizeof(float4)*_nbparticlesocl, _particlepositions, &_err);
		checkErr( "Buffer::Buffer() 1");
		}
	{

		/*_particlevelocitiesvbo=createVBO(&_particlevelocities[0], _nbparticlesocl*sizeof(float4), GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
		_inoutVelocityBuffer= cl::BufferGL(_context, CL_MEM_READ_WRITE, _particlevelocitiesvbo,&_err);
		_allvbos.push_back(_inoutVelocityBuffer);
		checkErr( "Buffer::Buffer() 2");

		_particleforcesvbo=createVBO(&_particleforces[0], _nbparticlesocl*sizeof(float4), GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
		_inoutForceBuffer= cl::BufferGL(_context, CL_MEM_READ_WRITE, _particleforcesvbo,&_err);
		_allvbos.push_back(_inoutForceBuffer);
		checkErr( "Buffer::Buffer() 3");

		_springsvbo=createVBO(&_springsocl[0], sizeof(Springocl)*_nbspringsocl, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
		_inSpringBuffer=cl::BufferGL(_context,CL_MEM_READ_ONLY,_springsvbo,&_err);
		_allvbos.push_back(_inSpringBuffer);
		checkErr( "Buffer::Buffer() 4");


		_particletospringindexesvbo=createVBO(&_particletospringindexes[0], sizeof(Springocl)*_nbspringsocl, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
		_inSpringIndexesBuffer=cl::BufferGL( _context,CL_MEM_READ_ONLY,_particletospringindexesvbo, &_err);
		_allvbos.push_back(_inSpringIndexesBuffer);
		checkErr( "Buffer::Buffer() 5");


		_springsvbo=createVBO(&_particleexternalforces[0], sizeof(Springocl)*_nbspringsocl, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
		_inExternalForceBuffer=cl::BufferGL(_context,CL_MEM_READ_ONLY,_particleexternalforcesvbo,&_err);
		_allvbos.push_back(_inExternalForceBuffer);

		checkErr( "Buffer::Buffer() 6");*/
	}
	#else
		_inoutPositionBuffer=cl::Buffer(
							   _context,
							   CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
							   sizeof(float4)*_nbparticlesocl,
							   _particlepositions,
							   &_err);
		checkErr( "Buffer::Buffer() 1");

	#endif
		_inoutVelocityBuffer= cl::Buffer (
								   _context,
								   CL_MEM_READ_WRITE| CL_MEM_USE_HOST_PTR,
								   sizeof(float4)*_nbparticlesocl,
								   _particlevelocities,
								   &_err);
		checkErr( "Buffer::Buffer() 2");

		_inoutForceBuffer=cl::Buffer (
							_context,
							CL_MEM_READ_WRITE| CL_MEM_USE_HOST_PTR,
							sizeof(float4)*_nbparticlesocl,
							_particleforces,
							&_err);
		checkErr( "Buffer::Buffer() 3");


		// A network with no spring at all is legitimate -- several examples are
		// pure steric or electrostatic -- but OpenCL rejects a zero-sized
		// buffer, and the C++ wrapper turns that into an exception nobody
		// catches. The backend used to abort() on those, with
		// "cl::Error: clCreateBuffer" and nothing to say which buffer.
		if (_nbspringsocl > 0)
			{
			_inSpringBuffer=cl::Buffer(
							  _context,
							  CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
							  sizeof(Springocl)*_nbspringsocl,
							  _springsocl,
							  &_err);
			checkErr( "Buffer::Buffer() 4");
			}

		_inMassBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particlemasses,
								 &_err);
		checkErr( "Buffer::Buffer() mass");

			_inChargeBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particlecharges,
								 &_err);
		checkErr( "Buffer::Buffer() charge");

		_inRadiusBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particleradii,
								 &_err);
		checkErr( "Buffer::Buffer() radius");

		_inEpsilonBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particleepsilons,
								 &_err);
		checkErr( "Buffer::Buffer() epsilon");

		_inHydrophobicityBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particlehydrophobicities,
								 &_err);
		checkErr( "Buffer::Buffer() hydrophobicity");

		_inSurfaceBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particlesurfaces,
								 &_err);
		checkErr( "Buffer::Buffer() surface");

		_inTransferBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particletransfers,
								 &_err);
		checkErr( "Buffer::Buffer() transfer");

		_inBuryingBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float)*_nbparticlesocl,
								 _particleburyings,
								 &_err);
		checkErr( "Buffer::Buffer() burying");

		// Device-only, so no host pointer: probegather reads it back on the
		// device and nothing else ever looks at it.
		_probeForceBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_WRITE,
								 sizeof(float)*4*_nbparticlesocl,
								 NULL,
								 &_err);
		checkErr( "Buffer::Buffer() probe force");

		// Only for a map that was actually read: OpenCL rejects a zero-sized
		// buffer, and a run without the term has no map at all.
		for (MapOnDevice * map : {&_electrostaticmap, &_densitymap})
			if (map->cellcount > 0)
				{
				map->buffer=cl::Buffer(
										 _context,
										 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
										 sizeof(cl_float4)*map->cellcount,
										 map->cells,
										 &_err);
				checkErr( "Buffer::Buffer() map");
				}

		_springenergyper = new float[_nbparticlesocl];
		_torsionenergyper = new float[_nbparticlesocl];
		for (unsigned i = 0; i < _nbparticlesocl; ++i)
			_springenergyper[i] = _torsionenergyper[i] = 0.0f;
		_springEnergyBuffer = cl::Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		                                 sizeof(float) * _nbparticlesocl, _springenergyper, &_err);
		checkErr("Buffer::Buffer() spring energy");
		_torsionEnergyBuffer = cl::Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		                                  sizeof(float) * _nbparticlesocl, _torsionenergyper, &_err);
		checkErr("Buffer::Buffer() torsion energy");

		// OpenCL rejects a zero-sized buffer, and a network with no torsion at
		// all is the common case -- every example without --dihedral.
		if (_nbtorsionsocl > 0)
			{
			const unsigned samples = _torsionbins + 1;
			const unsigned ntables = _torsiontableenergy == nullptr ? 0
			    : static_cast<unsigned>(getTorsionTables().size());
			_inTorsionAtomsBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                   sizeof(cl_uint4) * _nbtorsionsocl, _torsionatoms, &_err);
			checkErr("Buffer::Buffer() torsion atoms");
			_inTorsionTableBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                   sizeof(unsigned) * _nbtorsionsocl, _torsiontable, &_err);
			checkErr("Buffer::Buffer() torsion table");
			_inTorsionFamilyBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                    sizeof(unsigned) * _nbtorsionsocl, _torsionfamily, &_err);
			checkErr("Buffer::Buffer() torsion family");
			_inTorsionEnergyBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                    sizeof(float) * ntables * samples, _torsiontableenergy, &_err);
			checkErr("Buffer::Buffer() torsion energy");
			_inTorsionTorqueBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                    sizeof(float) * ntables * samples, _torsiontabletorque, &_err);
			checkErr("Buffer::Buffer() torsion torque");
			_inTorsionOffsetsBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                     sizeof(int) * (_nbparticlesocl + 1), _torsionoffsets, &_err);
			checkErr("Buffer::Buffer() torsion offsets");
			_inTorsionEntriesBuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			                                     sizeof(unsigned) * (_nbtorsionentries == 0 ? 1 : _nbtorsionentries),
			                                     _torsionentries, &_err);
			checkErr("Buffer::Buffer() torsion entries");
			}

	_inDynamicBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(int)*_nbparticlesocl,
								 _particledynamic,
								 &_err);
		checkErr( "Buffer::Buffer() dynamic state");

		_inSpringIndexesBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(int)*(_nbparticlesocl+1),
								 _particletospringindexes,
								 &_err);



		checkErr( "Buffer::Buffer() 5");


		_inExternalForceBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(float4)*_nbparticlesocl,
								 _particleexternalforces,
								 &_err);
		checkErr( "Buffer::Buffer() 6");













	//#ifdef OPENGL_SUPPORT
	//	_platforms[0].getDevices( CL_DEVICE_TYPE_ALL, &_devices) ;

	//#else
		_devices = _context.getInfo<CL_CONTEXT_DEVICES>();
	//#endif

	_err=_devices.size() > 0 ? CL_SUCCESS : -1;
	checkErr("devices.size() > 0");

	_queue=cl::CommandQueue(_context, _devices[0], CL_QUEUE_PROFILING_ENABLE, &_err);
	checkErr("CommandQueue::CommandQueue()");

	// See _profilingtickns: APPLE'S profiling timestamps are mach ticks where
	// OpenCL says nanoseconds, a factor of 41.667 that made the kernels look
	// like 0.9% of the run instead of 86%.
	//
	// Asked of the PLATFORM and not of the operating system, because the two
	// are no longer the same question: a build that dispatches through an ICD
	// loader can be running on POCL on this very machine, and POCL's timestamps
	// are the nanoseconds the standard asks for. Correcting them would divide
	// every measurement by 41.667 in the other direction.
	#if defined(__APPLE__)
		{
		std::string platformname;
		try { platformname = _devices[0].getInfo<CL_DEVICE_PLATFORM>()
		                     ? cl::Platform(_devices[0].getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_NAME>()
		                     : std::string(); }
		catch (...) { platformname.clear(); }

		if (platformname.find("Apple") != std::string::npos)
			{
			mach_timebase_info_data_t timebase;
			if (mach_timebase_info(&timebase) == KERN_SUCCESS && timebase.denom != 0)
				_profilingtickns = static_cast<double>(timebase.numer) / timebase.denom;
			}
		}
	#endif

	// The kernel text is compiled into the binary (see KernelSource.h), so
	// there is nothing to find on disk and nothing that depends on the
	// directory biospring was launched from.
	//
	// This used to be ifstream("biospring.cl") guarded by
	//     _err = _file.is_open() ? CL_SUCCESS : -1;
	// and -1 is the value of CL_DEVICE_NOT_FOUND, so a kernel file the process
	// could not open announced itself as a machine with no GPU. The device had
	// been found and reported by name two screens earlier.
	std::string prog(biospring::opencl::KERNEL_SOURCE);

	// Still allow a file, for editing the kernel without rebuilding. If one is
	// named and cannot be read, that is what gets reported.
	if (const char * kernelpath = getenv("BIOSPRING_OPENCL_KERNEL"))
	{
		std::ifstream _file(kernelpath);
		if (!_file.is_open())
		{
			std::cerr << "ERROR: BIOSPRING_OPENCL_KERNEL names '" << kernelpath
			          << "', which cannot be read." << std::endl;
			exit(EXIT_FAILURE);
		}
		prog.assign(std::istreambuf_iterator<char>(_file), std::istreambuf_iterator<char>());
		std::cerr << "Using OpenCL kernel source from " << kernelpath << std::endl;
	}

	try
		{

		_source=cl::Program::Sources(1,std::make_pair(prog.c_str(), prog.length()+1));
		_program=cl::Program(_context, _source);
		//_err = _file.is_open() ? CL_SUCCESS : -1;
		//checkErr("Open program name");
		}
	catch (cl::Error er)
		{
		printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
		}
	try
		{
		_err = _program.build(_devices,"");
		}

	catch (cl::Error er)
		{
		printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
		}


	// The build log is what tells you why a kernel would not compile, so it is
	// printed when the build fails -- and only then. It used to go to stderr on
	// every successful run, alongside the status, the options and a bare
	// "globalsize" integer, which is how a real diagnostic gets ignored.
	if (_program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(_devices[0]) != CL_BUILD_SUCCESS)
		{
		cerr << "ERROR: could not build the OpenCL kernels." << endl;
		cerr << "Build status: " << _program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(_devices[0]) << endl;
		cerr << "Build options: " << _program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(_devices[0]) << endl;
		cerr << "Build log:" << endl
		     << _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_devices[0]) << endl;
		exit(EXIT_FAILURE);
		}

	checkErr("Program::build()");
	unsigned workgroupsize=WORK_GROUP_SIZE;
	unsigned globalsize=(_nbparticlesocl/workgroupsize)*(workgroupsize)+workgroupsize;

	_kernelspring=cl::Kernel(_program, "spring", &_err);
	_kernelfunctorspring = _kernelspring.bind(_queue, cl::NDRange(globalsize), cl::NDRange(WORK_GROUP_SIZE));
	checkErr("Kernel::Kernel()");

	_kerneldamping=cl::Kernel(_program, "damping", &_err);
	_kernelfunctordamping = _kerneldamping.bind(_queue, cl::NDRange(globalsize), cl::NDRange(WORK_GROUP_SIZE));
	checkErr("Kernel::Kernel()");


	_kernelexternal=cl::Kernel(_program, "external", &_err);
	_kernelfunctorexternal = _kernelexternal.bind(_queue, cl::NDRange(globalsize), cl::NDRange(WORK_GROUP_SIZE));
	checkErr("Kernel::Kernel()");

	_kernelintegration=cl::Kernel(_program, "integration", &_err);
	_kernelfunctorintegration = _kernelintegration.bind(_queue, cl::NDRange(globalsize), cl::NDRange(WORK_GROUP_SIZE));
	checkErr("Kernel::Kernel()");

	// The cell-list kernels are enqueued by hand rather than through a bound
	// functor: their global size follows the cell count, which changes with the
	// box, while a functor fixes its NDRange when it is bound.
	_kernelblankcells = cl::Kernel(_program, "blankCells", &_err);
	checkErr("Kernel::Kernel()");
	_kernelbinparticles = cl::Kernel(_program, "binParticles", &_err);
	checkErr("Kernel::Kernel()");
	_kernelcountneighbours = cl::Kernel(_program, "countneighbours", &_err);
	checkErr("Kernel::Kernel()");
	_kernelfillneighbours = cl::Kernel(_program, "fillneighbours", &_err);
	checkErr("Kernel::Kernel()");
	_kernelelectrostatic = cl::Kernel(_program, "electrostatic", &_err);
	checkErr("Kernel::Kernel()");
	_kernelsteric = cl::Kernel(_program, "steric", &_err);
	checkErr("Kernel::Kernel()");
	_kernelhbondbreak   = cl::Kernel(_program, "hbondBreak", &_err);
	checkErr("Kernel(hbondBreak)");
	_kernelhbondscore   = cl::Kernel(_program, "hbondScore", &_err);
	checkErr("Kernel(hbondScore)");
	_kernelhbondconfirm = cl::Kernel(_program, "hbondConfirm", &_err);
	checkErr("Kernel(hbondConfirm)");
	_kernelhbondforce   = cl::Kernel(_program, "hbondForce", &_err);
	checkErr("Kernel(hbondForce)");
	_kernelhbondrepulsion = cl::Kernel(_program, "hbondCoreRepulsion", &_err);
	checkErr("Kernel(hbondCoreRepulsion)");
	_kernelbaoabdrift = cl::Kernel(_program, "baoabDriftBath", &_err);
	checkErr("Kernel(baoabDriftBath)");
	_kernelbaoabkick  = cl::Kernel(_program, "baoabFinalKick", &_err);
	checkErr("Kernel(baoabFinalKick)");
	_kernelscancounts    = cl::Kernel(_program, "scanCounts", &_err);
	checkErr("Kernel(scanCounts)");
	_kernelscanblocksums = cl::Kernel(_program, "scanBlockSums", &_err);
	checkErr("Kernel(scanBlockSums)");
	_kerneladdblocksums  = cl::Kernel(_program, "addBlockSums", &_err);
	checkErr("Kernel(addBlockSums)");
	_kernelhydrophobic = cl::Kernel(_program, "hydrophobic", &_err);
	checkErr("Kernel::Kernel()");
	_kernelelectrostaticfield = cl::Kernel(_program, "electrostaticfield", &_err);
	checkErr("Kernel::Kernel()");
	_kerneldensityfield = cl::Kernel(_program, "densityfield", &_err);
	checkErr("Kernel::Kernel()");
	_kernelprobe = cl::Kernel(_program, "probe", &_err);
	checkErr("Kernel::Kernel()");
	_kernelprobegather = cl::Kernel(_program, "probegather", &_err);
	checkErr("Kernel::Kernel()");
	_kernelimpala = cl::Kernel(_program, "impala", &_err);
	checkErr("Kernel::Kernel()");
	_kerneltorsion = cl::Kernel(_program, "torsion", &_err);
	checkErr("Kernel::Kernel()");

	}

void SpringNetworkOpenCL::InitOcl()
	{


	#ifdef OPENGL_SUPPORT
		#if defined (__APPLE__) || defined(MACOSX)
			// Sharing with OpenGL needs an OpenGL context to share WITH, and a
			// run without the viewer has none: CGLGetCurrentContext() returns
			// null, CGLGetShareGroup(null) returns null, and the context is
			// then asked for a share group of zero -- which is
			// CL_INVALID_VALUE, followed by exit(1).
			//
			// That made a build with the viewer compiled in unable to run
			// --opencl at all from the command line, and made all eleven
			// OpenCL tests fail, on a binary whose kernels were perfectly
			// fine. Whether the viewer is COMPILED IN is a build option;
			// whether it is RUNNING is a property of this particular run, and
			// only the second one decides whether there is anything to share.
			CGLContextObj kCGLContext = CGLGetCurrentContext();
			CGLShareGroupObj kCGLShareGroup = kCGLContext ? CGLGetShareGroup(kCGLContext) : NULL;

			if (kCGLShareGroup == NULL)
				{
				biospring::logging::info("OpenCL: no OpenGL context to share with, "
				                         "using a plain device context");
				_contextproperties = new cl_context_properties[3];
				_contextproperties[0] = CL_CONTEXT_PLATFORM;
				_contextproperties[1] = (cl_context_properties)(_platforms[0])();
				_contextproperties[2] = 0;
				// The same GPU-then-anything fallback as the no-viewer build,
				// and for the same reason: an ICD loader may expose no GPU.
				try
					{
					_context = cl::Context(CL_DEVICE_TYPE_GPU, _contextproperties, NULL, NULL, &_err);
					}
				catch (const cl::Error &)
					{
					_context = cl::Context(CL_DEVICE_TYPE_ALL, _contextproperties, NULL, NULL, &_err);
					}
				checkErr("Context::Context()");
				}
			else
				{
				_sharingwithgl = true;

			_contextproperties=new cl_context_properties[3];
			_contextproperties[0] =CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE;
			_contextproperties[1] = (cl_context_properties) kCGLShareGroup;
			_contextproperties[2] =0;

			//Apple's implementation is weird, and the default values assumed by cl.hpp don't work
			//this works
			//cl_context cxGPUContext = clCreateContext(props, 0, 0, NULL, NULL, &err);
			//these dont
			//cl_context cxGPUContext = clCreateContext(props, 1,(cl_device_id*)&devices.front(), NULL, NULL, &err);
			//cl_context cxGPUContext = clCreateContextFromType(props, CL_DEVICE_TYPE_GPU, NULL, NULL, &err);

			try
				{
			    	_context = cl::Context(_contextproperties);   //had to edit line 1448 of cl.hpp to add this constructor
				}
			catch (cl::Error er)
				{
			   	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
				exit(1);
				}
				}





		#else
			#if defined(_WIN32) // Win32
			    _contextproperties=new cl_context_properties[7];
			    _contextproperties[0]=CL_GL_CONTEXT_KHR;
			    _contextproperties[1]=(cl_context_properties)wglGetCurrentContext();
			    _contextproperties[2]=CL_WGL_HDC_KHR;
			    _contextproperties[3]=(cl_context_properties)wglGetCurrentDC();
			    _contextproperties[4]=CL_CONTEXT_PLATFORM;
			    _contextproperties[5]=(cl_context_properties)(_platforms[0])();
			    _contextproperties[6]=0;


			    //cl_context cxGPUContext = clCreateContext(props, 1, &cdDevices[uiDeviceUsed], NULL, NULL, &err);
			    try
			    	{
        			_context = cl::Context(CL_DEVICE_TYPE_GPU, _contextproperties);
			    	}
			    catch (cl::Error er)
			    	{
        			printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
			    	}
			#else
			    _contextproperties=new cl_context_properties[7];
			    _contextproperties[0]=CL_GL_CONTEXT_KHR;
			    _contextproperties[1]=(cl_context_properties)glXGetCurrentContext();
			    _contextproperties[2]=CL_GLX_DISPLAY_KHR;
			    _contextproperties[3]=(cl_context_properties)glXGetCurrentDisplay();
			    _contextproperties[4]=CL_CONTEXT_PLATFORM;
			    _contextproperties[5]=(cl_context_properties)(_platforms[0])();
			    _contextproperties[6]=0;
			    //cl_context cxGPUContext = clCreateContext(props, 1, &cdDevices[uiDeviceUsed], NULL, NULL, &err);
			    try{
        			_context = cl::Context(CL_DEVICE_TYPE_GPU, _contextproperties);
			    }
			    catch (cl::Error er) {
        			printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
			    }
			#endif
		#endif
	#else
		 _contextproperties=new cl_context_properties[3];
		_contextproperties[0] = CL_CONTEXT_PLATFORM;
		_contextproperties[1] = (cl_context_properties)(_platforms[0])();
		_contextproperties[2] = 0;

		// A GPU by preference, and whatever the platform has otherwise.
		//
		// Asking for CL_DEVICE_TYPE_GPU outright is right for the
		// implementations that ship with a machine, and wrong for an ICD loader:
		// POCL exposes the CPU and nothing else, so the context creation failed
		// before it could say why. Falling back keeps the GPU first where there
		// is one -- it is an order of magnitude faster per work item, measured
		// at 1.4 ns against 18 ns here -- while letting a CPU-only platform run
		// the same kernels.
		// Caught rather than tested: this binding is built with exceptions, so
		// the constructor throws before it can report through _err.
		try
			{
			_context = cl::Context(CL_DEVICE_TYPE_GPU, _contextproperties, NULL, NULL, &_err);
			}
		catch (const cl::Error &)
			{
			biospring::logging::info("OpenCL: no GPU on platform '%s', falling back to any device",
			                         _platforms[0].getInfo<CL_PLATFORM_NAME>().c_str());
			_context = cl::Context(CL_DEVICE_TYPE_ALL, _contextproperties, NULL, NULL, &_err);
			}
		checkErr( "Context::Context()");
	#endif





	}

void SpringNetworkOpenCL::wrappingOcl()
	{
	computeOpenCLPositions();
	computeOpenCLVelocities();
	computeOpenCLForces();
	computeOpenCLMasses();
	computeOpenCLCharges();
	computeOpenCLStericParameters();
	computeOpenCLHydrophobicity();
	if (isElectrostaticFieldEnabled())
		computeOpenCLMap(_electrostaticmap, getElectrostaticGrid(), "electrostatic");
	if (isDensityGridEnabled())
		computeOpenCLMap(_densitymap, getDensityGrid(), "density");
	computeOpenCLBuryings();
	computeOpenCLSurfaces();
	computeOpenCLTorsions();
	computeOpenCLDynamicState();
	computeOpenCLSprings();
	computeParticleToSpringIndexes();
	}

double springtime=0.0;
double electrostatictime=0.0;
double sterictime=0.0;
double torsiontime=0.0;
double hydrophobictime=0.0;
double hbondtime=0.0;
double electrostaticfieldtime=0.0;
double densityfieldtime=0.0;
double probetime=0.0;
double celllisttime=0.0;
// Building the stored neighbours: countneighbours + fillneighbours, over every
// term. Separate from the cell list, which is only the blanking and the binning.
double listbuildtime=0.0;
double listbuildsterictime=0.0;
double listbuildcoulombtime=0.0;
// Which term's list is being built, so the two kernels can be charged to it.
double * listbuildinto = &listbuildtime;
double impalatime=0.0;
double dampingtime=0.0;
double integrationtime=0.0;
double externalforcetime=0.0;
double totaltime=0.0;



void SpringNetworkOpenCL::idleRun()
	{
	_pendingevents.clear();
	// The spring and torsion kernels write their energy as they go, from the
	// pre-integration positions -- which is the state the CPU reports for this
	// step. What is NOT free is reading it back, so the flag says whether this
	// step's energies will actually be looked at.
	//
	// _nbiter is still the previous step's here: SpringNetwork::idleRun()
	// increments it below, and run() tests _isTimeToLogData() after that.
	_measuringthisstep = _willLogAfterThisStep();

	// The CPU gives its interactors their turn at the TOP of a step --
	// SpringNetwork::computeStep() calls idleRun() before computeForces() --
	// while this backend calls SpringNetwork::idleRun() at the BOTTOM, after
	// the kernels. Every step is therefore fed by the previous step's sync,
	// which saw exactly the positions the CPU's sync would have seen, so the
	// two agree -- except for the very first step, which has no previous sync
	// at all and so ran with no external force where the CPU had one.
	//
	// Measured on a constant pull over 100 steps: the device applied 99 of the
	// CPU's 100 and landed at 4950 against 5050. One sync here, at the same
	// point of the same step the CPU uses, closes it exactly.
	if (_nbiter == 0)
		for (Interactor * interactor : getInteractors())
			if (interactor != nullptr)
				interactor->syncSystemStateData();

	// BAOAB's first half, when the thermostat is on: half kick on the force
	// left standing by the previous step, half drift, the bath, half drift.
	// It moves the positions BEFORE anything is evaluated on them, so the
	// grids below have to see the moved ones -- which costs the one read this
	// path has that the ordinary one does not.
	if (isThermostatEnabled())
		{
		const unsigned wg = WORK_GROUP_SIZE;
		const unsigned global = (_nbparticlesocl / wg) * wg + wg;
		unsigned a = 0;
		_kernelbaoabdrift.setArg(a++, _inoutPositionBuffer);
		_kernelbaoabdrift.setArg(a++, _inoutVelocityBuffer);
		_kernelbaoabdrift.setArg(a++, _inoutForceBuffer);
		_kernelbaoabdrift.setArg(a++, _inMassBuffer);
		_kernelbaoabdrift.setArg(a++, _inDynamicBuffer);
		_kernelbaoabdrift.setArg(a++, getTimeStep());
		_kernelbaoabdrift.setArg(a++, isViscosityEnabled() ? getViscosity() : 0.0f);
		_kernelbaoabdrift.setArg(a++, getBoltzmannTemperature());
		_kernelbaoabdrift.setArg(a++, ++_thermostatstep);
		_kernelbaoabdrift.setArg(a++, THERMOSTAT_SEED);
		_kernelbaoabdrift.setArg(a++, _nbparticlesocl);
		_err = _queue.enqueueNDRangeKernel(_kernelbaoabdrift, cl::NullRange,
		                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
		checkErr("enqueueNDRangeKernel(baoabDriftBath)");
		_pendingevents.emplace_back(_event, &integrationtime);
		_readPositionsBack();
		}

	// The neighbour structure the non-bonded terms need. Built from the box the
	// positions read back last step occupy, which is one step stale -- hence the
	// cell of margin _buildCellList adds, the same reason the CPU's grid carries
	// a skin.
	_updateCellLists();
	// And the stored neighbours the force kernels read instead of walking the
	// cells again. A no-op without a skin -- see _updateNeighbourLists.
	_updateNeighbourLists();

	#ifdef OPENGL_SUPPORT
		if (_sharingwithgl)
			{
			glFinish();
			// map OpenGL buffer object for writing from OpenCL
			_err = _queue.enqueueAcquireGLObjects(&_allvbos, NULL, &_event);
			_queue.finish();
			}
	#endif

    // Exactly what the CPU folds into a spring force: the force field's spring
    // scale and the kJ.mol-1.A-1 -> Da.A.fs-2 conversion. The kernel calls the
    // same biospring_spring_force_module() with it, so the two sides cannot
    // disagree about the magnitude. Only the conversion used to be passed,
    // which made the GPU spring.scale times too weak.
    // Nothing to gather from when there is no spring, and the buffer the kernel
    // would read does not exist.
    if (_nbspringsocl > 0)
    {
    const float springForceScale =
        getForceField()->getSpringScale() *
        static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);
    _event = _kernelfunctorspring(_inoutPositionBuffer, _inSpringBuffer,
                                 _inSpringIndexesBuffer, _inoutForceBuffer,
                                 _springEnergyBuffer,
                                 _nbparticlesocl, springForceScale,
                                 getForceField()->getSpringScale());
	_pendingevents.emplace_back(_event, &springtime);
    }


    // Tabulated torsions, over the spring term's own gate: the CPU computes
    // them inside `if (isSpringEnabled())`, right after the springs.
    if (isSpringEnabled() && _nbtorsionsocl > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const float unit = getForceField()->getSpringScale()
                         * static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT);

        unsigned a = 0;
        _kerneltorsion.setArg(a++, _inoutPositionBuffer);
        _kerneltorsion.setArg(a++, _inoutForceBuffer);
        _kerneltorsion.setArg(a++, _torsionEnergyBuffer);
        _kerneltorsion.setArg(a++, _inTorsionAtomsBuffer);
        _kerneltorsion.setArg(a++, _inTorsionTableBuffer);
        _kerneltorsion.setArg(a++, _inTorsionFamilyBuffer);
        _kerneltorsion.setArg(a++, _inTorsionEnergyBuffer);
        _kerneltorsion.setArg(a++, _inTorsionTorqueBuffer);
        _kerneltorsion.setArg(a++, _torsionbins);
        _kerneltorsion.setArg(a++, _inTorsionOffsetsBuffer);
        _kerneltorsion.setArg(a++, _inTorsionEntriesBuffer);
        _kerneltorsion.setArg(a++, _torsionFamilyMask());
        _kerneltorsion.setArg(a++, static_cast<float>(M_PI));
        _kerneltorsion.setArg(a++, unit);
        _kerneltorsion.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kerneltorsion, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(torsion)");
	_pendingevents.emplace_back(_event, &torsiontime);
        }

    // Steric, over its own cell list. Before Coulomb only because that is the
    // order SpringNetwork::computeForces uses; the two accumulate into the same
    // buffer and the sum is the same either way, to the last bit of a float.
    if (isStericEnabled() && _cells.ncellstotal > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const int springsenabled = (isSpringEnabled() && _nbspringsocl > 0) ? 1 : 0;

        unsigned a = 0;
        _kernelsteric.setArg(a++, _inoutPositionBuffer);
        _kernelsteric.setArg(a++, _inRadiusBuffer);
        _kernelsteric.setArg(a++, _inEpsilonBuffer);
        _kernelsteric.setArg(a++, _inoutForceBuffer);
        _kernelsteric.setArg(a++, _cells.headbuffer);
        _kernelsteric.setArg(a++, _cells.nextbuffer);
        _kernelsteric.setArg(a++, _cells.origin);
        _kernelsteric.setArg(a++, _cells.width);
        _kernelsteric.setArg(a++, _cells.ncells);
        // How far out of its own cell this term has to look, in cells of the
        // width its own grid built. One with cells the size of the search
        // radius, which is the default, so this is 1.
        _kernelsteric.setArg(a++, biospring_stencil_radius(getStericCutoff(), _cells.width));
        // The stored neighbours, or nothing: with no list the kernel falls back
        // to walking the cells, which is what it did before there was one.
        if (_stericlist.valid)
            {
            _kernelsteric.setArg(a++, _stericlist.offsetsbuffer);
            _kernelsteric.setArg(a++, _stericlist.itemsbuffer);
            }
        else
            {
            _kernelsteric.setArg(a++, sizeof(cl_mem), NULL);
            _kernelsteric.setArg(a++, sizeof(cl_mem), NULL);
            }
        _kernelsteric.setArg(a++, springsenabled ? _inSpringBuffer : _inMassBuffer);
        _kernelsteric.setArg(a++, _inSpringIndexesBuffer);
        _kernelsteric.setArg(a++, springsenabled);
        _kernelsteric.setArg(a++, _stericMode());
        _kernelsteric.setArg(a++, getStericCutoff());
        _kernelsteric.setArg(a++, biospring::forcefield::STERIC_LINEAR_STIFFNESS);
        _kernelsteric.setArg(a++, static_cast<float>(
            biospring::forcefield::MINIMAL_DISTANCE_VDW_CUTOFF));
        _kernelsteric.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
        _kernelsteric.setArg(a++, getForceField()->getStericScale());
        _kernelsteric.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelsteric, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(steric)");
	_pendingevents.emplace_back(_event, &sterictime);
        }

    // Coulomb, over its own cell list. Skipped when the grid could not be
    // built rather than run over a stale one: _buildCellList invalidates on
    // failure precisely so this test means what it says.
    // isElectrostaticCoulombEnabled(), not isAnyElectrostaticEnabled(): the
    // latter is also true when only the potential GRID is on, and this kernel is
    // the pairwise Coulomb one. Gated on "any", the device applied a pairwise
    // force to five examples (011, 012, 013, 021, 022) whose .msp sets
    // coulomb.enable = 0, and the CPU applied none.
    if (isElectrostaticCoulombEnabled() && _chargedcells.ncellstotal > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const int springsenabled = (isSpringEnabled() && _nbspringsocl > 0) ? 1 : 0;

        unsigned a = 0;
        _kernelelectrostatic.setArg(a++, _inoutPositionBuffer);
        _kernelelectrostatic.setArg(a++, _inChargeBuffer);
        _kernelelectrostatic.setArg(a++, _inoutForceBuffer);
        _kernelelectrostatic.setArg(a++, _chargedcells.headbuffer);
        _kernelelectrostatic.setArg(a++, _chargedcells.nextbuffer);
        _kernelelectrostatic.setArg(a++, _chargedcells.origin);
        _kernelelectrostatic.setArg(a++, _chargedcells.width);
        _kernelelectrostatic.setArg(a++, _chargedcells.ncells);
        // How far out of its own cell this term has to look, in cells of the
        // width its own grid built. One with cells the size of the search
        // radius, which is the default, so this is 1.
        _kernelelectrostatic.setArg(a++, biospring_stencil_radius(getElectrostaticCutoff(), _chargedcells.width));
        // The stored neighbours, or nothing: with no list the kernel falls back
        // to walking the cells, which is what it did before there was one.
        if (_electrostaticlist.valid)
            {
            _kernelelectrostatic.setArg(a++, _electrostaticlist.offsetsbuffer);
            _kernelelectrostatic.setArg(a++, _electrostaticlist.itemsbuffer);
            }
        else
            {
            _kernelelectrostatic.setArg(a++, sizeof(cl_mem), NULL);
            _kernelelectrostatic.setArg(a++, sizeof(cl_mem), NULL);
            }
        // A network with no spring has no spring buffer at all (OpenCL rejects
        // a zero-sized one), so hand the kernel something valid and tell it not
        // to look: the exclusion is meaningless without springs anyway.
        _kernelelectrostatic.setArg(a++, springsenabled ? _inSpringBuffer : _inMassBuffer);
        _kernelelectrostatic.setArg(a++, _inSpringIndexesBuffer);
        _kernelelectrostatic.setArg(a++, springsenabled);
        _kernelelectrostatic.setArg(a++, getElectrostaticCutoff());
        _kernelelectrostatic.setArg(a++, getForceField()->getDielectric());
        _kernelelectrostatic.setArg(a++, static_cast<float>(
            biospring::forcefield::MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF));
        _kernelelectrostatic.setArg(a++, static_cast<float>(4.0 * biospring::forcefield::PI));
        _kernelelectrostatic.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_ELECTROSTATIC_FORCE_CONVERT));
        _kernelelectrostatic.setArg(a++, getForceField()->getCoulombScale());
        _kernelelectrostatic.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelelectrostatic, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(electrostatic)");
	_pendingevents.emplace_back(_event, &electrostatictime);
        }

    // Hydrophobic attraction, over its own cell list.
    if (isHydrophobicityEnabled() && _hydrophobiccells.ncellstotal > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const int springsenabled = (isSpringEnabled() && _nbspringsocl > 0) ? 1 : 0;

        unsigned a = 0;
        _kernelhydrophobic.setArg(a++, _inoutPositionBuffer);
        _kernelhydrophobic.setArg(a++, _inHydrophobicityBuffer);
        _kernelhydrophobic.setArg(a++, _inoutForceBuffer);
        _kernelhydrophobic.setArg(a++, _hydrophobiccells.headbuffer);
        _kernelhydrophobic.setArg(a++, _hydrophobiccells.nextbuffer);
        _kernelhydrophobic.setArg(a++, _hydrophobiccells.origin);
        _kernelhydrophobic.setArg(a++, _hydrophobiccells.width);
        _kernelhydrophobic.setArg(a++, _hydrophobiccells.ncells);
        // How far out of its own cell this term has to look, in cells of the
        // width its own grid built. One with cells the size of the search
        // radius, which is the default, so this is 1.
        _kernelhydrophobic.setArg(a++, biospring_stencil_radius(getHydrophobicCutoff(), _hydrophobiccells.width));
        // The stored neighbours, or nothing: with no list the kernel falls back
        // to walking the cells, which is what it did before there was one.
        if (_hydrophobiclist.valid)
            {
            _kernelhydrophobic.setArg(a++, _hydrophobiclist.offsetsbuffer);
            _kernelhydrophobic.setArg(a++, _hydrophobiclist.itemsbuffer);
            }
        else
            {
            _kernelhydrophobic.setArg(a++, sizeof(cl_mem), NULL);
            _kernelhydrophobic.setArg(a++, sizeof(cl_mem), NULL);
            }
        _kernelhydrophobic.setArg(a++, springsenabled ? _inSpringBuffer : _inMassBuffer);
        _kernelhydrophobic.setArg(a++, _inSpringIndexesBuffer);
        _kernelhydrophobic.setArg(a++, springsenabled);
        _kernelhydrophobic.setArg(a++, getHydrophobicCutoff());
        _kernelhydrophobic.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
        _kernelhydrophobic.setArg(a++, getForceField()->getHydrophobicityDecayLength());
        _kernelhydrophobic.setArg(a++, getForceField()->getHydrophobicityScale());
        _kernelhydrophobic.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelhydrophobic, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(hydrophobic)");
	_pendingevents.emplace_back(_event, &hydrophobictime);
        }

    // The two precomputed .dx maps. No cell list and no neighbour walk: one
    // lookup per particle, in a buffer uploaded once and never touched again.
    //
    // One kernel each. They do the same lookup, but they are independent terms
    // with independent scales -- electrostaticgrid.scale for one (which is
    // getForceFieldScale(), NOT the coulomb scale: that belongs to the pairwise
    // term) and densitygrid.scale for the other -- and keeping them apart is
    // what lets the timings below tell them apart.
    if (isElectrostaticFieldEnabled() && _electrostaticmap.cellcount > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;

        unsigned a = 0;
        _kernelelectrostaticfield.setArg(a++, _inoutPositionBuffer);
        _kernelelectrostaticfield.setArg(a++, _inChargeBuffer);
        _kernelelectrostaticfield.setArg(a++, _inoutForceBuffer);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.buffer);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.origin);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.invstep);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.shape);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.boxmin);
        _kernelelectrostaticfield.setArg(a++, _electrostaticmap.boxmax);
        _kernelelectrostaticfield.setArg(a++, getForceField()->getForceFieldScale());
        _kernelelectrostaticfield.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelelectrostaticfield, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(electrostaticfield)");
        _pendingevents.emplace_back(_event, &electrostaticfieldtime);
        }

    if (isDensityGridEnabled() && _densitymap.cellcount > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;

        unsigned a = 0;
        _kerneldensityfield.setArg(a++, _inoutPositionBuffer);
        _kerneldensityfield.setArg(a++, _inBuryingBuffer);
        _kerneldensityfield.setArg(a++, _inoutForceBuffer);
        _kerneldensityfield.setArg(a++, _densitymap.buffer);
        _kerneldensityfield.setArg(a++, _densitymap.origin);
        _kerneldensityfield.setArg(a++, _densitymap.invstep);
        _kerneldensityfield.setArg(a++, _densitymap.shape);
        _kerneldensityfield.setArg(a++, _densitymap.boxmin);
        _kerneldensityfield.setArg(a++, _densitymap.boxmax);
        _kerneldensityfield.setArg(a++, getDensityGridScale());
        _kerneldensityfield.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kerneldensityfield, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(densityfield)");
        _pendingevents.emplace_back(_event, &densityfieldtime);
        }

    // IMPALA. One body, no neighbour walk, and no per-step transfer: the
    // surface was uploaded once. Skipped entirely rather than approximated when
    // a client has curved or doubled the membrane, which is a different model
    // and not a parameterisation of this one -- see the kernel.
    if (isIMPEnabled() && _membraneIsFlat() && !isFreeSASADynamic())
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;

        unsigned a = 0;
        _kernelimpala.setArg(a++, _inoutPositionBuffer);
        _kernelimpala.setArg(a++, _inSurfaceBuffer);
        _kernelimpala.setArg(a++, _inTransferBuffer);
        _kernelimpala.setArg(a++, _inoutForceBuffer);
        _kernelimpala.setArg(a++, biospring::forcefield::ALIP);
        _kernelimpala.setArg(a++, biospring::forcefield::ALPHA);
        _kernelimpala.setArg(a++, biospring::forcefield::Z0);
        _kernelimpala.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_IMP_FORCE_CONVERT));
        _kernelimpala.setArg(a++, getForceField()->getIMPScale());
        _kernelimpala.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelimpala, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(impala)");
        _pendingevents.emplace_back(_event, &impalatime);
        }

    // The interactive probe. All pairs against one partner, so no cell list;
    // the reaction on the probe is scattered per particle and folded by a
    // second, one-group kernel, because OpenCL 1.2 has no float atomic and the
    // alternative -- reading N float4 back to sum them on the host -- would
    // cost a transfer this avoids entirely.
    if (isProbeEnabled() && _probeparticule.getId() >= 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const unsigned probeid = static_cast<unsigned>(_probeparticule.getId());

        unsigned a = 0;
        _kernelprobe.setArg(a++, _inoutPositionBuffer);
        _kernelprobe.setArg(a++, _inRadiusBuffer);
        _kernelprobe.setArg(a++, _inEpsilonBuffer);
        _kernelprobe.setArg(a++, _inChargeBuffer);
        _kernelprobe.setArg(a++, _inDynamicBuffer);
        _kernelprobe.setArg(a++, _inoutForceBuffer);
        _kernelprobe.setArg(a++, _probeForceBuffer);
        _kernelprobe.setArg(a++, probeid);
        _kernelprobe.setArg(a++, isProbeStericEnabled() ? 1 : 0);
        _kernelprobe.setArg(a++, isProbeElectrostaticEnabled() ? 1 : 0);
        _kernelprobe.setArg(a++, _stericMode());
        _kernelprobe.setArg(a++, _probeparticule.getRadius());
        _kernelprobe.setArg(a++, _probeparticule.getEpsilon());
        _kernelprobe.setArg(a++, _probeparticule.getCharge());
        _kernelprobe.setArg(a++, getForceField()->getDielectric());
        _kernelprobe.setArg(a++, static_cast<float>(
            biospring::forcefield::STERIC_LINEAR_STIFFNESS));
        _kernelprobe.setArg(a++, static_cast<float>(
            biospring::forcefield::MINIMAL_DISTANCE_VDW_CUTOFF));
        _kernelprobe.setArg(a++, static_cast<float>(
            biospring::forcefield::MINIMAL_DISTANCE_ELECTROSTATIC_CUTOFF));
        _kernelprobe.setArg(a++, static_cast<float>(4.0 * biospring::forcefield::PI));
        _kernelprobe.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
        _kernelprobe.setArg(a++, static_cast<float>(
            biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
        _kernelprobe.setArg(a++, getForceField()->getStericScale());
        _kernelprobe.setArg(a++, getForceField()->getCoulombScale());
        _kernelprobe.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelprobe, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(probe)");
        _pendingevents.emplace_back(_event, &probetime);

        a = 0;
        _kernelprobegather.setArg(a++, _probeForceBuffer);
        _kernelprobegather.setArg(a++, _inoutForceBuffer);
        _kernelprobegather.setArg(a++, cl::__local(sizeof(float) * 4 * wg));
        _kernelprobegather.setArg(a++, probeid);
        _kernelprobegather.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelprobegather, cl::NullRange,
                                           cl::NDRange(wg), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(probegather)");
        _pendingevents.emplace_back(_event, &probetime);
        }

	// Hydrogen bonds. The assignment is four rounds of kernels and the force
	// is one more; nothing crosses to the host, so an interactive step pays
	// nothing for a term that chooses its own pairs.
	if (isHydrogenBondEnabled())
		{
		_uploadHydrogenBondTopology();
		_assignHydrogenBondPairsOnDevice();

		const unsigned wg = WORK_GROUP_SIZE;
		const unsigned global = (_nbparticlesocl / wg) * wg + wg;
		unsigned a = 0;
		_kernelhbondforce.setArg(a++, _inoutPositionBuffer);
		_kernelhbondforce.setArg(a++, _inoutForceBuffer);
		_kernelhbondforce.setArg(a++, _hbond.donoroffsetbuffer);
		_kernelhbondforce.setArg(a++, _hbond.donorslotbuffer);
		_kernelhbondforce.setArg(a++, _hbond.antecedentbuffer);
		_kernelhbondforce.setArg(a++, getForceField()->getHydrogenBondWellDepth());
		_kernelhbondforce.setArg(a++, getForceField()->getHydrogenBondEquilibrium());
		_kernelhbondforce.setArg(a++, getForceField()->getHydrogenBondWidth());
		_kernelhbondforce.setArg(a++, getForceField()->getHydrogenBondScale());
		_kernelhbondforce.setArg(a++, static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
		_kernelhbondforce.setArg(a++, _hbond.energybuffer);
		_kernelhbondforce.setArg(a++, _nbparticlesocl);
		_err = _queue.enqueueNDRangeKernel(_kernelhbondforce, cl::NullRange,
		                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
		checkErr("enqueueNDRangeKernel(hbondForce)");
		_pendingevents.emplace_back(_event, &hbondtime);

		// After hbondForce, which ASSIGNS the per-particle energy where this
		// one adds to it -- the CPU reports the two as one number.
		if (_hydrogenbondlist.valid && _hydrogenbondcells.ncellstotal > 0)
			{
			a = 0;
			_kernelhbondrepulsion.setArg(a++, _inoutPositionBuffer);
			_kernelhbondrepulsion.setArg(a++, _inoutForceBuffer);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondcells.headbuffer);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondcells.nextbuffer);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondcells.origin);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondcells.width);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondcells.ncells);
			_kernelhbondrepulsion.setArg(a++, biospring_stencil_radius(getHydrogenBondCutoff(),
			                                                           _hydrogenbondcells.width));
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondlist.offsetsbuffer);
			_kernelhbondrepulsion.setArg(a++, _hydrogenbondlist.itemsbuffer);
			_kernelhbondrepulsion.setArg(a++, _inSpringBuffer);
			_kernelhbondrepulsion.setArg(a++, _inSpringIndexesBuffer);
			_kernelhbondrepulsion.setArg(a++, static_cast<int>(isSpringEnabled() && _nbspringsocl > 0));
			_kernelhbondrepulsion.setArg(a++, _hbond.donoroffsetbuffer);
			_kernelhbondrepulsion.setArg(a++, _hbond.donorslotbuffer);
			_kernelhbondrepulsion.setArg(a++, _hbond.acceptoroffsetbuffer);
			_kernelhbondrepulsion.setArg(a++, _hbond.acceptorslotbuffer);
			_kernelhbondrepulsion.setArg(a++, getHydrogenBondCutoff());
			_kernelhbondrepulsion.setArg(a++, getForceField()->getHydrogenBondWellDepth());
			_kernelhbondrepulsion.setArg(a++, getForceField()->getHydrogenBondEquilibrium());
			_kernelhbondrepulsion.setArg(a++, getForceField()->getHydrogenBondWidth());
			_kernelhbondrepulsion.setArg(a++, getForceField()->getHydrogenBondScale());
			_kernelhbondrepulsion.setArg(a++, static_cast<float>(biospring::forcefield::GLOBAL_SPRING_FORCE_CONVERT));
			_kernelhbondrepulsion.setArg(a++, _hbond.energybuffer);
			_kernelhbondrepulsion.setArg(a++, _nbparticlesocl);
			_err = _queue.enqueueNDRangeKernel(_kernelhbondrepulsion, cl::NullRange,
			                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
			checkErr("enqueueNDRangeKernel(hbondCoreRepulsion)");
			_pendingevents.emplace_back(_event, &hbondtime);
			}
		}

    // Two paths, matching the CPU's for the same reason: with the thermostat
    // off the friction stays folded into the force by its own kernel, exactly
    // as before, so every existing .msp keeps its trajectory. With it on, the
    // friction moves into the integration kernel, between the kick and the
    // drift, where a bath has to act -- and where the static particles have
    // already been filtered out by the isdynamic guard.
    const float viscosity = isViscosityEnabled() ? getViscosity() : 0.0f;
    if (!isThermostatEnabled())
        {
        _event = _kernelfunctordamping(_inoutForceBuffer, _inoutVelocityBuffer,
                                      viscosity, _nbparticlesocl);
        _pendingevents.emplace_back(_event, &dampingtime);
        }


	// Nothing to add and nothing to send when no interactor pulled: this used to
	// upload the same buffer of zeros and launch a kernel to add them, at every
	// step of every run, interactive or not.
	//
	// The upload sits HERE, immediately before the kernel that reads it, rather
	// than at the end of the step where it used to be. There it ran before
	// SpringNetwork::idleRun() had given the interactors their turn, so it
	// always carried the previous sync's values.
	if (_externalforcespending)
		{
		_err = _queue.enqueueWriteBuffer(_inExternalForceBuffer, CL_FALSE, 0,
		        sizeof(float4) * _nbparticlesocl, _particleexternalforces);
		checkErr("enqueueWriteBuffer(external forces)");

		_event=_kernelfunctorexternal(_inoutForceBuffer,_inExternalForceBuffer,_nbparticlesocl);
		_pendingevents.emplace_back(_event, &externalforcetime);

		// Cleared for the next step, as the CPU clears its accumulator. The
		// write above is non-blocking, so this cannot touch the array until
		// the finish() at the end of the step has let the transfer complete.
		_externalforcespending = false;
		_externalforcesneedclearing = true;
		}

    if (isThermostatEnabled())
        {
        // BAOAB's closing half kick, with the force just evaluated at the
        // position the drift reached. The positions are already where they
        // belong; only the velocity is completed here.
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        unsigned a = 0;
        _kernelbaoabkick.setArg(a++, _inoutVelocityBuffer);
        _kernelbaoabkick.setArg(a++, _inoutForceBuffer);
        _kernelbaoabkick.setArg(a++, _inMassBuffer);
        _kernelbaoabkick.setArg(a++, _inDynamicBuffer);
        _kernelbaoabkick.setArg(a++, getTimeStep());
        _kernelbaoabkick.setArg(a++, _nbparticlesocl);
        _err = _queue.enqueueNDRangeKernel(_kernelbaoabkick, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(baoabFinalKick)");
        _pendingevents.emplace_back(_event, &integrationtime);
        }
    else
        {
        _event = _kernelfunctorintegration(_inoutPositionBuffer, _inoutVelocityBuffer,
                                          _inoutForceBuffer, _inMassBuffer, _inDynamicBuffer,
                                          getTimeStep(), _nbparticlesocl);
        _pendingevents.emplace_back(_event, &integrationtime);
        }

	// Four transfers, ONE synchronisation. These were four BLOCKING calls with a
	// finish() after two of them: six points per step where the host stopped and
	// waited for the device, to do the work of one.
	//
	// Nothing is reordered and nothing is skipped. The queue is in-order, so
	// none of these can start before the integration kernel ends whatever we
	// ask for, they run in the same sequence as before, and no host code
	// between them touches the arrays. What CL_FALSE removes is only the stall
	// between one transfer and the next. The finish() is what makes the arrays
	// valid -- and the profiling counters read just below, which need every
	// event complete.
	//
	// The cost being removed is latency, not bandwidth: the write below was
	// measured at 0.280 ms for 25 069 particles and 0.267 ms for 37 200, i.e.
	// flat in the size. On 013.GLIC, which has no pairwise term at all, this
	// block was 0.79 ms of a 0.81 ms step.
	_err = _queue.enqueueReadBuffer(_inoutVelocityBuffer, CL_FALSE, 0,
        sizeof(float4) * _nbparticlesocl, _particlevelocities);
	checkErr("enqueueReadBuffer(velocities)");

	_err = _queue.enqueueReadBuffer(_inoutPositionBuffer, CL_FALSE, 0,
        sizeof(float4) * _nbparticlesocl, _particlepositions);
	checkErr("enqueueReadBuffer(positions)");

	// The forces are no longer read back. Nothing on the host consumed them:
	// the one caller that wrote this array was setForce, which now writes the
	// external-force array instead.

	_queue.finish();

	// After the finish(), so the transfer that read it has completed.
	if (_externalforcesneedclearing)
		{
		std::memset(_particleexternalforces, 0, sizeof(float4) * _nbparticlesocl);
		_externalforcesneedclearing = false;
		}




	// Copy the device's answer back into the Particle objects before the base
	// class runs its bookkeeping, because everything downstream of here reads
	// particles, not our float4 arrays: _writeNextStep() and its PDB/XTC/CSV
	// writers, the energies, the interactors.
	//
	// Without this the GPU path computed correctly and reported the structure
	// it started from -- every frame of every trajectory identical to the
	// input, on a run that was doing real work.
	// Every kernel of this step has completed -- the reads above are blocking
	// -- so their profiling counters are all readable now.
	for (auto & pending : _pendingevents)
		{
		const cl_ulong s = pending.first.getProfilingInfo<CL_PROFILING_COMMAND_START>();
		const cl_ulong e = pending.first.getProfilingInfo<CL_PROFILING_COMMAND_END>();
		*pending.second += (e - s) * _profilingtickns * 1.0E-9;
		}
	_pendingevents.clear();

	_syncParticlesFromDevice();

	// The insertion vector is a MEASUREMENT, not a force: it reads two
	// particles and reports how deep and at what angle the structure sits in
	// the membrane. The CPU updates it at the end of SpringNetwork::
	// computeStep(), which this backend replaces, so on the device it was never
	// updated at all and every frame reported the angle and depth of the
	// initial structure. Here, right after the positions come back, is the same
	// point of the step the CPU uses.
	//
	// It is also the observable IMPALA is judged on: bead positions say whether
	// two backends agree, the insertion angle and depth say whether the membrane
	// term did what it is for.
	if (isInsertionVectorEnabled())
		_updateInsertionVector();

	// SpringNetwork::idleRun() calls _resetEnergies(), so the energies have to
	// be filled after it, not before.
	SpringNetwork::idleRun();
	_computeEnergiesFromDeviceState();

	// After the interactors, because that is where FreeSASA publishes a freshly
	// computed surface into the Particle objects. Doing it here means the next
	// step's kernel reads the surface the CPU would read on that same step.
	//
	// NOT gated on isFreeSASADynamic(). Even in its static mode FreeSASA
	// publishes AFTER setup -- a worker thread computes the areas and flips
	// _sasaValid, and until then syncParticleStateData deliberately leaves the
	// particles alone. So the surfaces uploaded at init are not the ones the CPU
	// ends up using, and a device that never re-read them ran the whole
	// trajectory on the .nc's values. Measured on 052: 7.86 A apart after 300
	// steps. The flag also only becomes true at the first interactor sync, so
	// testing it here would miss the very publication it is meant to catch.
	if (isIMPEnabled())
		_refreshSurfacesIfChanged();

	#ifdef OPENGL_SUPPORT
		if (_sharingwithgl)
			{
			//Release the VBOs so OpenGL can play with them
			_err = _queue.enqueueReleaseGLObjects(&_allvbos, NULL, &_event);
			_queue.finish();
			}
	#endif
	}

void SpringNetworkOpenCL::initRun()
	{
	SpringNetwork::initRun();
	InitOcl();
	wrappingOcl();
	createBuffer();
	}

// One summary for the run instead of one line of raw timings on stdout for
// every step, which is what idleRun() used to do. The four counters are
// accumulated there rather than overwritten, so these are run totals.
void SpringNetworkOpenCL::endRun()
	{
	totaltime=springtime+torsiontime+sterictime+electrostatictime+electrostaticfieldtime+densityfieldtime+probetime+celllisttime+listbuildtime+listbuildsterictime+listbuildcoulombtime+impalatime+hydrophobictime+dampingtime+integrationtime+externalforcetime;
	std::cout<<"OpenCL kernel time: "<<totaltime<<" s ( spring: "<<springtime
	         <<", torsion: "<<torsiontime
	         <<", steric: "<<sterictime
	         <<", electrostatic: "<<electrostatictime
	         <<", electrostaticfield: "<<electrostaticfieldtime<<", densityfield: "<<densityfieldtime<<", probe: "<<probetime<<", listes de cellules: "<<celllisttime<<", construction steric: "<<listbuildsterictime<<", construction coulomb: "<<listbuildcoulombtime<<", construction autre: "<<listbuildtime<<", impala: "<<impalatime
	         <<", hydrophobic: "<<hydrophobictime<<", hbond: "<<hbondtime
	         <<", damping: "<<dampingtime<<", integration: "<<integrationtime
	         <<", external: "<<externalforcetime<<" )"<<std::endl;
	if (getNeighborSkin() > 0.0f)
		biospring::logging::info("Neighbour lists rebuilt %u times over %d steps", _listrebuilds, getNbIterations());
	SpringNetwork::endRun();
	}


void SpringNetworkOpenCL::run()
	{
	initRun();
	_warnAboutTermsTheDeviceIgnores();
	while (!isEnd())
		{
		idleRun();

		// The same cadence as SpringNetwork::run(), and for the same reason:
		// without it a --opencl run printed its kernel timings and nothing
		// else, so neither the energies nor the framerate could be compared
		// against a CPU run of the same .msp.
		if (_isTimeToLogData())
			{
			_updateFrameRate();
			_displayFrameData();
			}
		}
	endRun();
	}


// The longest reach of any enabled pairwise term, in A. Zero when none is on.
float SpringNetworkOpenCL::_largestPairwiseCutoff() const
	{
	float longest = 0.0f;
	if (isStericEnabled())
		longest = std::max(longest, getStericCutoff());
	if (isElectrostaticCoulombEnabled())
		longest = std::max(longest, getElectrostaticCutoff());
	if (isHydrophobicityEnabled())
		longest = std::max(longest, getHydrophobicCutoff());
	return longest;
	}



// Refreshes every enabled term's grid.
//
// One grid per term, each with cells as wide as that term's own search radius
// (cutoff + skin), so every kernel reads a stencil of radius 1. A single shared
// grid with a per-term stencil radius was tried and reverted: it saves two
// rangings per step and loses more than that to the cells it makes every kernel
// visit. See getCellWidthFor().
//
// Nothing here if no pairwise term is on: the cheapest neighbour search is the
// one that does not run.
void SpringNetworkOpenCL::_updateCellLists()
	{
	if (!isStericEnabled() && !isElectrostaticCoulombEnabled() && !isHydrophobicityEnabled()
	    && !isHydrogenBondEnabled())
		return;

	// One grid per term, each with cells the size of that term's own search
	// radius, and each holding only the particles that term can interact with.
	// They are all rebuilt here, every step, so the cell-walk fallback never
	// reads a frame the list build happened to skip.
	_buildTermMasks();
	if (isStericEnabled())
		_buildCellList(_cells, getCellWidthFor(getStericCutoff()));
	if (isElectrostaticCoulombEnabled())
		_binSubsetIntoCells(_chargedcells, _masks.charged,
		                    getCellWidthFor(getElectrostaticCutoff()));
	if (isHydrophobicityEnabled())
		_binSubsetIntoCells(_hydrophobiccells, _masks.hydrophobic,
		                    getCellWidthFor(getHydrophobicCutoff()));
	// One grid for the hydrogen bond term, holding donors AND acceptors
	// together: a hydroxyl is both, and the kernels filter roles themselves.
	if (isHydrogenBondEnabled())
		_binSubsetIntoCells(_hydrogenbondcells, _masks.hydrogenbond,
		                    getCellWidthFor(getHydrogenBondCutoff()));
	}


std::vector<unsigned> SpringNetworkOpenCL::neighboursFromList(const NeighbourList & list, unsigned i)
	{
	std::vector<unsigned> neighbours;
	if (!list.valid || i + 1 > _nbparticlesocl)
		return neighbours;

	// The two bounds come from the device, because that is where the scan
	// leaves them. list.offsets, the host mirror, stopped being filled when
	// the prefix sum moved onto the card -- and nothing in the simulation
	// loop reads it any more. This accessor exists for the parity tests and
	// pays for its own data rather than making every step maintain a copy
	// nobody else wants.
	unsigned bounds[2] = {0u, 0u};
	_err = _queue.enqueueReadBuffer(list.offsetsbuffer, CL_TRUE, sizeof(unsigned) * i,
	                                sizeof(unsigned) * 2, bounds);
	checkErr("enqueueReadBuffer(list bounds)");
	const unsigned first = bounds[0];
	const unsigned last = bounds[1];
	if (last <= first)
		return neighbours;

	neighbours.resize(last - first);
	_err = _queue.enqueueReadBuffer(list.itemsbuffer, CL_TRUE, sizeof(unsigned) * first,
	                                sizeof(unsigned) * (last - first), neighbours.data());
	checkErr("enqueueReadBuffer(list items)");
	return neighbours;
	}


// Measures a subset grid its own frame at `width`, and bins into it only the
// particles the mask keeps.
//
// Its own frame, not a copy of _cells': the two hold different populations AND
// want different cell sizes, since a term's cells are the size of its own
// search radius. What it costs is a measuring pass and a head array per subset
// -- 1.2 MB on 034 -- and what it buys is that a Coulomb query walks past
// charged particles only, where it used to walk past every bead in the cell to
// reject four out of five.
void SpringNetworkOpenCL::_binSubsetIntoCells(CellGrid & subset, const std::vector<unsigned char> & mask,
                                              float width)
	{
	if (mask.empty() || width <= 0.0f)
		{
		subset.ncellstotal = 0;
		return;
		}

	// Its own frame, at its own width. _measureCellGrid allocates the head and
	// next arrays and their buffers, and invalidates the grid if it cannot
	// measure one, so a term that fails here sees "no grid" rather than last
	// step's cells over this step's positions.
	if (subset.ncellstotal == 0 || subset.requestedwidth != width || !_frameStillHolds(subset))
		{
		if (!_measureCellGrid(subset, width))
			return;
		}

	subset.includedbuffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	                                   sizeof(unsigned char) * _nbparticlesocl,
	                                   const_cast<unsigned char *>(mask.data()), &_err);
	checkErr("Buffer(subset mask)");
	subset.restricted = true;

	_binParticlesIntoCells(subset);
	}


// Whether the stored lists still describe the current positions.
//
// A pair absent from a list was farther than cutoff + skin when the list was
// built, so it cannot have reached the cutoff before the two closed by skin --
// and two particles each drifting skin/2 towards each other close exactly that.
// So the question is whether ANY particle has moved half the skin.
//
// Asked on the host over the positions that came back at the end of the last
// step, which cost 0.05 to 0.09 ms per step measured; the walk it avoids costs
// milliseconds.
bool SpringNetworkOpenCL::_listsNeedRebuilding()
	{
	if (!_listsarebuilt || _listreference.size() != _nbparticlesocl)
		return true;

	const float half = 0.5f * getNeighborSkin();
	const float halfsquared = half * half;
	for (unsigned i = 0; i < _nbparticlesocl; ++i)
		{
		// The probe moves under someone's hand, by as much as they like, and
		// it is not what the lists are about: its interactions are the probe
		// kernel's, pair by pair against every particle, with no cell list at
		// all. Counting its drift would rebuild every list at every step of an
		// interactive session for nothing. The CPU has always skipped it --
		// see nsearch.hpp's has_drifted -- and the device did not.
		if (SpringNetwork::isProbeParticle(i))
			continue;

		const float dx = _particlepositions[i].x - _listreference[i].x;
		const float dy = _particlepositions[i].y - _listreference[i].y;
		const float dz = _particlepositions[i].z - _listreference[i].z;
		const float moved = dx * dx + dy * dy + dz * dz;
		if (!(moved <= halfsquared))   // catches NaN too, which invalidates it
			return true;
		}
	return false;
	}


// Builds one term's list: count, pack, fill.
//
// Two passes over the same walk rather than one pass into a fixed capacity,
// because a capacity has to be sized for the densest particle and paid for by
// every other one. The prefix sum happens on the host: it is over N unsigned
// ints, a rebuild is rare by construction, and a device scan would be more
// code than it saves here.
void SpringNetworkOpenCL::_buildNeighbourList(NeighbourList & list, float cutoff,
                                              const std::vector<unsigned char> & targets,
                                              const std::vector<unsigned char> & candidates,
                                              const CellGrid & grid)
	{
	list.valid = false;
	if (_nbparticlesocl == 0 || grid.ncellstotal == 0 || cutoff <= 0.0f)
		return;

	list.radius = cutoff + getNeighborSkin();
	const auto upload = [&](const std::vector<unsigned char> & mask, cl::Buffer & buffer) {
		if (mask.empty())
			return false;
		buffer = cl::Buffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		                    sizeof(unsigned char) * _nbparticlesocl,
		                    const_cast<unsigned char *>(mask.data()), &_err);
		checkErr("Buffer(mask)");
		return true;
	};
	// Once, not once per step: _buildTermMasks clears the flag when the masks
	// are rebuilt, which only happens when the network changes size.
	if (!list.targetsuploaded)
		{
		list.hastargets = upload(targets, list.targetsbuffer);
		list.targetsuploaded = true;
		}
	list.hascandidates = upload(candidates, list.candidatesbuffer);

	// Sized once per particle count. The two host mirrors that used to live
	// beside these -- one N uints, one N+1 -- are gone with the prefix sum
	// they existed for: 900 kB on the capsid, across three lists, that nothing
	// read any more.
	if (list.buffersfor != _nbparticlesocl)
		{
		list.buffersfor = _nbparticlesocl;
		list.countsbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(unsigned) * _nbparticlesocl, NULL, &_err);
		checkErr("Buffer(counts)");
		// READ_WRITE, not READ_ONLY: the scan kernels write these. A kernel
		// writing to a read-only buffer does not fail -- the write simply does
		// not happen, and every neighbour list comes back empty.
		list.offsetsbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(unsigned) * (_nbparticlesocl + 1), NULL, &_err);
		checkErr("Buffer(offsets)");
		}

	const unsigned wg = WORK_GROUP_SIZE;
	const unsigned global = (_nbparticlesocl / wg) * wg + wg;

	const auto setwalkargs = [&](cl::Kernel & kernel) {
		unsigned a = 0;
		kernel.setArg(a++, _inoutPositionBuffer);
		kernel.setArg(a++, grid.headbuffer);
		kernel.setArg(a++, grid.nextbuffer);
		kernel.setArg(a++, grid.origin);
		kernel.setArg(a++, grid.width);
		kernel.setArg(a++, grid.ncells);
		kernel.setArg(a++, biospring_stencil_radius(list.radius, grid.width));
		if (list.hastargets)
			kernel.setArg(a++, list.targetsbuffer);
		else
			kernel.setArg(a++, sizeof(cl_mem), NULL);
		if (list.hascandidates)
			kernel.setArg(a++, list.candidatesbuffer);
		else
			kernel.setArg(a++, sizeof(cl_mem), NULL);
		kernel.setArg(a++, list.radius);
		return a;
	};

	unsigned a = setwalkargs(_kernelcountneighbours);
	_kernelcountneighbours.setArg(a++, list.countsbuffer);
	_kernelcountneighbours.setArg(a++, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kernelcountneighbours, cl::NullRange,
	                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(countneighbours)");
	_pendingevents.emplace_back(_event, listbuildinto);

	// The counts become offsets on the device: three kernels instead of a read
	// of the whole counts array, a serial pass over it, and a write back.
	// Four bytes still cross -- the grand total, which decides whether the
	// items buffer is big enough, and that is a host allocation.
	const unsigned nblocks = global / wg;
	if (list.blocksumsfor != nblocks)
		{
		list.blocksumsbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(unsigned) * nblocks, NULL, &_err);
		checkErr("Buffer(block sums)");
		list.totalbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(unsigned), NULL, &_err);
		checkErr("Buffer(scan total)");
		list.blocksumsfor = nblocks;
		}

	a = 0;
	_kernelscancounts.setArg(a++, list.countsbuffer);
	_kernelscancounts.setArg(a++, list.offsetsbuffer);
	_kernelscancounts.setArg(a++, list.blocksumsbuffer);
	_kernelscancounts.setArg(a++, cl::__local(sizeof(unsigned) * wg));
	_kernelscancounts.setArg(a++, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kernelscancounts, cl::NullRange,
	                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(scanCounts)");
	_pendingevents.emplace_back(_event, listbuildinto);

	a = 0;
	_kernelscanblocksums.setArg(a++, list.blocksumsbuffer);
	_kernelscanblocksums.setArg(a++, list.totalbuffer);
	_kernelscanblocksums.setArg(a++, nblocks);
	_err = _queue.enqueueNDRangeKernel(_kernelscanblocksums, cl::NullRange,
	                                   cl::NDRange(1), cl::NDRange(1), NULL, &_event);
	checkErr("enqueueNDRangeKernel(scanBlockSums)");
	_pendingevents.emplace_back(_event, listbuildinto);

	a = 0;
	_kerneladdblocksums.setArg(a++, list.offsetsbuffer);
	_kerneladdblocksums.setArg(a++, list.blocksumsbuffer);
	_kerneladdblocksums.setArg(a++, list.totalbuffer);
	_kerneladdblocksums.setArg(a++, wg);
	_kerneladdblocksums.setArg(a++, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kerneladdblocksums, cl::NullRange,
	                                   cl::NDRange(global + wg), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(addBlockSums)");
	_pendingevents.emplace_back(_event, listbuildinto);

	_err = _queue.enqueueReadBuffer(list.totalbuffer, CL_TRUE, 0, sizeof(unsigned), &list.total);
	checkErr("enqueueReadBuffer(scan total)");

	// Grown, never shrunk: the size a structure needs swings from step to step
	// and reallocating on every rebuild would cost more than the slack.
	if (list.total > list.capacity)
		{
		list.capacity = list.total + list.total / 8 + 1;
		list.itemsbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(unsigned) * list.capacity, NULL, &_err);
		checkErr("Buffer(items)");
		}

	if (list.total > 0)
		{
		a = setwalkargs(_kernelfillneighbours);
		_kernelfillneighbours.setArg(a++, list.offsetsbuffer);
		_kernelfillneighbours.setArg(a++, list.itemsbuffer);
		_kernelfillneighbours.setArg(a++, _nbparticlesocl);
		_err = _queue.enqueueNDRangeKernel(_kernelfillneighbours, cl::NullRange,
		                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
		checkErr("enqueueNDRangeKernel(fillneighbours)");
		_pendingevents.emplace_back(_event, listbuildinto);
		}

	list.valid = true;
	}


// Refreshes the stored neighbours of every enabled pairwise term.
//
// Nothing happens without a skin. A list built at the bare cutoff holds exactly
// this step's neighbours and is stale the moment anything moves, so it would be
// rebuilt every step -- two walks to save one, measured on the CPU at up to
// 2.4x slower on a structure whose fastest beads move 3 A in a step. So
// simulation.neighborskin = 0, the default, keeps the cell walk.
// Built once: nothing that feeds them moves during a run. isStatic, isCharged,
// isHydrophobic and the probe's index are all fixed once the network is loaded.
//
// TARGETS say who a term computes a force FOR. A static particle's force is
// never read -- the CPU loops over _dynamicparticules alone -- and an uncharged
// particle has no Coulomb force to compute, so neither needs a neighbourhood
// and neither needs a work item.
//
// CANDIDATES say who may appear in someone else's neighbourhood. A STATIC
// charged particle still pushes the dynamic ones, so it stays a candidate even
// though it is not a target. That asymmetry is the whole point of having two
// masks rather than one.
void SpringNetworkOpenCL::_buildTermMasks()
	{
	const bool wantcoulomb = isElectrostaticCoulombEnabled();
	const bool wanthydrophobic = isHydrophobicityEnabled();
	const bool wanthbond = isHydrogenBondEnabled();
	if (_masks.builtfor == _nbparticlesocl &&
	    _masks.dynamic.size() == _nbparticlesocl &&
	    _masks.charged.empty() != wantcoulomb &&
	    _masks.hydrophobic.empty() != wanthydrophobic &&
	    _masks.hydrogenbond.empty() != wanthbond)
		return;

	const unsigned n = std::min<unsigned>(SpringNetwork::getNumberOfParticles(), _nbparticlesocl);
	_masks.dynamic.assign(_nbparticlesocl, 0);
	_masks.dynamiccharged.clear();
	_masks.charged.clear();
	_masks.dynamichydrophobic.clear();
	_masks.hydrophobic.clear();
	_masks.hydrogenbond.clear();
	if (wantcoulomb)
		{ _masks.dynamiccharged.assign(_nbparticlesocl, 0); _masks.charged.assign(_nbparticlesocl, 0); }
	if (wanthydrophobic)
		{ _masks.dynamichydrophobic.assign(_nbparticlesocl, 0); _masks.hydrophobic.assign(_nbparticlesocl, 0); }
	if (wanthbond)
		_masks.hydrogenbond.assign(_nbparticlesocl, 0);

	for (unsigned i = 0; i < n; ++i)
		{
		const Particle & particle = SpringNetwork::getParticle(i);
		// The probe is nobody's neighbour: its interactions are the probe
		// kernel's, and the CPU keeps it out of every searcher for the same
		// reason.
		if (SpringNetwork::isProbeParticle(i))
			continue;

		const bool isdynamic = particle.isDynamic();
		_masks.dynamic[i] = isdynamic ? 1 : 0;
		if (wantcoulomb && particle.isCharged())
			{
			_masks.charged[i] = 1;
			_masks.dynamiccharged[i] = isdynamic ? 1 : 0;
			}
		if (wanthydrophobic && particle.isHydrophobic())
			{
			_masks.hydrophobic[i] = 1;
			_masks.dynamichydrophobic[i] = isdynamic ? 1 : 0;
			}
		// Donors and acceptors are both targets AND candidates here, unlike
		// every other term: a bond needs a free slot at each end, and which
		// end donates is decided pair by pair rather than by particle. A
		// STATIC donor still anchors a bond, so it is not filtered out.
		if (wanthbond && (particle.donorCapacity() > 0 || particle.acceptorCapacity() > 0))
			_masks.hydrogenbond[i] = 1;
		}

	_masks.builtfor = _nbparticlesocl;
	_stericlist.targetsuploaded = false;
	_electrostaticlist.targetsuploaded = false;
	_hydrophobiclist.targetsuploaded = false;
	_hydrogenbondlist.targetsuploaded = false;
	}


void SpringNetworkOpenCL::_updateNeighbourLists()
	{
	if (!isStericEnabled() && !isElectrostaticCoulombEnabled() && !isHydrophobicityEnabled()
	    && !isHydrogenBondEnabled())
		{
		_stericlist.valid = _electrostaticlist.valid = _hydrophobiclist.valid = false;
		_hydrogenbondlist.valid = false;
		return;
		}

	// A list beats walking the cells -- by 34% on 023 and 45% on 034 before any
	// subset filtering, 62% and 92% with it -- because the force kernel stops
	// chasing `nextincell[p]`, whose every load waits for the previous one to
	// return, and reads a contiguous run of indices instead.
	//
	// It IS a Verlet list and there IS something to amortise, which is the part
	// this used to deny. Built at cutoff + skin, it survives until something
	// moves half the skin, and at the 1.0 A default that is 165 rebuilds in 1000
	// steps on 023 and 9 on the capsid. That matters because the build is not
	// cheap: countNeighbours + fillNeighbours are 45% of all device time on 023
	// and 58% on 024 when they run every step. See defaultNeighborSkin().
	//
	// _listsNeedRebuilding() is what keeps it exact -- no pair can be missed,
	// and an interactive pull that moves a bead far triggers a rebuild like any
	// other motion.
	if (getNeighborSkin() > 0.0f && !_listsNeedRebuilding())
		return;

	// Two masks per term, because the two ends of a pair are not the same
	// question -- see the kernels.
	//
	// TARGETS: who gets a list. A static particle's force is never read -- the
	// CPU only loops over _dynamicparticules -- so finding its neighbours is
	// work thrown away, and several examples here are 89% to 100% static.
	// Coulomb narrows it further to the charged ones.
	//
	// CANDIDATES: who may appear in a list. A static charged particle still
	// pushes the dynamic ones, so it stays a candidate. For the steric term
	// that is everybody, which is what an empty mask means.
	if (isStericEnabled())
		{
		// Everyone is a candidate, so the term's own grid already holds them all.
		listbuildinto = &listbuildsterictime;
		_buildNeighbourList(_stericlist, getStericCutoff(), _masks.dynamic, std::vector<unsigned char>(), _cells);
		}
	if (isElectrostaticCoulombEnabled())
		{
		// Its grid holds the charged particles only -- built in
		// _updateCellLists -- so the candidate mask is redundant here: there is
		// nothing else in there to reject, and dropping it takes a test out of
		// the innermost loop.
		listbuildinto = &listbuildcoulombtime;
		_buildNeighbourList(_electrostaticlist, getElectrostaticCutoff(), _masks.dynamiccharged,
		                    std::vector<unsigned char>(), _chargedcells);
		}
	listbuildinto = &listbuildtime;
	if (isHydrophobicityEnabled())
		{
		_buildNeighbourList(_hydrophobiclist, getHydrophobicCutoff(), _masks.dynamichydrophobic,
		                    std::vector<unsigned char>(), _hydrophobiccells);
		}
	// Targets are donors and acceptors alike, STATIC ONES INCLUDED: a static
	// donor still holds a bond and still pulls the dynamic partner, and the
	// assignment is symmetric, so restricting targets to the dynamic ones
	// would lose every bond a static side happens to donate.
	if (isHydrogenBondEnabled())
		{
		_buildNeighbourList(_hydrogenbondlist, getHydrogenBondCutoff(), _masks.hydrogenbond,
		                    std::vector<unsigned char>(), _hydrogenbondcells);
		}

	_listreference.assign(_particlepositions, _particlepositions + _nbparticlesocl);
	_listsarebuilt = true;
	_listrebuilds++;
	}


// Sizes the grid to the box the particles currently occupy.
//
// THIS is what "building the grid" means, and it is not what happens every
// step. The cells, their width and the origin they are counted from are a
// fixed frame; what moves is which cell each particle is in, and that is
// _binParticlesIntoCells below. Conflating the two costs a pass over every
// position, every step, to answer a question that is almost always no.
//
// Remeasured rather than fixed for the whole run, though: this is not a
// periodic simulation, the structure translates and swells, and a particle
// outside the grid is invisible to every neighbour walk. The device says when
// that has happened (see binParticles), so the pass below runs on the first
// step and then only when someone has actually left.
bool SpringNetworkOpenCL::_measureCellGrid(CellGrid & grid, float requestedwidth)
	{
	if (_nbparticlesocl == 0 || requestedwidth <= 0.0f)
		return false;

	float lo[3] = {_particlepositions[0].x, _particlepositions[0].y, _particlepositions[0].z};
	float hi[3] = {lo[0], lo[1], lo[2]};
	for (unsigned i = 1; i < _nbparticlesocl; ++i)
		{
		const float p[3] = {_particlepositions[i].x, _particlepositions[i].y, _particlepositions[i].z};
		for (int d = 0; d < 3; ++d)
			{
			if (!std::isfinite(p[d]))
				{
				// _syncParticlesFromDevice reports this properly, with the
				// particle's id, at the end of the step. Here there is nothing
				// to measure: invalidate rather than leave a frame that no
				// longer describes anything.
				grid.ncellstotal = 0;
				return false;
				}
			if (p[d] < lo[d]) lo[d] = p[d];
			if (p[d] > hi[d]) hi[d] = p[d];
			}
		}

	// Margin on each side, in angstroms rather than in cells. It is what buys
	// the frame its lifetime: a particle outside the grid is binned nowhere and
	// is invisible to every walk, so the device reports the first one to leave
	// and the pass above runs again. The stencil itself needs no margin -- it is
	// bounds-checked, and a cell outside the grid holds nothing by construction.
	//
	// Two cells, which is this grid's own term reaching twice. It used to be
	// twice the LONGEST pairwise cutoff of any term, from when one grid served
	// them all; a 6 A steric grid then carried a 16 A term's headroom and six
	// cells of empty margin on every side.
	const float MARGIN = 2.0f * requestedwidth;
	const long CELL_LIMIT = 8L * 1024L * 1024L;

	// WIDER cells than asked for are still correct: every walk derives its
	// stencil radius from the width it is given, so it simply takes fewer,
	// bigger steps and the distance test drops the surplus. Only cells narrower
	// than the caller believes would be wrong. So a box too big for the cell
	// count is answered by widening rather than by refusing to build -- which
	// degrades towards brute force, slowly and correctly, instead of leaving
	// every term with no neighbour structure and no way to know it.
	//
	// It is not reached by a healthy structure: the largest example here is
	// 034's capsid, a 292 A cube, which at 3 A cells asks for about a million
	// against the 8M below. It is reached by a diverging one, and then saying so
	// is worth more than the grid.
	float width = requestedwidth;
	cl_int4 ncells;
	long total = 0;
	for (int attempt = 0; attempt < 64; ++attempt)
		{
		const int margincells = static_cast<int>(std::ceil(MARGIN / width));
		total = 1;
		bool overflowed = false;
		for (int d = 0; d < 3; ++d)
			{
			const double span = (hi[d] - lo[d]) / width;
			if (!(span < static_cast<double>(CELL_LIMIT)))
				{
				overflowed = true;
				break;
				}
			const int n = static_cast<int>(std::floor(span)) + 1 + 2 * margincells;
			ncells.s[d] = n;
			total *= n;
			if (total > CELL_LIMIT)
				{
				overflowed = true;
				break;
				}
			}
		if (!overflowed)
			{
			for (int d = 0; d < 3; ++d)
				grid.origin.s[d] = lo[d] - margincells * width;
			break;
			}
		width *= 2.0f;
		total = 0;
		}

	if (total <= 0)
		{
		// 64 doublings and still too big means the coordinates are not a
		// structure any more.
		BIOSPRING_WARN_ONCE("OpenCL cell list: no grid fits the coordinates (extent %.3g x %.3g x %.3g A); "
		                    "the structure has almost certainly diverged.",
		                    hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
		grid.ncellstotal = 0;
		return false;
		}

	if (width != requestedwidth)
		BIOSPRING_WARN_ONCE("OpenCL cell list: the structure's box needs cells of %.2f A rather than the "
		                    "%.2f A asked for, to stay under %ld cells. Still exact, and slower: "
		                    "each walk sifts (%.1f)^3 times as many candidates.",
		                    width, requestedwidth, CELL_LIMIT, width / requestedwidth);

	ncells.s[3] = 0;
	grid.origin.s[3] = 0.0f;
	grid.requestedwidth = requestedwidth;
	grid.width = width;
	grid.ncells = ncells;

	const unsigned ncellstotal = static_cast<unsigned>(total);
	if (ncellstotal != grid.ncellstotal)
		{
		delete[] grid.head;
		grid.head = new unsigned[ncellstotal];
		grid.ncellstotal = ncellstotal;
		grid.headbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		                             sizeof(unsigned) * ncellstotal, grid.head, &_err);
		checkErr("Buffer::Buffer(cellhead)");
		}
	if (grid.next == nullptr)
		{
		grid.next = new unsigned[_nbparticlesocl];
		grid.nextbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
		                               sizeof(unsigned) * _nbparticlesocl, grid.next, &_err);
		checkErr("Buffer::Buffer(nextincell)");
		}

	return true;
	}


// Puts every particle in the cell it is in now. This is the per-step work, and
// all of it: two kernel launches over a frame that does not move.
void SpringNetworkOpenCL::_binParticlesIntoCells(CellGrid & grid)
	{
	const unsigned wg = WORK_GROUP_SIZE;

	const unsigned cellglobal = (grid.ncellstotal / wg) * wg + wg;
	_kernelblankcells.setArg(0, grid.ncellstotal);
	_kernelblankcells.setArg(1, grid.headbuffer);
	_err = _queue.enqueueNDRangeKernel(_kernelblankcells, cl::NullRange,
	                                   cl::NDRange(cellglobal), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(blankCells)");
	_pendingevents.emplace_back(_event, &celllisttime);

	const unsigned partglobal = (_nbparticlesocl / wg) * wg + wg;
	_kernelbinparticles.setArg(0, _inoutPositionBuffer);
	_kernelbinparticles.setArg(1, grid.origin);
	_kernelbinparticles.setArg(2, grid.width);
	_kernelbinparticles.setArg(3, grid.ncells);
	_kernelbinparticles.setArg(4, grid.headbuffer);
	_kernelbinparticles.setArg(5, grid.nextbuffer);
	if (grid.restricted)
		_kernelbinparticles.setArg(6, grid.includedbuffer);
	else
		_kernelbinparticles.setArg(6, sizeof(cl_mem), NULL);
	_kernelbinparticles.setArg(7, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kernelbinparticles, cl::NullRange,
	                                   cl::NDRange(partglobal), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(binParticles)");
	_pendingevents.emplace_back(_event, &celllisttime);
	}


// Whether the frame still contains every particle.
//
// Asked on the HOST, over the positions _syncParticlesFromDevice has already
// read back, because that is free and asking the device is not. The first
// version of this had binParticles raise a flag in a one-int buffer, which the
// host then read: correct, and slower than what it replaced -- two blocking
// four-byte transfers per grid per step cost more than a pass over 2811
// positions already in cache. Measured on 072: 0.40 ms a step per grid that
// way against 0.22 ms with the pass.
bool SpringNetworkOpenCL::_frameStillHolds(const CellGrid & grid) const
	{
	if (grid.ncellstotal == 0)
		return false;

	for (unsigned i = 0; i < _nbparticlesocl; ++i)
		{
		const float p[3] = {_particlepositions[i].x, _particlepositions[i].y, _particlepositions[i].z};
		for (int d = 0; d < 3; ++d)
			{
			const float local = (p[d] - grid.origin.s[d]) / grid.width;
			if (!(local >= 0.0f) || local >= static_cast<float>(grid.ncells.s[d]))
				return false;   // outside, or not a number at all
			}
		}
	return true;
	}


// Re-places the particles in the grid, measuring a new frame only when the one
// in hand no longer holds them.
bool SpringNetworkOpenCL::_buildCellList(CellGrid & grid, float width)
	{
	// A first call, a width that changed under us -- the .msp can be reloaded
	// mid-run, and a cutoff moving moves the width with it -- or a structure
	// that has outgrown its frame.
	if (grid.ncellstotal == 0 || grid.requestedwidth != width || !_frameStillHolds(grid))
		{
		// _measureCellGrid invalidates the grid when it cannot measure one, so
		// a term that consults it afterwards sees "no grid" rather than last
		// step's cells over this step's positions -- which is the one outcome
		// worse than having no grid at all, being wrong without saying so.
		if (!_measureCellGrid(grid, width))
			return false;
		}

	_binParticlesIntoCells(grid);
	return true;
	}


// The stencil around a particle's own cell, at the radius `cutoff` asks of the
// grid's width, filtered to `cutoff`.
//
// A stencil is a box and a cutoff is a ball, so the stencil is necessary and
// not sufficient: everything within the cutoff is in one of its cells, and some
// of what is in them is beyond it. The distance test is what makes the answer
// exact, and a force kernel has to make the same one -- which is why this walks
// the cells exactly as the kernels do, down to the corner test.
std::vector<unsigned> SpringNetworkOpenCL::neighborsFromCellList(const CellGrid & grid, unsigned i,
                                                                float cutoff)
	{
	std::vector<unsigned> neighbors;
	if (grid.ncellstotal == 0 || i >= _nbparticlesocl)
		return neighbors;

	// The device owns these between rebuilds; read them rather than trusting
	// the host mirror to be coherent.
	_err = _queue.enqueueReadBuffer(grid.headbuffer, CL_TRUE, 0,
	                                sizeof(unsigned) * grid.ncellstotal, grid.head);
	checkErr("enqueueReadBuffer(cellhead)");
	_err = _queue.enqueueReadBuffer(grid.nextbuffer, CL_TRUE, 0,
	                                sizeof(unsigned) * _nbparticlesocl, grid.next);
	checkErr("enqueueReadBuffer(nextincell)");

	const float4 & here = _particlepositions[i];
	const int cx = static_cast<int>(std::floor((here.x - grid.origin.s[0]) / grid.width));
	const int cy = static_cast<int>(std::floor((here.y - grid.origin.s[1]) / grid.width));
	const int cz = static_cast<int>(std::floor((here.z - grid.origin.s[2]) / grid.width));

	const float cutoffsq = cutoff * cutoff;
	const int k = biospring_stencil_radius(cutoff, grid.width);
	for (int dz = -k; dz <= k; ++dz)
		for (int dy = -k; dy <= k; ++dy)
			for (int dx = -k; dx <= k; ++dx)
				{
				if (!biospring_cell_in_range(dx, dy, dz, grid.width, cutoffsq))
					continue;

				const int x = cx + dx, y = cy + dy, z = cz + dz;
				// No periodicity: a stencil cell outside the grid is simply not
				// there, never the cell on the opposite face.
				if (x < 0 || y < 0 || z < 0 || x >= grid.ncells.s[0] || y >= grid.ncells.s[1] || z >= grid.ncells.s[2])
					continue;

				const unsigned cell = static_cast<unsigned>((z * grid.ncells.s[1] + y) * grid.ncells.s[0] + x);
				for (unsigned p = grid.head[cell]; p != EMPTY_CELL; p = grid.next[p])
					{
					if (p == i)
						continue;
					const float ddx = _particlepositions[p].x - here.x;
					const float ddy = _particlepositions[p].y - here.y;
					const float ddz = _particlepositions[p].z - here.z;
					if (ddx * ddx + ddy * ddy + ddz * ddz <= cutoffsq)
						neighbors.push_back(p);
					}
				}

	return neighbors;
	}


// Flattens the torsions and builds the CSR from each particle to the torsions
// it takes part in.
//
// The tables are stored one after another, each bins + 1 samples long, so a
// torsion's table is at torsiontable[ti] * (bins + 1). Every table in a file
// shares its bin count -- the generator emits them together -- which is what
// lets one stride serve them all.
void SpringNetworkOpenCL::computeOpenCLTorsions()
	{
	delete[] _torsionatoms;   _torsionatoms = nullptr;
	delete[] _torsiontable;   _torsiontable = nullptr;
	delete[] _torsionfamily;  _torsionfamily = nullptr;
	delete[] _torsiontableenergy; _torsiontableenergy = nullptr;
	delete[] _torsiontabletorque; _torsiontabletorque = nullptr;
	delete[] _torsionoffsets; _torsionoffsets = nullptr;
	delete[] _torsionentries; _torsionentries = nullptr;

	const std::vector<Torsion> & torsions = getTorsions();
	const std::vector<TorsionTable> & tables = getTorsionTables();
	_nbtorsionsocl = static_cast<unsigned>(torsions.size());
	_nbtorsionentries = 0;
	_torsionbins = tables.empty() ? 0 : tables[0].bins;
	if (_nbtorsionsocl == 0 || tables.empty() || _torsionbins == 0)
		return;

	_torsionatoms = new cl_uint4[_nbtorsionsocl];
	_torsiontable = new unsigned[_nbtorsionsocl];
	_torsionfamily = new unsigned[_nbtorsionsocl];
	for (unsigned i = 0; i < _nbtorsionsocl; ++i)
		{
		_torsionatoms[i].s[0] = torsions[i].atoms[0];
		_torsionatoms[i].s[1] = torsions[i].atoms[1];
		_torsionatoms[i].s[2] = torsions[i].atoms[2];
		_torsionatoms[i].s[3] = torsions[i].atoms[3];
		_torsiontable[i] = torsions[i].table;
		_torsionfamily[i] = torsions[i].family;
		}

	const unsigned samples = _torsionbins + 1;
	_torsiontableenergy = new float[tables.size() * samples];
	_torsiontabletorque = new float[tables.size() * samples];
	for (size_t t = 0; t < tables.size(); ++t)
		{
		// A table shorter than the stride would make the kernel read the next
		// one's first samples as its own last, which is a wrong force and not a
		// crash. Refuse rather than truncate.
		if (tables[t].bins != _torsionbins || tables[t].energy.size() < samples
		    || tables[t].torque.size() < samples)
			{
			biospring::logging::warning("OpenCL torsions: table %zu has %u bins against %u, or too few "
			                            "samples; torsions are not evaluated on the device.",
			                            t, tables[t].bins, _torsionbins);
			_nbtorsionsocl = 0;
			return;
			}
		for (unsigned k = 0; k < samples; ++k)
			{
			_torsiontableenergy[t * samples + k] = tables[t].energy[k];
			_torsiontabletorque[t * samples + k] = tables[t].torque[k];
			}
		}

	// The CSR: count, then fill.
	std::vector<unsigned> counts(_nbparticlesocl, 0u);
	for (unsigned i = 0; i < _nbtorsionsocl; ++i)
		for (unsigned k = 0; k < 4; ++k)
			{
			const unsigned a = torsions[i].atoms[k];
			if (a < _nbparticlesocl)
				counts[a]++;
			}

	_torsionoffsets = new int[_nbparticlesocl + 1];
	unsigned running = 0;
	for (unsigned p = 0; p < _nbparticlesocl; ++p)
		{
		_torsionoffsets[p] = static_cast<int>(running);
		running += counts[p];
		}
	_torsionoffsets[_nbparticlesocl] = static_cast<int>(running);
	_nbtorsionentries = running;

	_torsionentries = new unsigned[_nbtorsionentries == 0 ? 1 : _nbtorsionentries];
	std::vector<unsigned> cursor(_nbparticlesocl, 0u);
	for (unsigned i = 0; i < _nbtorsionsocl; ++i)
		for (unsigned k = 0; k < 4; ++k)
			{
			const unsigned a = torsions[i].atoms[k];
			if (a >= _nbparticlesocl)
				continue;
			const unsigned at = static_cast<unsigned>(_torsionoffsets[a]) + cursor[a]++;
			// (torsion << 2 | slot): the slot says which of the four forces is
			// this particle's share.
			_torsionentries[at] = (i << 2) | k;
			}
	}


// The eight family switches as the bitmask the kernel takes, in the order of
// SpringNetwork::DihedralFamilyIndex -- the same order computeTorsionForces
// builds its `enabled` array in.
int SpringNetworkOpenCL::_torsionFamilyMask() const
	{
	int mask = 0;
	const bool on[] = {isDihedralPhiEnabled(),          isDihedralPsiEnabled(),
	                   isDihedralOmegaEnabled(),        isDihedralChiEnabled(),
	                   isDihedralPlanarityEnabled(),    isDihedralNucleicBackboneEnabled(),
	                   isDihedralNucleicChiEnabled(),   isDihedralNucleicSugarEnabled()};
	for (int i = 0; i < 8; ++i)
		if (on[i])
			mask |= (1 << i);
	return mask;
	}


// The kinetic energy of the velocities the device just returned, by the
// formula Particle::_integrateForce uses -- and from the same side of the
// integration as the CPU, which writes it there.
void SpringNetworkOpenCL::_computeEnergiesFromDeviceState()
	{
	if (!_measuringthisstep)
		return;   // nothing will read them; see idleRun

	// Summed here rather than reduced on the device: 2811 floats is a fraction
	// of the transfer that brought them, so a reduction kernel would save
	// nothing worth a second launch.
	_energies.spring = 0.0f;
	_energies.dihedral = 0.0f;
	if (isSpringEnabled() && _nbspringsocl > 0)
		{
		_err = _queue.enqueueReadBuffer(_springEnergyBuffer, CL_TRUE, 0,
		                                sizeof(float) * _nbparticlesocl, _springenergyper);
		checkErr("enqueueReadBuffer(spring energy)");
		for (unsigned i = 0; i < _nbparticlesocl; ++i)
			_energies.spring += _springenergyper[i];
		}
	if (isSpringEnabled() && _nbtorsionsocl > 0)
		{
		_err = _queue.enqueueReadBuffer(_torsionEnergyBuffer, CL_TRUE, 0,
		                                sizeof(float) * _nbparticlesocl, _torsionenergyper);
		checkErr("enqueueReadBuffer(torsion energy)");
		for (unsigned i = 0; i < _nbparticlesocl; ++i)
			_energies.dihedral += _torsionenergyper[i];
		}

	// Summed on the host like the spring energy above, and for the same
	// reason: the transfer that brings it back costs more than the sum.
	_energies.hbond = 0.0f;
	if (isHydrogenBondEnabled() && _hbond.uploaded && _hbondenergyper != nullptr)
		{
		_err = _queue.enqueueReadBuffer(_hbond.energybuffer, CL_TRUE, 0,
		                                sizeof(float) * _nbparticlesocl, _hbondenergyper);
		checkErr("enqueueReadBuffer(hydrogen bond energy)");
		for (unsigned i = 0; i < _nbparticlesocl; ++i)
			_energies.hbond += _hbondenergyper[i];
		}

	float kinetic = 0.0f;
	for (size_t i = 0; i < _dynamicparticules.size(); ++i)
		{
		const Particle & p = getParticle(_dynamicparticules[i]);
		const float v = p.getVelocity().norm();
		kinetic += 0.5f * p.getMass() * v * v *
		           biospring::forcefield::GLOBAL_KINETIC_ENERGY_CONVERT;
		}
	_energies.kinetic = kinetic;
	}


// Uploaded once. Capacities, antecedents, residue numbers and chain identity
// are all fixed when the network is loaded, and the slots themselves start
// empty -- the base class sized them in setup(), so this only has to move them
// across.
//
// Chain NAMES are strings, which a kernel cannot compare. They are mapped to
// indices here, in first-appearance order: the kernel only ever tests two
// particles for the same chain, so any injective mapping does.
// The one transfer the thermostatted path adds. BAOAB drifts before anything
// is evaluated, and the grid frame, the list rebuild criterion and the
// bounding box are all host loops over _particlepositions, so they have to see
// the drifted ones or they would be a full step behind what the force kernels
// read.
void SpringNetworkOpenCL::_readPositionsBack()
	{
	_err = _queue.enqueueReadBuffer(_inoutPositionBuffer, CL_TRUE, 0,
	                                sizeof(float4) * _nbparticlesocl, _particlepositions);
	checkErr("enqueueReadBuffer(positions after drift)");
	}


void SpringNetworkOpenCL::_uploadHydrogenBondTopology()
	{
	if (_hbond.uploaded || _nbparticlesocl == 0)
		return;

	const unsigned n = std::min<unsigned>(SpringNetwork::getNumberOfParticles(), _nbparticlesocl);

	std::vector<unsigned> donoroffset(_nbparticlesocl + 1, 0);
	std::vector<unsigned> acceptoroffset(_nbparticlesocl + 1, 0);
	std::vector<cl_int2> antecedent(_nbparticlesocl);
	std::vector<int> resid(_nbparticlesocl, -1);
	std::vector<int> chain(_nbparticlesocl, -1);

	std::map<std::string, int> chainindex;
	for (unsigned i = 0; i < _nbparticlesocl; ++i)
		{
		antecedent[i].s[0] = -1;
		antecedent[i].s[1] = -1;
		if (i >= n)
			{
			donoroffset[i + 1] = donoroffset[i];
			acceptoroffset[i + 1] = acceptoroffset[i];
			continue;
			}
		const Particle & p = SpringNetwork::getParticle(i);
		donoroffset[i + 1] = donoroffset[i] + static_cast<unsigned>(p.donorCapacity());
		acceptoroffset[i + 1] = acceptoroffset[i] + static_cast<unsigned>(p.acceptorCapacity());
		antecedent[i].s[0] = p.antecedentIndex();
		antecedent[i].s[1] = p.antecedentIndex2();
		resid[i] = static_cast<int>(p.getResId());
		const auto inserted = chainindex.emplace(p.getChainName(), static_cast<int>(chainindex.size()));
		chain[i] = inserted.first->second;
		}

	const unsigned totaldonor = donoroffset[_nbparticlesocl];
	const unsigned totalacceptor = acceptoroffset[_nbparticlesocl];
	// Every slot starts free. A zero-sized buffer is rejected by OpenCL, so a
	// network with no donor at all still gets one entry that nothing reads.
	std::vector<int> donorslot(std::max(1u, totaldonor), -1);
	std::vector<int> acceptorslot(std::max(1u, totalacceptor), -1);

	const auto upload = [&](cl::Buffer & buffer, const void * data, size_t bytes, cl_mem_flags flags) {
		buffer = cl::Buffer(_context, flags | CL_MEM_COPY_HOST_PTR, bytes, const_cast<void *>(data), &_err);
		checkErr("Buffer(hydrogen bond)");
	};
	upload(_hbond.donoroffsetbuffer, donoroffset.data(), sizeof(unsigned) * donoroffset.size(), CL_MEM_READ_ONLY);
	upload(_hbond.acceptoroffsetbuffer, acceptoroffset.data(), sizeof(unsigned) * acceptoroffset.size(), CL_MEM_READ_ONLY);
	upload(_hbond.donorslotbuffer, donorslot.data(), sizeof(int) * donorslot.size(), CL_MEM_READ_WRITE);
	upload(_hbond.acceptorslotbuffer, acceptorslot.data(), sizeof(int) * acceptorslot.size(), CL_MEM_READ_WRITE);
	upload(_hbond.antecedentbuffer, antecedent.data(), sizeof(cl_int2) * antecedent.size(), CL_MEM_READ_ONLY);
	upload(_hbond.residbuffer, resid.data(), sizeof(int) * resid.size(), CL_MEM_READ_ONLY);
	upload(_hbond.chainbuffer, chain.data(), sizeof(int) * chain.size(), CL_MEM_READ_ONLY);

	_hbond.nearestbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(int) * _nbparticlesocl, NULL, &_err);
	checkErr("Buffer(hbond nearest)");
	_hbond.strengthbuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(float) * _nbparticlesocl, NULL, &_err);
	checkErr("Buffer(hbond strength)");
	_hbond.energybuffer = cl::Buffer(_context, CL_MEM_READ_WRITE, sizeof(float) * _nbparticlesocl, NULL, &_err);
	checkErr("Buffer(hbond energy)");
	delete[] _hbondenergyper;
	_hbondenergyper = new float[_nbparticlesocl];

	biospring::logging::info("Hydrogen bond slots on device: %u donor, %u acceptor.", totaldonor, totalacceptor);
	_hbond.uploaded = true;
	}


// The CPU's four rounds, as kernels. The host neither sees nor decides any of
// it -- which is the whole point: an assignment computed on the host would have
// to come down, be walked serially, and go back up, on the critical path of
// every step.
//
// Four rounds fixed rather than "until a round confirms nothing". The CPU stops
// early by counting confirmations, which on the device would mean reading a
// counter back between rounds -- four round trips to save at most two kernel
// launches that cost microseconds. One round fills one slot per particle, and
// the largest capacity in the data is two, so four is already generous.
void SpringNetworkOpenCL::_assignHydrogenBondPairsOnDevice()
	{
	const unsigned wg = WORK_GROUP_SIZE;
	const unsigned global = (_nbparticlesocl / wg) * wg + wg;
	const float cutoff = getHydrogenBondCutoff();
	const int probeid = isProbeEnabled() ? static_cast<int>(SpringNetwork::getNumberOfParticles()) : -1;

	unsigned a = 0;
	_kernelhbondbreak.setArg(a++, _inoutPositionBuffer);
	_kernelhbondbreak.setArg(a++, _hbond.donoroffsetbuffer);
	_kernelhbondbreak.setArg(a++, _hbond.donorslotbuffer);
	_kernelhbondbreak.setArg(a++, _hbond.acceptoroffsetbuffer);
	_kernelhbondbreak.setArg(a++, _hbond.acceptorslotbuffer);
	_kernelhbondbreak.setArg(a++, cutoff);
	_kernelhbondbreak.setArg(a++, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kernelhbondbreak, cl::NullRange,
	                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
	checkErr("enqueueNDRangeKernel(hbondBreak)");
	_pendingevents.emplace_back(_event, &hbondtime);

	if (!_hydrogenbondlist.valid)
		return;   // no list, no candidates to propose from

	for (unsigned round = 0; round < 4; ++round)
		{
		a = 0;
		_kernelhbondscore.setArg(a++, _inoutPositionBuffer);
		_kernelhbondscore.setArg(a++, _hbond.donoroffsetbuffer);
		_kernelhbondscore.setArg(a++, _hbond.donorslotbuffer);
		_kernelhbondscore.setArg(a++, _hbond.acceptoroffsetbuffer);
		_kernelhbondscore.setArg(a++, _hbond.acceptorslotbuffer);
		_kernelhbondscore.setArg(a++, _hbond.antecedentbuffer);
		_kernelhbondscore.setArg(a++, _hbond.residbuffer);
		_kernelhbondscore.setArg(a++, _hbond.chainbuffer);
		_kernelhbondscore.setArg(a++, _hydrogenbondlist.offsetsbuffer);
		_kernelhbondscore.setArg(a++, _hydrogenbondlist.itemsbuffer);
		_kernelhbondscore.setArg(a++, _inSpringBuffer);
		_kernelhbondscore.setArg(a++, _inSpringIndexesBuffer);
		_kernelhbondscore.setArg(a++, static_cast<int>(isSpringEnabled() && _nbspringsocl > 0));
		_kernelhbondscore.setArg(a++, probeid);
		_kernelhbondscore.setArg(a++, cutoff);
		_kernelhbondscore.setArg(a++, getForceField()->getHydrogenBondWellDepth());
		_kernelhbondscore.setArg(a++, getForceField()->getHydrogenBondEquilibrium());
		_kernelhbondscore.setArg(a++, getForceField()->getHydrogenBondWidth());
		_kernelhbondscore.setArg(a++, getForceField()->getHydrogenBondScale());
		_kernelhbondscore.setArg(a++, _hbond.nearestbuffer);
		_kernelhbondscore.setArg(a++, _hbond.strengthbuffer);
		_kernelhbondscore.setArg(a++, _nbparticlesocl);
		_err = _queue.enqueueNDRangeKernel(_kernelhbondscore, cl::NullRange,
		                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
		checkErr("enqueueNDRangeKernel(hbondScore)");
		_pendingevents.emplace_back(_event, &hbondtime);

		a = 0;
		_kernelhbondconfirm.setArg(a++, _hbond.nearestbuffer);
		_kernelhbondconfirm.setArg(a++, _hbond.donoroffsetbuffer);
		_kernelhbondconfirm.setArg(a++, _hbond.donorslotbuffer);
		_kernelhbondconfirm.setArg(a++, _hbond.acceptoroffsetbuffer);
		_kernelhbondconfirm.setArg(a++, _hbond.acceptorslotbuffer);
		_kernelhbondconfirm.setArg(a++, _nbparticlesocl);
		_err = _queue.enqueueNDRangeKernel(_kernelhbondconfirm, cl::NullRange,
		                                   cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
		checkErr("enqueueNDRangeKernel(hbondConfirm)");
		_pendingevents.emplace_back(_event, &hbondtime);
		}
	}


void SpringNetworkOpenCL::_warnAboutTermsTheDeviceIgnores() const
	{
	// biospring.cl evaluates springs, steric, Coulomb, damping, external forces
	// and the integration. Everything else a .msp can switch on is not computed
	// on this path, and used to be so without a word -- the run looked like the
	// CPU's and was a different model.
	std::string ignored;
	const auto add = [&ignored](const char * name) {
		if (!ignored.empty())
			ignored += ", ";
		ignored += name;
	};
	// Evaluated now, except in the two cases the kernel does not cover.
	if (isIMPEnabled() && !_membraneIsFlat())
		add("impala (the membrane is curved or doubled: a different model)");
	if (isInsertionVectorEnabled()) add("insertionvector");
	if (isRigidBodyEnabled())       add("rigidbody");
	// biospring.cl has no hydrogen bond code at all: not the Morse well, not
	// the two-sided angular weight, not the per-step re-pairing. A .msp with
	// hbond.enable = 1 run through --opencl is therefore a DIFFERENT model,
	// and was silently so until this line -- the very thing this warning was
	// written to stop.
	// hbond is no longer listed: biospring.cl now carries the Morse well, the
	// two-sided angular weight, the per-step re-pairing AND the core
	// repulsion.

	if (!ignored.empty())
		biospring::logging::warning(
		    "OpenCL backend: %s enabled in the configuration but NOT evaluated on the device -- "
		    "it computes springs, steric, Coulomb, hydrophobicity, damping, external forces and "
		    "the integration, nothing else. The reported energies cover only what it evaluated.",
		    ignored.c_str());
	}


// Only the terms biospring.cl evaluates, so that an enabled-but-unevaluated
// term reads as absent rather than as a measured zero.
void SpringNetworkOpenCL::_displayFrameData()
	{
	biospring::logging::info("Step: %5d", _nbiter);
	biospring::logging::info("Framerate: %5.2f", _framerate);
	biospring::logging::info("Kinetic energy: %5.2f kJ.mol-1", _energies.kinetic);
	if (isSpringEnabled())
		{
		biospring::logging::info("Spring energy: %5.2f kJ.mol-1", _energies.spring);
		biospring::logging::info("Dihedral energy: %5.2f kJ.mol-1", _energies.dihedral);
		}
	// Reported without the bond COUNT the CPU prints beside it: the count
	// lives in the device's slot arrays, and bringing them back every logged
	// step to print a number would cost a transfer the energy does not need.
	if (isHydrogenBondEnabled())
		biospring::logging::info("Hydrogen bond energy: %5.2f kJ.mol-1", _energies.hbond);
	// Measurements rather than energies, and the ones an IMPALA run is read on.
	if (isInsertionVectorEnabled() && _insertionVector)
		{
		biospring::logging::info("Insertion angle: %5.2lf °", _insertionVector->getAngle());
		biospring::logging::info("Roll angle: %5.2lf °", _insertionVector->getRollAngle());
		biospring::logging::info("Insertion depth: %5.2lf.", _insertionVector->getInsertionDepth());
		}
	}



void SpringNetworkOpenCL::convertSpringtoSpringocl(const Spring & spin, Springocl & spout)
	{
	spout.id1=spin.getParticle1().getId();
	spout.id2=spin.getParticle2().getId();
	spout.stiffness=spin.getStiffness();
	spout.equilibrium=spin.getEquilibrium();
	}
// Device -> host. The reverse of computeOpenCLPositions()/Velocities().
void SpringNetworkOpenCL::_syncParticlesFromDevice()
	{
	const unsigned n = SpringNetwork::getNumberOfParticles();
	for (unsigned i = 0; i < n; ++i)
		{
		Particle & particle = SpringNetwork::getParticle(i);
		const float x = _particlepositions[i].x;
		const float y = _particlepositions[i].y;
		const float z = _particlepositions[i].z;

		// The same check SpringNetwork::updateParticlePositions makes after
		// integrating, and for the same reason. It lives there, inside the CPU
		// integrator, which this path replaces wholesale: the kernel integrates
		// instead, so without this the two backends disagree on what a diverged
		// run does -- the CPU stops and says which particle went, the GPU
		// carries on silently. Measured on 074.DNADuplex at dt = 8 fs, the CPU
		// exited 1 with the message and the GPU exited 0.
		//
		// Only the dynamic particles, as on the CPU: a static one is never
		// integrated, so it cannot be sent non-finite by the integrator.
		//
		// It costs nothing worth counting -- measured at 0.98x to 1.05x on
		// systems from 1525 to 37200 dynamic particles, i.e. inside the noise.
		// isfinite on a float is a comparison against the exponent mask, the
		// value has just been read into a register, and the branch is never
		// taken.
		if (particle.isDynamic() && (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)))
			biospring::logging::die("Found non-finite position for particle %d.", particle.getId());

		particle.setPosition(Vector3f(x, y, z));
		particle.setVelocity(Vector3f(_particlevelocities[i].x, _particlevelocities[i].y,
		                              _particlevelocities[i].z));
		}
	}

void SpringNetworkOpenCL::computeParticleToSpringIndexes()
    {
    // A CSR offset array: N+1 entries, where entry i is where particle i's
    // springs begin and entry i+1 is where they end. A particle with no spring
    // has start == end and its loop simply does not run.
    //
    // This used to be N entries with -1 meaning "no spring", and the kernel
    // read the NEXT particle's entry to find where to stop. So a particle
    // whose successor had no springs got an end index of -1 and contributed
    // nothing, and the last particle fell back to the PARTICLE count used as
    // a SPRING index. Offsets remove both cases rather than guarding them.
    delete[] _particletospringindexes;
    _particletospringindexes = new int[_nbparticlesocl + 1];

    unsigned springIndex = 0;
    for (unsigned particleId = 0; particleId < _nbparticlesocl; ++particleId)
        {
        _particletospringindexes[particleId] = static_cast<int>(springIndex);
        while (springIndex < _nbspringsocl && _springsocl[springIndex].id1 == particleId)
            ++springIndex;
        }
    _particletospringindexes[_nbparticlesocl] = static_cast<int>(springIndex);
    }

float SpringNetworkOpenCL::distance (const float4 p1,const float4 p2)
		{
		float diffx=p1.x-p2.x;
		float diffy=p1.y-p2.y;
		float diffz=p1.z-p2.z;
		float diffw=p1.w-p2.w;

		return sqrt(diffx*diffx+diffy*diffy+diffz*diffz+diffw*diffw);
		}


void SpringNetworkOpenCL::computeOpenCLSprings()
    {

    const unsigned particleCount = SpringNetwork::getNumberOfParticles();
    const unsigned springCount = SpringNetwork::getNumberOfSprings();
    vector<vector<Springocl>> adjacentList(particleCount);

    delete[] _springparticlesindexes;
    _springparticlesindexes = springCount == 0 ? nullptr : new unsigned[springCount * 2];

    for (unsigned i = 0; i < springCount; ++i)
        {
        const Spring & spring = SpringNetwork::getSpring(i);
        const unsigned id1 = static_cast<unsigned>(spring.getParticle1().getId());
        const unsigned id2 = static_cast<unsigned>(spring.getParticle2().getId());
        if (id1 >= particleCount || id2 >= particleCount)
            throw std::runtime_error("OpenCL spring references an invalid particle id");

        _springparticlesindexes[i * 2] = id1;
        _springparticlesindexes[i * 2 + 1] = id2;

        Springocl forward{};
        convertSpringtoSpringocl(spring, forward);
        adjacentList[id1].push_back(forward);

        Springocl reverse = forward;
        std::swap(reverse.id1, reverse.id2);
        adjacentList[id2].push_back(reverse);
        }

    _nbspringsocl = springCount * 2;
    delete[] _springsocl;
    _springsocl = _nbspringsocl == 0 ? nullptr : new Springocl[_nbspringsocl];

    unsigned springIndex = 0;
    for (const auto & particleSprings : adjacentList)
        for (const Springocl & spring : particleSprings)
            _springsocl[springIndex++] = spring;

    }

void SpringNetworkOpenCL::computeOpenCLPositions()
		{

		_nbparticlesocl=SpringNetwork::getNumberOfParticles();
		delete[] _particlepositions;
		_particlepositions = _nbparticlesocl == 0 ? nullptr : new float4[_nbparticlesocl];


		Vector3f v;
		for(unsigned i=0;i<_nbparticlesocl;i++)
			{
			v=SpringNetwork::getParticle(i).getPosition();
			_particlepositions[i].x=v.getX();
			_particlepositions[i].y=v.getY();
			_particlepositions[i].z=v.getZ();
			_particlepositions[i].w=1.0f;
			}



		}

// The integration kernel divides the force by the mass, exactly as
// Particle::_integrateForce does. Before this buffer existed it did not, so
// the GPU integrated every particle as if it weighed 1 Da -- invisible on a
// toy system where that is true, wrong on any real structure.
void SpringNetworkOpenCL::computeOpenCLCharges()
	{
	delete[] _particlecharges;
	_particlecharges = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particlecharges[i] = SpringNetwork::getParticle(i).getCharge();
	}

void SpringNetworkOpenCL::computeOpenCLStericParameters()
	{
	delete[] _particleradii;
	delete[] _particleepsilons;
	_particleradii = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	_particleepsilons = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		{
		_particleradii[i] = SpringNetwork::getParticle(i).getRadius();
		_particleepsilons[i] = SpringNetwork::getParticle(i).getEpsilon();
		}
	}

// The .msp names the law; the kernel takes a number. Resolved from the same
// string SpringNetwork::_setupForceField switches on, so the two cannot pick
// different laws from one configuration.
int SpringNetworkOpenCL::_stericMode() const
	{
	const std::string mode = _config.steric.mode;
	if (mode == "lennard-jones-12-6Amber")
		return BIOSPRING_STERIC_AMBER_12_6;
	if (mode == "lennard-jones-8-6Lewitt")
		return BIOSPRING_STERIC_LEWITT_8_6;
	if (mode == "lennard-jones-8-6Zacharias")
		return BIOSPRING_STERIC_ZACHARIAS_8_6;
	return BIOSPRING_STERIC_LINEAR;
	}

void SpringNetworkOpenCL::computeOpenCLHydrophobicity()
	{
	delete[] _particlehydrophobicities;
	_particlehydrophobicities = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particlehydrophobicities[i] = SpringNetwork::getParticle(i).getHydrophobicity();
	}

// A .dx map, flattened row-major into one float4 per cell: the scalar in .x and
// the field PotentialGrid::compute_gradient already derived in .yzw. Nothing is
// recomputed here and nothing is interpolated -- the host differentiated once,
// when OpenDXReader read the file.
//
// Called once per map, because neither changes during a run. That is what makes
// these terms cheap on the device: one upload, then no transfer at all.
//
// The frame stores one INVERSE STEP PER AXIS rather than a cell width, because
// the maps are anisotropic: eleven of the twelve in the examples have a
// different step on each axis. boxmax carries GridCoordinatesSystem's -1e-6 so
// the kernel's bounds test is the one the CPU makes.
void SpringNetworkOpenCL::computeOpenCLMap(MapOnDevice & map,
                                           const biospring::grid::PotentialGrid & grid,
                                           const char * what)
	{
	delete[] map.cells;
	map.cells = nullptr;
	map.cellcount = 0;

	const std::array<size_t, 3> shape = grid.shape();
	const size_t total = shape[0] * shape[1] * shape[2];
	if (total == 0)
		return;

	const std::array<double, 3> step = grid.cell_size();
	const std::array<double, 3> origin = grid.origin();
	const biospring::Box & box = grid.boundaries();

	map.cells = new cl_float4[total];
	for (size_t i = 0; i < shape[0]; ++i)
		for (size_t j = 0; j < shape[1]; ++j)
			for (size_t k = 0; k < shape[2]; ++k)
				{
				// An explicit discrete_coordinates, not a braced list: at() is
				// overloaded on discrete_ and real_coordinates and a braced list
				// matches both.
				const biospring::grid::discrete_coordinates cell(
				    static_cast<int>(i), static_cast<int>(j), static_cast<int>(k));
				const biospring::grid::PotentialCell & c = grid.at(cell);
				cl_float4 & out = map.cells[(i * shape[1] + j) * shape[2] + k];
				out.s[0] = c.scalar;
				out.s[1] = c.vector.getX();
				out.s[2] = c.vector.getY();
				out.s[3] = c.vector.getZ();
				}
	map.cellcount = total;

	map.origin = {{static_cast<float>(origin[0]), static_cast<float>(origin[1]),
	               static_cast<float>(origin[2]), 0.0f}};
	map.invstep = {{static_cast<float>(1.0 / step[0]), static_cast<float>(1.0 / step[1]),
	                static_cast<float>(1.0 / step[2]), 0.0f}};
	map.boxmin = {{static_cast<float>(box.min_x()), static_cast<float>(box.min_y()),
	               static_cast<float>(box.min_z()), 0.0f}};
	map.boxmax = {{static_cast<float>(box.max_x() - 1e-6), static_cast<float>(box.max_y() - 1e-6),
	               static_cast<float>(box.max_z() - 1e-6), 0.0f}};
	map.shape = {{static_cast<cl_int>(shape[0]), static_cast<cl_int>(shape[1]),
	              static_cast<cl_int>(shape[2]), 0}};

	biospring::logging::info("OpenCL: %s map on the device, %d x %d x %d cells (%.1f MB), uploaded once",
	                         what, map.shape.s[0], map.shape.s[1], map.shape.s[2],
	                         total * sizeof(cl_float4) / 1048576.0);
	}

// The density map's per-particle weight. Unlike the charge it is not read from
// anything: ParticleProperty sets it to 1.0 and setBurying() has no caller, so
// today this array is all ones. Filled from the particles anyway rather than
// assumed, so that wiring setBurying() up does not silently leave the device
// behind.
void SpringNetworkOpenCL::computeOpenCLBuryings()
	{
	delete[] _particleburyings;
	_particleburyings = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particleburyings[i] = SpringNetwork::getParticle(i).getBurying();
	}

void SpringNetworkOpenCL::computeOpenCLSurfaces()
	{
	delete[] _particlesurfaces;
	delete[] _particletransfers;
	_particlesurfaces = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	_particletransfers = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		{
		const Particle & p = SpringNetwork::getParticle(i);
		_particlesurfaces[i] = p.getSolventAccessibilitySurface();
		_particletransfers[i] = p.getTransferEnergyByAccessibleSurface();
		}
	}

// FreeSASA's worker thread recomputes the surfaces every _freesasaState.step
// steps (1000 by default) and syncParticleStateData copies them into the
// Particle objects. There is no generation counter to ask, so this compares --
// one pass over N floats, against a transfer of the same N floats, so the
// comparison is worth making: on the default cadence it avoids 999 uploads out
// of 1000.
void SpringNetworkOpenCL::_refreshSurfacesIfChanged()
	{
	if (_particlesurfaces == nullptr)
		return;

	bool changed = false;
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		{
		const float surface = SpringNetwork::getParticle(i).getSolventAccessibilitySurface();
		if (surface != _particlesurfaces[i])
			{
			_particlesurfaces[i] = surface;
			changed = true;
			}
		}
	if (!changed)
		return;

	// An explicit write even though the buffer was created CL_MEM_USE_HOST_PTR:
	// on unified memory the device may well be reading the same pages already,
	// but OpenCL only guarantees that after a map or a write.
	_err = _queue.enqueueWriteBuffer(_inSurfaceBuffer, CL_TRUE, 0,
	                                 sizeof(float) * _nbparticlesocl, _particlesurfaces);
	checkErr("enqueueWriteBuffer(surface)");
	}

// The four membrane parameters default to 0.0 and no .msp key reaches them;
// only an MDDriver client can set them, through the "dmou", "dmol" and "dmtc"
// custom data. So this is cheap to re-read every step and has to be: a client
// can curve the membrane while the run is going.
bool SpringNetworkOpenCL::_membraneIsFlat() const
	{
	const biospring::forcefield::ForceField * ff = getForceField();
	return ff->getImpDoubleMembraneUpperMembOffset() == 0.0f &&
	       ff->getImpDoubleMembraneLowerMembOffset() == 0.0f &&
	       ff->getImpDoubleMembraneUpperMembTubeCurv() == 0.0f &&
	       ff->getImpDoubleMembraneLowerMembTubeCurv() == 0.0f;
	}

void SpringNetworkOpenCL::computeOpenCLMasses()
	{
	delete[] _particlemasses;
	_particlemasses = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particlemasses[i] = SpringNetwork::getParticle(i).getMass();
	}

// The CPU integrates only the particles on its dynamic list, and a static one
// is left entirely alone -- not moved, and not even force-reset, since
// resetForce() sits inside that same loop. The kernel had no notion of any of
// this: it walked every particle and used the mass as its only filter, which
// reproduces the CPU only because every static particle in today's examples is
// a massless ghost. pdb2spn --static freezes particles WITHOUT touching their
// mass, so a network built that way would have the GPU moving what the CPU
// holds.
void SpringNetworkOpenCL::computeOpenCLDynamicState()
	{
	delete[] _particledynamic;
	_particledynamic = _nbparticlesocl == 0 ? nullptr : new int[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particledynamic[i] = SpringNetwork::getParticle(i).isDynamic() ? 1 : 0;
	}

void SpringNetworkOpenCL::computeOpenCLForces()
{
	delete[] _particleforces;
	delete[] _hbondenergyper;
	delete[] _particleexternalforces;
	_particleforces = _nbparticlesocl == 0 ? nullptr : new float4[_nbparticlesocl];
	_particleexternalforces = _nbparticlesocl == 0 ? nullptr : new float4[_nbparticlesocl];
	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		_particleforces[i].x=0.0f;
		_particleforces[i].y=0.0f;
		_particleforces[i].z=0.0f;
		_particleforces[i].w=0.0f;
		_particleexternalforces[i].x=0.0f;
		_particleexternalforces[i].y=0.0f;
		_particleexternalforces[i].z=0.0f;
		_particleexternalforces[i].w=0.0f;

		}
}

void SpringNetworkOpenCL::computeOpenCLVelocities()
{
	delete[] _particlevelocities;
	_particlevelocities = _nbparticlesocl == 0 ? nullptr : new float4[_nbparticlesocl];

	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		_particlevelocities[i].x=0.0f;
		_particlevelocities[i].y=0.0f;
		_particlevelocities[i].z=0.0f;
		_particlevelocities[i].w=0.0f;
		}
}


void SpringNetworkOpenCL::computeRandomSet()
	{
	_nbparticlesocl=10000;
	_particlepositions=new float4[_nbparticlesocl];

	std::srand(static_cast<unsigned>(std::time(nullptr)));
	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		_particlepositions[i].x=((float) rand())/((float) RAND_MAX)-0.5;
		_particlepositions[i].y=((float) rand())/((float) RAND_MAX)-0.5;
		_particlepositions[i].z=((float) rand())/((float) RAND_MAX)-0.5;
		_particlepositions[i].w=0.0f;
		}
	_nbspringsocl=(_nbparticlesocl*(_nbparticlesocl-1));
	_springsocl=new Springocl[_nbspringsocl];
	unsigned springindex=0;
	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		for(unsigned j=0;j<_nbparticlesocl;j++)
			{
			if(i!=j)
				{
				_springsocl[springindex].id1=i;
				_springsocl[springindex].id2=j;
				_springsocl[springindex].equilibrium=distance(_particlepositions[i], _particlepositions[j]);
				_springsocl[springindex].stiffness=1.0;
				springindex++;
				}
			}
		}
	}

void SpringNetworkOpenCL::getParticlePosition(unsigned i, float position[3]) const
	{
	position[0]=_particlepositions[i].x;
	position[1]=_particlepositions[i].y;
	position[2]=_particlepositions[i].z;
	}

unsigned SpringNetworkOpenCL::getNumberOfParticles() const
	{
	return _nbparticlesocl;
	}

unsigned SpringNetworkOpenCL::getNumberOfSprings() const
	{
	return _nbspringsocl;
	}


// An interactor's force for the step about to run. This is MDDriver's pull:
// InteractorMDDriver::syncParticleStateData calls it once per particle per
// sync, and the CPU's version adds it to the particle's force accumulator,
// which the integrator then uses and resets.
//
// It used to write _particleforces, which is the host array _inoutForceBuffer
// was created over with CL_MEM_USE_HOST_PTR. On a unified-memory device that
// buffer IS this array, so the pull did reach the kernels -- by aliasing, not
// by any transfer, and therefore only where the runtime happens to alias.
// Nothing uploaded it, and on a discrete card it would have done nothing at
// all. _particleexternalforces and the `external` kernel exist for exactly
// this and had never been wired to anything.
void SpringNetworkOpenCL::setForce(unsigned i, float force[3])
	{
	if (_particleexternalforces == nullptr || i >= _nbparticlesocl)
		return;
	// A zero force is not a pull. biospring-cli registers an MDDriver
	// interactor whether or not anyone ever connects, and it calls this once
	// per particle per sync with whatever its force array holds -- zeros, until
	// somebody grabs something. Raising the flag for those meant uploading the
	// whole array and launching the kernel at every step of every run to add
	// nothing: 12.2 us a step on 013.GLIC, whose whole step is 420 us.
	if (force[0] == 0.0f && force[1] == 0.0f && force[2] == 0.0f)
		return;
	// Accumulated, like Particle::addForce, and cleared once the kernel has
	// read it -- the CPU's forces are reset at the end of every step too, so a
	// pull lasts exactly the step it was set for.
	_particleexternalforces[i].x += force[0];
	_particleexternalforces[i].y += force[1];
	_particleexternalforces[i].z += force[2];
	_externalforcespending = true;
	}

const char* SpringNetworkOpenCL::oclErrorString(cl_int error)
{
    static const char* errorString[] = {
        "CL_SUCCESS",
        "CL_DEVICE_NOT_FOUND",
        "CL_DEVICE_NOT_AVAILABLE",
        "CL_COMPILER_NOT_AVAILABLE",
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",
        "CL_OUT_OF_RESOURCES",
        "CL_OUT_OF_HOST_MEMORY",
        "CL_PROFILING_INFO_NOT_AVAILABLE",
        "CL_MEM_COPY_OVERLAP",
        "CL_IMAGE_FORMAT_MISMATCH",
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",
        "CL_BUILD_PROGRAM_FAILURE",
        "CL_MAP_FAILURE",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "CL_INVALID_VALUE",
        "CL_INVALID_DEVICE_TYPE",
        "CL_INVALID_PLATFORM",
        "CL_INVALID_DEVICE",
        "CL_INVALID_CONTEXT",
        "CL_INVALID_QUEUE_PROPERTIES",
        "CL_INVALID_COMMAND_QUEUE",
        "CL_INVALID_HOST_PTR",
        "CL_INVALID_MEM_OBJECT",
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
        "CL_INVALID_IMAGE_SIZE",
        "CL_INVALID_SAMPLER",
        "CL_INVALID_BINARY",
        "CL_INVALID_BUILD_OPTIONS",
        "CL_INVALID_PROGRAM",
        "CL_INVALID_PROGRAM_EXECUTABLE",
        "CL_INVALID_KERNEL_NAME",
        "CL_INVALID_KERNEL_DEFINITION",
        "CL_INVALID_KERNEL",
        "CL_INVALID_ARG_INDEX",
        "CL_INVALID_ARG_VALUE",
        "CL_INVALID_ARG_SIZE",
        "CL_INVALID_KERNEL_ARGS",
        "CL_INVALID_WORK_DIMENSION",
        "CL_INVALID_WORK_GROUP_SIZE",
        "CL_INVALID_WORK_ITEM_SIZE",
        "CL_INVALID_GLOBAL_OFFSET",
        "CL_INVALID_EVENT_WAIT_LIST",
        "CL_INVALID_EVENT",
        "CL_INVALID_OPERATION",
        "CL_INVALID_GL_OBJECT",
        "CL_INVALID_BUFFER_SIZE",
        "CL_INVALID_MIP_LEVEL",
        "CL_INVALID_GLOBAL_WORK_SIZE",
    };

    const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

    const int index = -error;

    return (index >= 0 && index < errorCount) ? errorString[index] : "";

}

#endif
