#ifdef OPENCL_SUPPORT



#include "SpringNetworkOpenCL.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/CSVSampleWriter.h"
#include "KernelSource.h"
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
    for (CellGrid * g : {&_stericcells, &_electrostaticcells, &_hydrophobiccells})
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
		cout<<"viewer"<<_viewer<<endl;
		cout<<"vbo"<<_viewer->getParticlesPositionVBO()<<endl;


		_inoutPositionBuffer=cl::BufferGL(_context, CL_MEM_READ_WRITE,_viewer->getParticlesPositionVBO() , &_err);
		_allvbos.push_back(_inoutPositionBuffer);
		checkErr( "Buffer::Buffer() 1");
		glFinish();

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

		// Only when a .dx was actually read: OpenCL rejects a zero-sized buffer,
		// and a run without electrostaticgrid.enable has no map at all.
		if (_electrostaticgridcellcount > 0)
			{
			_inElectrostaticGridBuffer=cl::Buffer(
									 _context,
									 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
									 sizeof(cl_float4)*_electrostaticgridcellcount,
									 _electrostaticgridcells,
									 &_err);
			checkErr( "Buffer::Buffer() electrostatic grid");
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
	_kernelelectrostatic = cl::Kernel(_program, "electrostatic", &_err);
	checkErr("Kernel::Kernel()");
	_kernelsteric = cl::Kernel(_program, "steric", &_err);
	checkErr("Kernel::Kernel()");
	_kernelhydrophobic = cl::Kernel(_program, "hydrophobic", &_err);
	checkErr("Kernel::Kernel()");
	_kernelelectrostaticfield = cl::Kernel(_program, "electrostaticfield", &_err);
	checkErr("Kernel::Kernel()");
	_kerneltorsion = cl::Kernel(_program, "torsion", &_err);
	checkErr("Kernel::Kernel()");

	}

void SpringNetworkOpenCL::InitOcl()
	{


	#ifdef OPENGL_SUPPORT
		#if defined (__APPLE__) || defined(MACOSX)
			CGLContextObj kCGLContext = CGLGetCurrentContext();
			CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);

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

		_context=cl::Context(CL_DEVICE_TYPE_GPU,_contextproperties,NULL,NULL,&_err);
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
	computeOpenCLElectrostaticGrid();
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
double electrostaticfieldtime=0.0;
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

	// The neighbour structure the non-bonded terms need. Built from the box the
	// positions read back last step occupy, which is one step stale -- hence the
	// cell of margin _buildCellList adds, the same reason the CPU's grid carries
	// a skin.
	_updateCellLists();

	#ifdef OPENGL_SUPPORT
		glFinish();
		// map OpenGL buffer object for writing from OpenCL
		//this passes in the vector of VBO buffer objects (position and color)
		_err = _queue.enqueueAcquireGLObjects(&_allvbos, NULL, &_event);
		//printf("acquire: %s\n", oclErrorString(err));
		_queue.finish();
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
                                 _nbparticlesocl, springForceScale);
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
    if (isStericEnabled() && _stericcells.ncellstotal > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const int springsenabled = (isSpringEnabled() && _nbspringsocl > 0) ? 1 : 0;

        unsigned a = 0;
        _kernelsteric.setArg(a++, _inoutPositionBuffer);
        _kernelsteric.setArg(a++, _inRadiusBuffer);
        _kernelsteric.setArg(a++, _inEpsilonBuffer);
        _kernelsteric.setArg(a++, _inoutForceBuffer);
        _kernelsteric.setArg(a++, _stericcells.headbuffer);
        _kernelsteric.setArg(a++, _stericcells.nextbuffer);
        _kernelsteric.setArg(a++, _stericcells.origin);
        _kernelsteric.setArg(a++, _stericcells.width);
        _kernelsteric.setArg(a++, _stericcells.ncells);
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
    if (isElectrostaticCoulombEnabled() && _electrostaticcells.ncellstotal > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;
        const int springsenabled = (isSpringEnabled() && _nbspringsocl > 0) ? 1 : 0;

        unsigned a = 0;
        _kernelelectrostatic.setArg(a++, _inoutPositionBuffer);
        _kernelelectrostatic.setArg(a++, _inChargeBuffer);
        _kernelelectrostatic.setArg(a++, _inoutForceBuffer);
        _kernelelectrostatic.setArg(a++, _electrostaticcells.headbuffer);
        _kernelelectrostatic.setArg(a++, _electrostaticcells.nextbuffer);
        _kernelelectrostatic.setArg(a++, _electrostaticcells.origin);
        _kernelelectrostatic.setArg(a++, _electrostaticcells.width);
        _kernelelectrostatic.setArg(a++, _electrostaticcells.ncells);
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

    // The precomputed potential map. No cell list and no neighbour walk: one
    // lookup per particle, in a buffer that was uploaded once and is never
    // touched again.
    if (isElectrostaticFieldEnabled() && _electrostaticgridcellcount > 0)
        {
        const unsigned wg = WORK_GROUP_SIZE;
        const unsigned global = (_nbparticlesocl / wg) * wg + wg;

        unsigned a = 0;
        _kernelelectrostaticfield.setArg(a++, _inoutPositionBuffer);
        _kernelelectrostaticfield.setArg(a++, _inChargeBuffer);
        _kernelelectrostaticfield.setArg(a++, _inoutForceBuffer);
        _kernelelectrostaticfield.setArg(a++, _inElectrostaticGridBuffer);
        _kernelelectrostaticfield.setArg(a++, _gridorigin);
        _kernelelectrostaticfield.setArg(a++, _gridinvstep);
        _kernelelectrostaticfield.setArg(a++, _gridshape);
        _kernelelectrostaticfield.setArg(a++, _gridboxmin);
        _kernelelectrostaticfield.setArg(a++, _gridboxmax);
        // getForceFieldScale(), which is what electrostaticgrid.scale sets --
        // NOT the coulomb scale, which belongs to the pairwise term. Particle::
        // addElectrostaticFieldForce reads the same one.
        _kernelelectrostaticfield.setArg(a++, getForceField()->getForceFieldScale());
        _kernelelectrostaticfield.setArg(a++, _nbparticlesocl);

        _err = _queue.enqueueNDRangeKernel(_kernelelectrostaticfield, cl::NullRange,
                                           cl::NDRange(global), cl::NDRange(wg), NULL, &_event);
        checkErr("enqueueNDRangeKernel(electrostaticfield)");
	_pendingevents.emplace_back(_event, &electrostaticfieldtime);
        }

    const float viscosity = isViscosityEnabled() ? getViscosity() : 0.0f;
    _event = _kernelfunctordamping(_inoutForceBuffer, _inoutVelocityBuffer,
                                  viscosity, _nbparticlesocl);
	_pendingevents.emplace_back(_event, &dampingtime);


	_event=_kernelfunctorexternal(_inoutForceBuffer,_inExternalForceBuffer,_nbparticlesocl);
	_pendingevents.emplace_back(_event, &integrationtime);

    _event = _kernelfunctorintegration(_inoutPositionBuffer, _inoutVelocityBuffer,
                                      _inoutForceBuffer, _inMassBuffer, _inDynamicBuffer,
                                      getTimeStep(), _nbparticlesocl);
	_pendingevents.emplace_back(_event, &externalforcetime);

	_err = _queue.enqueueReadBuffer(_inoutVelocityBuffer, CL_TRUE, 0,
        sizeof(float4) * _nbparticlesocl, _particlevelocities);
	checkErr("ComamndQueue::enqueueReadBuffer()");

	_err = _queue.enqueueReadBuffer(_inoutForceBuffer, CL_TRUE, 0,
        sizeof(float4) * _nbparticlesocl, _particleforces);
	checkErr("ComamndQueue::enqueueReadBuffer()");

	_err = _queue.enqueueReadBuffer(_inoutPositionBuffer, CL_TRUE, 0,
        sizeof(float4) * _nbparticlesocl, _particlepositions);
	_queue.finish();
	checkErr("CommandQueue::enqueueReadBuffer()");


	_err = _queue.enqueueWriteBuffer(_inExternalForceBuffer, CL_TRUE, 0,
        sizeof(float4) * _nbparticlesocl, _particleexternalforces);
	_queue.finish();
	checkErr("CommandQueue::enqueueWriteBuffer()");




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
		*pending.second += (e - s) * 1.0E-9;
		}
	_pendingevents.clear();

	_syncParticlesFromDevice();

	// SpringNetwork::idleRun() calls _resetEnergies(), so the energies have to
	// be filled after it, not before.
	SpringNetwork::idleRun();
	_computeEnergiesFromDeviceState();

	#ifdef OPENGL_SUPPORT
		//Release the VBOs so OpenGL can play with them
		_err = _queue.enqueueReleaseGLObjects(&_allvbos, NULL, &_event);
		_queue.finish();
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
	totaltime=springtime+torsiontime+sterictime+electrostatictime+electrostaticfieldtime+hydrophobictime+dampingtime+integrationtime+externalforcetime;
	std::cout<<"OpenCL kernel time: "<<totaltime<<" s ( spring: "<<springtime
	         <<", torsion: "<<torsiontime
	         <<", steric: "<<sterictime
	         <<", electrostatic: "<<electrostatictime
	         <<", electrostaticfield: "<<electrostaticfieldtime
	         <<", hydrophobic: "<<hydrophobictime
	         <<", damping: "<<dampingtime<<", integration: "<<integrationtime
	         <<", external: "<<externalforcetime<<" )"<<std::endl;
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


// Refreshes the grid of every enabled non-bonded term, each at its own cutoff.
//
// A term that is off gets no grid: the cheapest neighbour search is the one
// that does not run.
void SpringNetworkOpenCL::_updateCellLists()
	{
	if (isStericEnabled())
		_buildCellList(_stericcells, getStericCutoff());
	// Same distinction as at the dispatch: the cell list serves the PAIRWISE
	// kernel, and a grid-only configuration has no use for it.
	if (isElectrostaticCoulombEnabled())
		_buildCellList(_electrostaticcells, getElectrostaticCutoff());
	if (isHydrophobicityEnabled())
		_buildCellList(_hydrophobiccells, getHydrophobicCutoff());
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
bool SpringNetworkOpenCL::_measureCellGrid(CellGrid & grid, float cutoff)
	{
	if (_nbparticlesocl == 0 || cutoff <= 0.0f)
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

	// Margin on each side, in whole cells. It is what buys the frame its
	// lifetime: a structure has to expand by this much before anything leaves
	// and the pass above runs again. Two cells rather than one because a cell
	// is also the stencil's reach, so a particle in the outermost ring still
	// has its full neighbourhood inside the grid.
	const int MARGIN_CELLS = 2;
	const long CELL_LIMIT = 8L * 1024L * 1024L;

	// A cell WIDER than the cutoff is still correct: the 3x3x3 stencil then
	// covers more than the cutoff asks for, and the distance test each term
	// makes anyway drops the surplus. Only narrower would be wrong. So a box
	// too big for the cell count is answered by widening the cells rather than
	// by refusing to build -- which degrades towards brute force, slowly and
	// correctly, instead of leaving a term with no neighbour structure and no
	// way to know it.
	//
	// It is not reached by a healthy structure: the largest example here is
	// 034's capsid, a 292 A cube, which asks for 68921 cells of 8 A against the
	// 8M below. It is reached by a diverging one, and then saying so is worth
	// more than the grid.
	float width = cutoff;
	cl_int4 ncells;
	long total = 0;
	for (int attempt = 0; attempt < 64; ++attempt)
		{
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
			const int n = static_cast<int>(std::floor(span)) + 1 + 2 * MARGIN_CELLS;
			ncells.s[d] = n;
			total *= n;
			if (total > CELL_LIMIT)
				{
				overflowed = true;
				break;
				}
			}
		if (!overflowed)
			break;
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

	if (width != cutoff)
		BIOSPRING_WARN_ONCE("OpenCL cell list: the structure's box needs cells of %.2f A rather than the "
		                    "%.2f A cutoff asks for, to stay under %ld cells. Still exact, and slower: "
		                    "each walk sifts (%.1f)^3 times as many candidates.",
		                    width, cutoff, CELL_LIMIT, width / cutoff);

	ncells.s[3] = 0;
	for (int d = 0; d < 3; ++d)
		grid.origin.s[d] = lo[d] - MARGIN_CELLS * width;
	grid.origin.s[3] = 0.0f;
	grid.cutoff = cutoff;
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
	                                   cl::NDRange(cellglobal), cl::NDRange(wg));
	checkErr("enqueueNDRangeKernel(blankCells)");

	const unsigned partglobal = (_nbparticlesocl / wg) * wg + wg;
	_kernelbinparticles.setArg(0, _inoutPositionBuffer);
	_kernelbinparticles.setArg(1, grid.origin);
	_kernelbinparticles.setArg(2, grid.width);
	_kernelbinparticles.setArg(3, grid.ncells);
	_kernelbinparticles.setArg(4, grid.headbuffer);
	_kernelbinparticles.setArg(5, grid.nextbuffer);
	_kernelbinparticles.setArg(6, _nbparticlesocl);
	_err = _queue.enqueueNDRangeKernel(_kernelbinparticles, cl::NullRange,
	                                   cl::NDRange(partglobal), cl::NDRange(wg));
	checkErr("enqueueNDRangeKernel(binParticles)");
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
bool SpringNetworkOpenCL::_buildCellList(CellGrid & grid, float cutoff)
	{
	// A first call, a cutoff that changed under us (the .msp can be reloaded
	// mid-run), or a structure that has outgrown its frame.
	if (grid.ncellstotal == 0 || grid.cutoff != cutoff || !_frameStillHolds(grid))
		{
		// _measureCellGrid invalidates the grid when it cannot measure one, so
		// a term that consults it afterwards sees "no grid" rather than last
		// step's cells over this step's positions -- which is the one outcome
		// worse than having no grid at all, being wrong without saying so.
		if (!_measureCellGrid(grid, cutoff))
			return false;
		}

	_binParticlesIntoCells(grid);
	return true;
	}


// The 3x3x3 stencil around a particle's own cell, filtered to `cutoff`.
//
// A cell is a box of side `cutoff` and a cutoff is a sphere of radius `cutoff`,
// so the stencil is necessary and not sufficient: everything within the cutoff
// is in one of the 27 cells, and some of what is in them is beyond it. The
// distance test is what makes the answer exact, and a force kernel has to make
// the same one.
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
	for (int dz = -1; dz <= 1; ++dz)
		for (int dy = -1; dy <= 1; ++dy)
			for (int dx = -1; dx <= 1; ++dx)
				{
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
	if (isIMPEnabled())             add("impala");
	if (isInsertionVectorEnabled()) add("insertionvector");
	if (isRigidBodyEnabled())       add("rigidbody");

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

// The .dx potential map, flattened row-major into one float4 per cell:
// potential in .x and the field PotentialGrid::compute_gradient already derived
// in .yzw. Nothing is recomputed here and nothing is interpolated -- the host
// has done the differentiation once, at setup.
//
// Called once, because the map is read from the .dx in
// SpringNetwork::_setupElectrostatic and never changes afterwards. That is what
// makes this term cheap on the device: one upload, then no transfer at all for
// the rest of the run.
//
// The frame is stored as one INVERSE STEP PER AXIS rather than a cell width,
// because the maps are anisotropic: eleven of the twelve in the examples have a
// different step on each axis. _gridboxmax carries GridCoordinatesSystem's
// -1e-6 so the kernel's bounds test is the one the CPU makes.
void SpringNetworkOpenCL::computeOpenCLElectrostaticGrid()
	{
	delete[] _electrostaticgridcells;
	_electrostaticgridcells = nullptr;
	_electrostaticgridcellcount = 0;

	if (!isElectrostaticFieldEnabled())
		return;

	const biospring::grid::PotentialGrid & grid = getElectrostaticGrid();
	const std::array<size_t, 3> shape = grid.shape();
	const size_t total = shape[0] * shape[1] * shape[2];
	if (total == 0)
		return;

	const std::array<double, 3> step = grid.cell_size();
	const std::array<double, 3> origin = grid.origin();
	const biospring::Box & box = grid.boundaries();

	_electrostaticgridcells = new cl_float4[total];
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
				cl_float4 & out = _electrostaticgridcells[(i * shape[1] + j) * shape[2] + k];
				out.s[0] = c.scalar;
				out.s[1] = c.vector.getX();
				out.s[2] = c.vector.getY();
				out.s[3] = c.vector.getZ();
				}
	_electrostaticgridcellcount = total;

	_gridorigin = {{static_cast<float>(origin[0]), static_cast<float>(origin[1]),
	                static_cast<float>(origin[2]), 0.0f}};
	_gridinvstep = {{static_cast<float>(1.0 / step[0]), static_cast<float>(1.0 / step[1]),
	                 static_cast<float>(1.0 / step[2]), 0.0f}};
	_gridboxmin = {{static_cast<float>(box.min_x()), static_cast<float>(box.min_y()),
	                static_cast<float>(box.min_z()), 0.0f}};
	_gridboxmax = {{static_cast<float>(box.max_x() - 1e-6), static_cast<float>(box.max_y() - 1e-6),
	                static_cast<float>(box.max_z() - 1e-6), 0.0f}};
	_gridshape = {{static_cast<cl_int>(shape[0]), static_cast<cl_int>(shape[1]),
	               static_cast<cl_int>(shape[2]), 0}};

	biospring::logging::info("OpenCL: electrostatic map on the device, %d x %d x %d cells (%.1f MB), "
	                         "uploaded once",
	                         _gridshape.s[0], _gridshape.s[1], _gridshape.s[2],
	                         total * sizeof(cl_float4) / 1048576.0);
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


void SpringNetworkOpenCL::setForce(unsigned i, float force[3])
	{
	_particleforces[i].x=force[0];
	_particleforces[i].y=force[1];
	_particleforces[i].z=force[2];
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
