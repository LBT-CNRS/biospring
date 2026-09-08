#ifdef OPENCL_SUPPORT



#include "SpringNetworkOpenCL.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/CSVSampleWriter.h"
#include "KernelSource.h"

#include <fstream>


#include "Spring.h"
#include "Particle.h"
#include "Vector3f.h"
#include "forcefield/constants.hpp"

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
      _particleexternalforces(nullptr), _particlemasses(nullptr), _particletospringindexes(nullptr), _springsocl(nullptr),
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
    delete[] _springsocl;
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
	computeOpenCLSprings();
	computeParticleToSpringIndexes();
	}

double springtime=0.0;
double dampingtime=0.0;
double integrationtime=0.0;
double externalforcetime=0.0;
double totaltime=0.0;



void SpringNetworkOpenCL::idleRun()
	{
	cl_ulong startTime;
	cl_ulong endTime;
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
                                 _nbparticlesocl, springForceScale);
	_event.wait();


	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	springtime+=(endTime-startTime)*1.0E-9;
    }


    const float viscosity = isViscosityEnabled() ? getViscosity() : 0.0f;
    _event = _kernelfunctordamping(_inoutForceBuffer, _inoutVelocityBuffer,
                                  viscosity, _nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	dampingtime+=(endTime-startTime)*1.0E-9;


	_event=_kernelfunctorexternal(_inoutForceBuffer,_inExternalForceBuffer,_nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	integrationtime+=(endTime-startTime)*1.0E-9;

    _event = _kernelfunctorintegration(_inoutPositionBuffer, _inoutVelocityBuffer,
                                      _inoutForceBuffer, _inMassBuffer, getTimeStep(),
                                      _nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	externalforcetime+=(endTime-startTime)*1.0E-9;

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
	_syncParticlesFromDevice();

	SpringNetwork::idleRun();

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
	totaltime=springtime+dampingtime+integrationtime+externalforcetime;
	std::cout<<"OpenCL kernel time: "<<totaltime<<" s ( spring: "<<springtime
	         <<", damping: "<<dampingtime<<", integration: "<<integrationtime
	         <<", external: "<<externalforcetime<<" )"<<std::endl;
	SpringNetwork::endRun();
	}


void SpringNetworkOpenCL::run()
	{
	initRun();
	while (!isEnd())
		{
		idleRun();

		}
	endRun();
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
		particle.setPosition(Vector3f(_particlepositions[i].x, _particlepositions[i].y,
		                              _particlepositions[i].z));
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
void SpringNetworkOpenCL::computeOpenCLMasses()
	{
	delete[] _particlemasses;
	_particlemasses = _nbparticlesocl == 0 ? nullptr : new float[_nbparticlesocl];
	for (unsigned i = 0; i < _nbparticlesocl; i++)
		_particlemasses[i] = SpringNetwork::getParticle(i).getMass();
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
