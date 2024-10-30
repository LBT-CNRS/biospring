#ifdef OPENCL_SUPPORT



#include "SpringNetworkOpenCL.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/CSVSampleWriter.h"


#include "Spring.h"
#include "Particle.h"
#include "Vector3f.h"




#if defined __APPLE__ || defined(MACOSX)
#else
    #if defined WIN32
    #else
        //needed for context sharing functions
        #include <GL/glx.h>
    #endif
#endif



#include <stdlib.h>
#include <math.h>


#define WORK_GROUP_SIZE 256


SpringNetworkOpenCL::SpringNetworkOpenCL() : SpringNetwork()
	{
	_nbspringsocl=0;
	_nbparticlesocl=0;
	_springparticlesindexes=NULL;
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


		_inSpringBuffer=cl::Buffer(
						  _context,
						  CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
						  sizeof(Springocl)*_nbspringsocl,
						  _springsocl,
						  &_err);
		checkErr( "Buffer::Buffer() 4");

		_inSpringIndexesBuffer=cl::Buffer(
								 _context,
								 CL_MEM_READ_ONLY| CL_MEM_USE_HOST_PTR,
								 sizeof(int)*_nbparticlesocl,
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

	ifstream _file("biospring.cl");
	_err=_file.is_open() ? CL_SUCCESS:-1;
	checkErr("biospring.cl");

	std::string prog(std::istreambuf_iterator<char>(_file),(std::istreambuf_iterator<char>()));

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


	printf("done building program\n");
	cerr << "Build Status: " << _program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(_devices[0]) << endl;
	cerr << "Build Options:\t" << _program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(_devices[0]) << endl;
	cerr << "Build Log:\t " << _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_devices[0]) << endl;

	checkErr("Program::build()");
	unsigned workgroupsize=WORK_GROUP_SIZE;
	unsigned globalsize=(_nbparticlesocl/workgroupsize)*(workgroupsize)+workgroupsize;

	cout<<globalsize<<endl;
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

	cerr<<__FUNCTION__<<" end"<<endl;
	}

void SpringNetworkOpenCL::InitOcl()
	{
	cerr<<__FUNCTION__<<" start"<<endl;


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
			#if defined WIN32 // Win32
			    _contextproperties=new cl_context_properties[7];
			    _contextproperties[0]=CL_GL_CONTEXT_KHR;
			    _contextproperties[1]=(cl_context_properties)wglGetCurrentContext();
			    _contextproperties[2]=CL_WGL_HDC_KHR;
			    _contextproperties[3]=(cl_context_properties)wglGetCurrentDC();
			    _contextproperties[4]=CL_CONTEXT_PLATFORM;
			    _contextproperties[5]=(cl_context_properties)(platforms[0])();
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
			    _contextproperties[5]=(cl_context_properties)(platforms[0])();
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
	#else`
		 _contextproperties=new cl_context_properties[3];
		_contextproperties[0] = CL_CONTEXT_PLATFORM;
		_contextproperties[1] = (cl_context_properties)(_platforms[0])();
		_contextproperties[2] = 0;

		cerr<<__FUNCTION__<<" test"<<endl;
		_context=cl::Context(CL_DEVICE_TYPE_GPU,_contextproperties,NULL,NULL,&_err);
		checkErr( "Context::Context()");
	#endif





	}

void SpringNetworkOpenCL::wrappingOcl()
	{
	cerr<<__FUNCTION__<<" start"<<endl;
	computeOpenCLPositions();
	//computeRandomSet();
	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		cout<<"x "<<_particlepositions[i].x<<" y "<<_particlepositions[i].y<<" z " <<_particlepositions[i].z<<std::endl;
		}

		for(unsigned i=0;i<_nbspringsocl;i++)
		{
		cout<<"id1 "<<_springsocl[i].id1<<" id2 "<<_springsocl[i].id2<<" stiffness " <<_springsocl[i].stiffness<<" equilibrium " <<_springsocl[i].equilibrium<<std::endl;
		}

		computeOpenCLVelocities();
		for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		cout<<"vx "<<_particlevelocities[i].x<<" vy "<<_particlevelocities[i].y<<" vz " <<_particlevelocities[i].z<<std::endl;
		}

		computeOpenCLForces();
		for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		cout<<"fx "<<_particleforces[i].x<<" fy "<<_particleforces[i].y<<" fz " <<_particleforces[i].z<<std::endl;
		}

		computeOpenCLSprings();

		computeParticleToSpringIndexes();

		for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		std::cout<<"Indexes "<<i<<" "<<_particletospringindexes[i]<<std::endl;
		}
	cerr<<__FUNCTION__<<" end"<<endl;
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

	_event=_kernelfunctorspring(_inoutPositionBuffer,_inSpringBuffer,_inSpringIndexesBuffer,_inoutForceBuffer, _nbparticlesocl, ForceField::GLOBALSPRINGFORCECONVERT);
	_event.wait();


	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	springtime=(endTime-startTime)*1.0E-9;


	_event=_kernelfunctordamping(_inoutForceBuffer,_inoutVelocityBuffer,_viscosity, _nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	dampingtime=(endTime-startTime)*1.0E-9;


	_event=_kernelfunctorexternal(_inoutForceBuffer,_inExternalForceBuffer,_nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	integrationtime=(endTime-startTime)*1.0E-9;

	_event=_kernelfunctorintegration(_inoutPositionBuffer,_inoutVelocityBuffer,_inoutForceBuffer, _timestep,_nbparticlesocl);
	_event.wait();
	startTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
	endTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
	externalforcetime=(endTime-startTime)*1.0E-9;

	_err = _queue.enqueueReadBuffer(_inoutVelocityBuffer,CL_TRUE,0,_nbparticlesocl,_particlevelocities);
	checkErr("ComamndQueue::enqueueReadBuffer()");

	_err = _queue.enqueueReadBuffer(_inoutForceBuffer,CL_TRUE,0,_nbparticlesocl,_particleforces);
	checkErr("ComamndQueue::enqueueReadBuffer()");

	_err = _queue.enqueueReadBuffer(_inoutPositionBuffer,CL_TRUE,0,_nbparticlesocl,_particlepositions);
	_queue.finish();
	checkErr("CommandQueue::enqueueReadBuffer()");


	_err = _queue.enqueueWriteBuffer(_inExternalForceBuffer,CL_TRUE,0,_nbparticlesocl,_particleexternalforces);
	_queue.finish();
	checkErr("CommandQueue::enqueueWriteBuffer()");




	SpringNetwork::idleRun();
	totaltime=springtime+dampingtime+integrationtime+externalforcetime;

	std::cout <<"Output particles :" << std::endl;
	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		cout<<"x "<<_particlepositions[i].x<<" y "<<_particlepositions[i].y<<" z " <<_particlepositions[i].z<<std::endl;
		cout<<"vx "<<_particlevelocities[i].x<<" vy "<<_particlevelocities[i].y<<" vz " <<_particlevelocities[i].z<<std::endl;
		cout<<"fx "<<_particleforces[i].x<<" fy "<<_particleforces[i].y<<" fz " <<_particleforces[i].z<<std::endl;
		}

	cl_ulong queueTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_QUEUED>();
	cl_ulong submitTime=_event.getProfilingInfo<CL_PROFILING_COMMAND_SUBMIT>();

	double runTime=endTime-startTime;
	std::cout<<"Times: "<<totaltime<<"( spring : "<<springtime<<", damping : "<<dampingtime<<", integration : "<<integrationtime<<", external : "<<externalforcetime<<")"<< std::endl;

	//usleep(10);

	#ifdef OPENGL_SUPPORT
		//Release the VBOs so OpenGL can play with them
		_err = _queue.enqueueReleaseGLObjects(&_allvbos, NULL, &_event);
		//printf("release gl: %s\n", oclErrorString(err));
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


void SpringNetworkOpenCL::run()
	{
	initRun();
	while(_nbiter<_maxiter || _maxiter<=0)
		{
		idleRun();

		}
	endRun();
	}



void SpringNetworkOpenCL::convertSpringtoSpringocl(const Spring & spin, Springocl & spout)
	{
	spout.id1=spin.getParticle1()->getId();
	spout.id2=spin.getParticle2()->getId();
	spout.stiffness=spin.getStiffness();
	spout.equilibrium=spin.getEquilibrium();
	}
void SpringNetworkOpenCL::computeParticleToSpringIndexes()
	{
	cerr<<__FUNCTION__<<" start"<<endl;
	_particletospringindexes=new int[_nbparticlesocl];
	bool cut=true;

	for(unsigned i=0;i<_nbparticlesocl;i++)
		{
		_particletospringindexes[i]=-1;
		}

	for(unsigned i=0;i<_nbspringsocl-1;i++)
		{
		if(cut)
			{
			_particletospringindexes[_springsocl[i].id1]=i;
			}
		cut=_springsocl[i].id1!=_springsocl[i+1].id1;
		}

	cerr<<__FUNCTION__<<" endl"<<endl;
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
		cerr<<__FUNCTION__<<" start"<<endl;
		vector < vector <Spring *> > adjacentlist;
		vector <Spring *> v;

		cerr<<__FUNCTION__<<" test1 "<<getNumberOfParticles()<<endl;
		for(unsigned i=0;i<getNumberOfParticles();i++)
			{
			adjacentlist.push_back(v);
			}
		cerr<<__FUNCTION__<<" test1 "<<SpringNetwork::getNumberOfSprings()<<endl;

		for(unsigned i=0;i<SpringNetwork::getNumberOfSprings();i++)
			{
			Spring * s1 = new Spring(*(_springs[i]));
			Spring * s2 = new Spring(s1->getParticle2(), s1->getParticle1(), s1->getStiffness(), s1->getEquilibrium());
			adjacentlist[s1->getParticle1()->getId()].push_back(s1);
			adjacentlist[s2->getParticle1()->getId()].push_back(s2);
			}
		_springparticlesindexes=new unsigned[SpringNetwork::getNumberOfSprings()*2];

		for(unsigned i=0;i<SpringNetwork::getNumberOfSprings();i++)
			{
			_springparticlesindexes[i*2]=_springs[i]->getParticle1()->getId();
			_springparticlesindexes[i*2+1]=_springs[i]->getParticle2()->getId();
			}



		cerr<<__FUNCTION__<<" test2"<<endl;
		_nbspringsocl=SpringNetwork::getNumberOfSprings()*2;
		_springsocl=new Springocl[_nbspringsocl];




		unsigned springindex=0;
		for(unsigned i=0;i<getNumberOfParticles();i++)
			{
			for(unsigned j=0;j<adjacentlist[i].size();j++)
				{
				convertSpringtoSpringocl(*(adjacentlist[i][j]),_springsocl[springindex]);
				springindex++;
				}
			}

		cerr<<__FUNCTION__<<" test"<<endl;
		for(unsigned i=0;i<getNumberOfParticles();i++)
			{
			for(unsigned j=0;j<adjacentlist[i].size();j++)
				{
				delete(adjacentlist[i][j]);
				}
			adjacentlist[i].clear();
			}
		adjacentlist.clear();
		cerr<<__FUNCTION__<<" end"<<endl;
		}

void SpringNetworkOpenCL::computeOpenCLPositions()
		{
		cerr<<__FUNCTION__<<" start"<<endl;

		_nbparticlesocl=SpringNetwork::getNumberOfParticles();
		_particlepositions = new float4[_nbparticlesocl];


		Vector3f v;
		for(unsigned i=0;i<_nbparticlesocl;i++)
			{
			v=_particules[i]->getPosition();
			_particlepositions[i].x=v.getX();
			_particlepositions[i].y=v.getY();
			_particlepositions[i].z=v.getZ();
			_particlepositions[i].w=1.0f;
			}
		cerr<<__FUNCTION__<<" end"<<endl;



		}

void SpringNetworkOpenCL::computeOpenCLForces()
{
	_particleforces = new float4[_nbparticlesocl];
	_particleexternalforces = new float4[_nbparticlesocl];
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
	_particlevelocities = new float4[_nbparticlesocl];

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

	srand(NULL);
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
	cerr<<"getParticlePosition i"<<i<<endl;
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
