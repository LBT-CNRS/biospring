#include <iostream>

#include <GL/glew.h>
#if defined __APPLE__ || defined(MACOSX)
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif

#ifdef MDDRIVER_SUPPORT
	#include "interactor/mddriver/InteractorMDDriver.h"
#endif
#ifdef OPENCL_SUPPORT
	#include "SpringNetworkOpenCL.h"
#endif

#include "SpringNetwork.h"
#include "SpringNetworkViewer.h"


#include "IO/NetCDFReader.h"
#include "IO/MSParameterReader.h"
#define FILENAMEMAXSIZE 256
#define HOSTMAXSIZE 256
#define FILTERSIZE 256

using namespace std;


void usage()
        {
        cout<<"biospring"<<endl;
        cout<<"\tOptions : "<<endl;
        cout<<"\t\t -nc filename			(molecular system descriptions filename)"<<endl;
        cout<<"\t\t -msp filename 			(molecular simulation settings filename)"<<endl;
        cout<<""<<endl;
	#ifdef MDDRIVER_SUPPORT
	cout<<"\tMDDriver support     MDDriver support is enabled : "<<endl;
	cout<<"\t\tMDDriver is a library which extends the IMD protocol used in VMD for simulation and visualisation coupling."<<endl;
	cout<<"\t\t[-port portnumber] 			(Listening and outcoming port for MDDriver with default = 8888)"<<endl;
	cout<<"\t\t[-wait 0 or 1]     			(Wait a connection before starting the simulation with default = 1 means wait)"<<endl;
	cout<<"\t\t[-debug 0, 1 or 2] 	 		(Debug level of the MDDriver simulation with default = 0 means no debug message)"<<endl;
	cout<<"\t\t[-log filename]     			(Filename for MDDriver debug messages with default = stdout)"<<endl;
	//cout<<"\t\t[-opencl]     			(Use GPU/OpenCL for simulation)"<<endl;

	#endif

	cout<<"\t\t[-help] 			Print this message"<<endl;

        exit(1);
        }

int main(int argc, char ** argv)
	{
	bool openclok=true;
	#ifdef MDDRIVER_SUPPORT
   		unsigned port=8888;
		unsigned wait=1;
		unsigned debug=0;
        	char logfilename[FILENAMEMAXSIZE]="";
		float forcescale=1.0;
	#endif
       char ncname[FILENAMEMAXSIZE]="";
       char mspname[FILENAMEMAXSIZE]="";


	bool ncok=false;
//	bool dxok=false;
	bool mspok=false;

	#ifdef MDDRIVER_SUPPORT
		InteractorMDDriver * interactor=new InteractorMDDriver();
		interactor->setDebug(debug);
		interactor->setPort(port);
		interactor->setLog("");
		interactor->setWait(wait);
		interactor->setForceScale(forcescale);
	#endif

	for(unsigned i=1;i<(unsigned)argc;i++)
		{
		if(strcmp("-help",argv[i])==0)
			{
			usage();
			}
		else if(strcmp("-nc",argv[i])==0)
			{
			i++;
			if(i<(unsigned)argc)
				{
				strcpy(ncname,argv[i]);
				ncok=true;
				}
			else
				usage();
			}
		else if(strcmp("-msp",argv[i])==0)
			{
			i++;
			if(i<(unsigned)argc)
				{
				strcpy(mspname,argv[i]);
				mspok=true;
				}
			else
				usage();
			}
		#ifdef OPENCL_SUPPORT
			else if(strcmp("-opencl",argv[i])==0)
				{
				openclok=true;
				}

		#endif

		#ifdef MDDRIVER_SUPPORT

			else if(strcmp("-forcescale",argv[i])==0)
				{
				i++;
				if(i<(unsigned)argc)
					{
					forcescale=atof(argv[i]);
					interactor->setForceScale(forcescale);
					}
				else
					usage();
				}


			else if(strcmp("-port",argv[i])==0)
				{
				i++;
				if(i<(unsigned)argc)
					{
					port=atoi(argv[i]);
					interactor->setPort(port);
					}
				else
					usage();
				}

			else if(strcmp("-wait",argv[i])==0)
				{
				i++;
				if(i<(unsigned)argc)
					{
					wait=atoi(argv[i]);
					interactor->setWait(wait);
					}
				else
					usage();
				}
			else if(strcmp("-debug",argv[i])==0)
				{
				i++;
				if(i<(unsigned)argc)
					{
					debug=atoi(argv[i]);
					interactor->setDebug(debug);
					}
				else
					usage();
				}
			else if(strcmp("-log",argv[i])==0)
				{
				i++;
				if(i<(unsigned)argc)
					{
					strcpy(logfilename,argv[i]);
					interactor->setLog(logfilename);
					}
				else
					usage();
				}
			else
				{
				usage();
				}


		#endif
		}





	//cout<<"Creating grids (springcutoff : "<<springcutoff<<"  electrostaticcutoff : "<<electrostaticcutoff<<" vanderwaalscutoff : "<<vanderwaalscutoff<< ")..."<<endl;
	cout<<"Run simulation..."<<endl;

	#ifdef OPENCL_SUPPORT
		SpringNetwork * spn;
		if(openclok)
			spn= new SpringNetworkOpenCL();
		else
			spn=new SpringNetwork();
	#else
		SpringNetwork * spn=new SpringNetwork();
	#endif
	#if defined(MDDRIVER_SUPPORT))
		interactor->setSpringNetwork(spn);
		spn->setMDDriverInteractor(interactor);
	#endif




	NetCDFReader * netcdfreader=new NetCDFReader();
	netcdfreader->setSpringNetwork(spn);


	MSParameterReader * mspreader=new MSParameterReader();
	mspreader->setSpringNetwork(spn);

	if(ncok && mspok)
		{
		cout<<"Reading MSP file "<<mspname<<"..."<<endl;
		mspreader->setFileName(mspname);
		mspreader->read();
		cout<<"Reading NC file "<<ncname<<"..."<<endl;
		netcdfreader->setFileName(ncname);
		netcdfreader->read();
		}
	else
		{
		usage();
		}


	SpringNetworkViewer * spv=new SpringNetworkViewer(spn);
	spn->setSpringNetworkViewer(spv);
	spv->init_gl(argc, argv);






	glutMainLoop();

	delete netcdfreader;
	delete  mspreader;
	delete spn;
	return 0;
	}
