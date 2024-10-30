#ifdef OPENGL_SUPPORT

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "SpringNetworkViewer.h"
#include "Quaternion.h"

#ifdef OPENCL_SUPPORT
	#include "SpringNetworkOpenCL.h"
#endif

#include <GL/glew.h>
#if defined __APPLE__ || defined(MACOSX)
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif



int SpringNetworkViewer::window_width = 800;
int SpringNetworkViewer::window_height = 600;
int SpringNetworkViewer::glutWindowHandle = 0;
// mouse controls
int SpringNetworkViewer::mouse_old_x=0;
int SpringNetworkViewer::mouse_old_y=0;

int SpringNetworkViewer::mouse_lastclic_x=0;
int SpringNetworkViewer::mouse_lastclic_y=0;
bool SpringNetworkViewer::firsttime_clic=true;


int SpringNetworkViewer::mouse_buttons = 0;
SpringNetwork * SpringNetworkViewer::_springnetwork=NULL;

GLuint SpringNetworkViewer::particlespositionvbo=0;
GLuint SpringNetworkViewer::particlescolorvbo=0;
GLuint SpringNetworkViewer::particlesforcevbo=0;

GLfloat * SpringNetworkViewer::particlespos;
GLfloat * SpringNetworkViewer::particlescol;
GLfloat * SpringNetworkViewer::particlesforce;

unsigned SpringNetworkViewer::num=0;
GLfloat SpringNetworkViewer::camerainitposition[3];
GLfloat SpringNetworkViewer::barycentre[3];
GLfloat SpringNetworkViewer::camera[16];
GLfloat SpringNetworkViewer::radius=100.0f;
TrackBall SpringNetworkViewer::trackball=TrackBall(SpringNetworkViewer::window_width,SpringNetworkViewer::window_height,0.5);
int SpringNetworkViewer::mouseoverid=-1;
int SpringNetworkViewer::mouseoverlastid=-1;

Vector3f SpringNetworkViewer::force=Vector3f();
set <unsigned> SpringNetworkViewer::selection;
bool SpringNetworkViewer::interactionmode=false;
bool SpringNetworkViewer::navigationmode=false;




std::vector<Particle>::const_reference SpringNetworkViewer::getNearestParticleByCoords(float x, float y, float z, float radius) const
{
	const auto particles = _springnetwork->getParticles();
    float distance = 0.0;
    float distancemax = 0.0;
    Vector3f pos = Vector3f(x, y, z);
    unsigned target_id = 0;
    for (unsigned i = 0; i < particles.size(); i++)
    {
        distance = Vector3f::distance(particles[i].getPosition(), pos);
        if (distance > distancemax && distance <= radius)
        {
            distancemax = distance;
            target_id = i;
        }
    }
    return particles[target_id];
}



SpringNetworkViewer::SpringNetworkViewer()
	{
	}

SpringNetworkViewer::SpringNetworkViewer(SpringNetwork * springnetwork)
	{
	_springnetwork=springnetwork;
	}


void SpringNetworkViewer::init_gl(int argc, char** argv)
	{

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(SpringNetworkViewer::window_width, SpringNetworkViewer::window_height);
	glutInitWindowPosition (glutGet(GLUT_SCREEN_WIDTH)/2 - SpringNetworkViewer::window_width/2,
        glutGet(GLUT_SCREEN_HEIGHT)/2 - SpringNetworkViewer::window_height/2);

	std::stringstream ss;
	ss << "BioSpring OpenCL/OpenGL bindings"<< std::ends;
	SpringNetworkViewer::glutWindowHandle = glutCreateWindow(ss.str().c_str());

	glutDisplayFunc(SpringNetworkViewer::appRender); //main rendering function
	glutTimerFunc(10, SpringNetworkViewer::timerCB, 10); //determin a minimum time between frames
	glutKeyboardFunc(SpringNetworkViewer::appKeyboard);

	glutMouseFunc(SpringNetworkViewer::appMouse);
	glutMotionFunc(SpringNetworkViewer::appMotion);
	glutPassiveMotionFunc(SpringNetworkViewer::appPassiveMotion);

	glewInit();

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_DEPTH_TEST);
	// viewport
	glViewport(0, 0, SpringNetworkViewer::window_width, SpringNetworkViewer::window_height);

	// projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0, (GLfloat)SpringNetworkViewer::window_width / (GLfloat) SpringNetworkViewer::window_height, 0.1, 1000.0);

	num=_springnetwork->getNumberOfParticles();
	particlespos=new float[4*num];
	particlescol=new float[4*num];
	particlesforce=new float[4*num];



	for(unsigned i=0;i< num;i++)
		{
		particlescol[i*4+0]=1.0f;
		particlescol[i*4+1]=0.0f;
		particlescol[i*4+2]=0.0f;
		particlescol[i*4+3]=0.5f;


		float positions[3];
		_springnetwork->getParticlePosition(i,&(positions[0]));

		particlespos[i*4+0]=positions[0];
		particlespos[i*4+1]=positions[1];
		particlespos[i*4+2]=positions[2];
		particlespos[i*4+3]=1.0f;

		particlesforce[i*4+0]=0.0f;
		particlesforce[i*4+1]=0.0f;
		particlesforce[i*4+2]=0.0f;
		particlesforce[i*4+3]=0.0f;
		}

	auto centre = _springnetwork->getCentroid();
	barycentre[0]=centre.getX();
	barycentre[1]=centre.getY();
	barycentre[2]=centre.getZ();

	camerainitposition[0]=-barycentre[0];
	camerainitposition[1]=-barycentre[1];
	camerainitposition[2]=-barycentre[2]-radius;

	cout<<"Bary :"<<barycentre[0]<<" "<<barycentre[1]<<" "<<barycentre[2]<<endl;

	// set view matrix
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();


	glTranslatef(camerainitposition[0],camerainitposition[1], camerainitposition[2]);
	glGetFloatv(GL_MODELVIEW, camera);
	glPopMatrix();




	particlespositionvbo=createVBO( particlespos, num*sizeof(float)*4, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
	particlescolorvbo = createVBO( particlescol, num*sizeof(float)*4, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);
	particlesforcevbo = createVBO( particlesforce, num*sizeof(float)*4, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW);

	glFinish();
	cout<<"Particle vbo "<<particlespositionvbo<<" "<<particlescolorvbo<<" "<<particlesforcevbo<<endl;

	_springnetwork->initRun();
	}


void SpringNetworkViewer::mousePositionTo3DPoint(int x, int y, GLdouble &pX, GLdouble &pY, GLdouble &pZ)
	{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLfloat winX, winY, winZ;
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);

	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

	winX = (float)x;
	winY = (float)viewport[3]-(float)y;
	glReadPixels(x, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);


	gluUnProject(winX, winY, winZ, modelview, projection, viewport, &pX, &pY, &pZ);
	const Particle & p = getNearestParticleByCoords(pX,pY,pZ, 1.0);

	mouseoverlastid=mouseoverid;
	mouseoverid=p.getId();

	}



void SpringNetworkViewer::appRender()
	{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//this updates the particle system by calling the kernel

	_springnetwork->idleRun();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	Quaternion rotation=trackball.getRotation();
	Vector3f axis;
	float angle;
	rotation.toAxisAngle(&axis, &angle);

	glTranslatef(0.0f,0.0f, -radius);
	glRotatef(angle*180.0/3.1415,axis.getX(),axis.getY(),axis.getZ());
	glTranslatef(-barycentre[0],-barycentre[1], -barycentre[2]);

	glGetFloatv(GL_MODELVIEW, camera);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(8.0);


	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, particlespositionvbo);
	glVertexPointer(4, GL_FLOAT, 0, NULL);

	glEnableClientState(GL_COLOR_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER,particlescolorvbo);
	glColorPointer(4, GL_FLOAT, 0, NULL);

	updateColor();

	glDisableClientState(GL_NORMAL_ARRAY);

	glDrawArrays(GL_POINTS, 0, num);

#ifdef OPENCL_SUPPORT
	glDrawElements(GL_LINES, _springnetwork->getNumberOfSprings()*2,GL_UNSIGNED_INT,((SpringNetworkOpenCL *)_springnetwork)->_springparticlesindexes);
#endif

	glDisableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, particlesforcevbo);
	glVertexPointer(4, GL_FLOAT, 0, NULL);

	updateForce();

	//glDrawArrays(GL_LINES, 0, num);

	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glFlush();

	glutSwapBuffers();
	}


void SpringNetworkViewer::updateColor()
	{
	particlescol = (GLfloat*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_READ_WRITE_ARB);
	if(particlescol)
		{
		set <unsigned >::iterator it;
		for(it=selection.begin();it!=selection.end();++it)
			{
			unsigned id=*it;
			//cout<<" "<<id;
			particlescol[4*id]=0.0f;
			particlescol[4*id+1]=1.0f;
			particlescol[4*id+2]=0.0f;
			particlescol[4*id+3]=0.5f;
			}
		if(mouseoverid!=-1)
			{
			Particle * p=_springnetwork->getParticle(mouseoverid);
			Vector3f pos=p->getPosition();
			/*glPushMatrix();
			glTranslatef(p->getX(),p->getY(),p->getZ());
			glutSolidSphere(1.0,20,20);
			glPopMatrix();*/
			// wobble vertex in and out along normal
			particlescol[4*mouseoverid]=1.0f;
			particlescol[4*mouseoverid+1]=1.0f;
			particlescol[4*mouseoverid+2]=1.0f;
			particlescol[4*mouseoverid+3]=0.5f;
			}
		if(mouseoverlastid!=-1 && mouseoverid!=mouseoverlastid)
			{
			particlescol[4*mouseoverlastid]=1.0f;
			particlescol[4*mouseoverlastid+1]=0.0f;
			particlescol[4*mouseoverlastid+2]=0.0f;
			particlescol[4*mouseoverlastid+3]=0.5f;
			}
		glUnmapBufferARB(GL_ARRAY_BUFFER_ARB); // release pointer to mapping buffer
		}
	}
void SpringNetworkViewer::updateForce()
	{
	particlesforce = (GLfloat*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_READ_WRITE_ARB);
	if(particlesforce)
		{
		memset(particlesforce,0, sizeof(float)*num*4);
		set <unsigned >::iterator it;
		for(it=selection.begin();it!=selection.end();++it)
			{
			unsigned id=*it;
			//cout<<" "<<id;
			particlesforce[4*id]=force.getX();
			particlesforce[4*id+1]=force.getY();
			particlesforce[4*id+2]=force.getZ();
			particlesforce[4*id+3]=0.0f;
			}
		glUnmapBufferARB(GL_ARRAY_BUFFER_ARB); // release pointer to mapping buffer
		}
	}



void SpringNetworkViewer::appDestroy()
	{
	//this makes sure we properly cleanup our OpenCL context
	//delete example;
	deleteVBO(particlespositionvbo);
	deleteVBO(particlescolorvbo);

	if(glutWindowHandle)glutDestroyWindow(glutWindowHandle);
	printf("about to exit!\n");

	exit(0);
	}


void SpringNetworkViewer::timerCB(int ms)
	{
	//this makes sure the appRender function is called every ms miliseconds
	glutTimerFunc(ms, SpringNetworkViewer::timerCB, ms);
	glutPostRedisplay();
	}


void SpringNetworkViewer::appKeyboard(unsigned char key, int x, int y)
	{
	//this way we can exit the program cleanly
	switch(key)
		{
		case '\033': // escape quits
		case '\015': // Enter quits
		case 'Q':    // Q quits
		case 'q':    // q (or escape) quits
		// Cleanup up and quit
		//_springnetwork->endRun();
		SpringNetworkViewer::appDestroy();
		break;
		}
	}

void SpringNetworkViewer::appMouse(int button, int state, int x, int y)
	{
	if(glutGetModifiers()==GLUT_ACTIVE_SHIFT)
		{
		navigationmode=false;
		interactionmode=true;
		}
	else
		{
		navigationmode=true;
		interactionmode=false;
		}

	if (state == GLUT_DOWN)
		{
		if(mouseoverid!=-1)
			{
			updateSelection(mouseoverid);
			}
		else
			{
			SpringNetworkViewer::mouse_buttons |= 1<<button;
			trackball.mouseClic(x,y);
			}
		}
	else if (state == GLUT_UP)
		{
		SpringNetworkViewer::mouse_buttons = 0;
		}

	SpringNetworkViewer::mouse_old_x = x;
	SpringNetworkViewer::mouse_old_y = y;
	}

void SpringNetworkViewer::updateSelection(unsigned id)
	{
	set <unsigned>::iterator it;
	it=selection.find(id);
	if(it==selection.end())
		selection.insert(mouseoverid);
	else
		selection.erase(it);
	}



void SpringNetworkViewer::appMotion(int x, int y)
	{
	float dx, dy;

	if(firsttime_clic)
		{
		mouse_lastclic_x=x;
		mouse_lastclic_y=y;
		firsttime_clic=false;
		}
	dx = x - SpringNetworkViewer::mouse_old_x;
	dy = y - SpringNetworkViewer::mouse_old_y;

	if(navigationmode)
		{
		if (SpringNetworkViewer::mouse_buttons & 1)
			{
			trackball.mouseDrag(x,y);
			}
		else if (SpringNetworkViewer::mouse_buttons & 4)
			{
			SpringNetworkViewer::radius += dy * 0.1;
			}
		}
	if(navigationmode)
		{
		if (SpringNetworkViewer::mouse_buttons & 1)
			{
			trackball.mouseDrag(x,y);
			}
		else if (SpringNetworkViewer::mouse_buttons & 4)
			{
			SpringNetworkViewer::radius += dy * 0.1;
			}
		}
	if(interactionmode)
		{
		force=Vector3f(x-mouse_lastclic_x,y-mouse_lastclic_y,0.0);
		}
	SpringNetworkViewer::mouse_old_x = x;
	SpringNetworkViewer::mouse_old_y = y;

	}

void SpringNetworkViewer::appPassiveMotion(int x, int y)
	{
	float dx, dy;
	GLdouble px,py,pz;
	mousePositionTo3DPoint(x,y,px,py,pz);
	}

GLuint SpringNetworkViewer::getParticlesPositionVBO()
	{
	return particlespositionvbo;
	}

GLuint SpringNetworkViewer::getParticlesColorVBO()
	{
	return particlescolorvbo;
	}

GLuint SpringNetworkViewer::createVBO(const void* data, int dataSize, GLenum target, GLenum usage)
	{
	GLuint id = 0;  // 0 is reserved, glGenBuffersARB() will return non-zero id if success

	glGenBuffers(1, &id);                        // create a vbo
	glBindBuffer(target, id);                    // activate vbo id to use
	glBufferData(target, dataSize, data, usage); // upload data to video card

	// check data size in VBO is same as input array, if not return 0 and delete VBO
	//int bufferSize = 0;
	GLint&& bufferSize = 0;
	//glGetBufferParameteriv(target, GL_BUFFER_SIZE, &((GLint) bufferSize));
	glGetBufferParameteriv(target, GL_BUFFER_SIZE, &bufferSize);
	if(dataSize != bufferSize)
		{
		glDeleteBuffers(1, &id);
		id = 0;
		//cout << "[createVBO()] Data size is mismatch with input array\n";
		printf("[createVBO()] Data size is mismatch with input array\n");
		}
	//this was important for working inside blender!
	glBindBuffer(target, 0);
	return id;      // return VBO id
	}


void SpringNetworkViewer::deleteVBO(const GLuint vboId)
	{
	glDeleteBuffers(1, &vboId);
	}


#endif



