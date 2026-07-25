#ifndef _SPRINGNTWORKVIEWER_H_
#define _SPRINGNTWORKVIEWER_H_

#ifdef OPENGL_SUPPORT

#if defined(_WIN32)
    #include <windows.h>
#endif
#include <GL/glew.h>

#if defined(__APPLE__) || defined(MACOSX)
    #if __has_include(<OpenGL/glu.h>)
        #include <OpenGL/glu.h>
    #else
        #include <GL/glu.h>
    #endif
    #if __has_include(<GLUT/glut.h>)
        #include <GLUT/glut.h>
    #elif __has_include(<GL/freeglut.h>)
        #include <GL/freeglut.h>
    #else
        #include <GL/glut.h>
    #endif
#else
    #include <GL/glu.h>
    #if __has_include(<GL/freeglut.h>)
        #include <GL/freeglut.h>
    #else
        #include <GL/glut.h>
    #endif
#endif

/*#ifdef OPENCL_SUPPORT
	#include "cl.hpp"
#endif*/

#include <set>
#include "SpringNetwork.h"
#include "TrackBall.h"

using biospring::spn::Particle;
using biospring::spn::SpringNetwork;

using namespace std;

class SpringNetworkViewer
	{
	public :
		static int window_width ;
		static int window_height ;
		static int glutWindowHandle ;

		static int mouse_buttons ;

		static int mouse_old_x;
		static int mouse_old_y;

		static int mouse_lastclic_x;
		static int mouse_lastclic_y;
		static bool firsttime_clic;

		static GLuint particlespositionvbo;
		static GLuint particlescolorvbo;
		static GLuint particlesforcevbo;

		static GLfloat * particlespos;
		static GLfloat * particlescol;
		static GLfloat * particlesforce;

		static unsigned num;
		static GLfloat camerainitposition[3];
		static GLfloat camera[16];
		static GLfloat barycentre[3];
		static GLfloat radius;
		static TrackBall trackball;
		static int mouseoverid;
		static int mouseoverlastid;
		static set<unsigned> selection;



		static bool interactionmode;
		static bool navigationmode;

		static Vector3f force;

		static GLuint getParticlesPositionVBO();
		static GLuint getParticlesColorVBO();



		static void init_gl(int argc, char** argv);
		static void appRender();
		static void appDestroy();
		static void timerCB(int ms);
		static void appKeyboard(unsigned char key, int x, int y);

		static void appMouse(int button, int state, int x, int y);
		static void appMotion(int x, int y);
		static void appPassiveMotion(int x, int y);

		static void updateSelection(unsigned id);
		static void updateColor();
		static void updateForce();



		SpringNetworkViewer();
		SpringNetworkViewer(SpringNetwork * springnetwork);
		static GLuint createVBO(const void* data, int dataSize, GLenum target, GLenum usage);
		static void deleteVBO(const GLuint vboId);
		static void mousePositionTo3DPoint(int x, int y, GLdouble &pX, GLdouble &pY, GLdouble &pZ);
	private :
		static SpringNetwork * _springnetwork;
		static Vector3f computeTrackBallPoint(float x,float y,float radius);

		// Returns the closest particle given some coordinates.
        static std::vector<Particle>::const_reference getNearestParticleByCoords(float x, float y, float z, float radius);

	};

#endif
#endif
