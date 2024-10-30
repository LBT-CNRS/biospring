#ifndef _TRACKBALL_H_
#define _TRACKBALL_H_

#include "Vector3f.h"
#include "Quaternion.h"


class TrackBall
	{
	public : 

		TrackBall(unsigned witdh, unsigned height, float radius);
				
		Quaternion getRotation() const;
		void mouseClic(unsigned x,unsigned y);
		void mouseDrag(unsigned x, unsigned y);
	private : 

		Quaternion _rotation;
		unsigned _width;
		unsigned _height;
		float _radius;
		Vector3f _newspherepoint;
		Vector3f _lastspherepoint;		
		
		Vector3f computeTrackBallPoint(float x, float y) const;
		
		
	};

#endif
