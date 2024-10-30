#include "TrackBall.h"
#include "Vector3f.h"
#include "Quaternion.h"

#include "math.h"


TrackBall::TrackBall(unsigned width, unsigned height, float radius)
	{
	_rotation=Quaternion();
	_width=width;
	_height=height;
	_radius=radius;
	}

Quaternion TrackBall::getRotation() const
	{
	return _rotation;
	}

void TrackBall::mouseClic(unsigned x,unsigned y)
	{
	if(x<=_width && y<=_height)
		{
		float rx = (((float) x)-_width/2.0)/(_height/2.0);
		float ry = -(((float) y)-_height/2.0)/(_height/2.0);
		_newspherepoint=computeTrackBallPoint(rx,ry);
		_lastspherepoint=_newspherepoint;
		}
	}
void TrackBall::mouseDrag(unsigned x, unsigned y)
	{
	if(x<=_width && y<=_height)
		{
		float rx = (((float) x)-_width/2.0)/(_height/2.0);
		float ry = -(((float) y)-_height/2.0)/(_height/2.0);
		_newspherepoint=computeTrackBallPoint(rx,ry);

		Quaternion q;
		q.rotationUsing2Vectors(_lastspherepoint, _newspherepoint);
		_rotation=q*_rotation;	

		_lastspherepoint=_newspherepoint;
		}
	}
	
Vector3f TrackBall::computeTrackBallPoint(float x, float y) const
	{
	if((x*x+y*y)<((_radius*_radius)/2.0f))
		{
		return Vector3f(x,y,sqrt(_radius*_radius-x*x-y*y));
		}
	else
		{
		return Vector3f(x,y,(_radius*_radius)/(2.0f*sqrt(x*x+y*y)));
		}
	}
