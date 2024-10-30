#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include "Vector3f.h"

class Selection;

class Constraint
	{
	public :
		Constraint(Selection * sel1,Selection * sel2, float module);
		void apply();

		void updateDistance(void);
		float getDistance() { return _distance; }

		Selection * getSelection1(void) { return _selection1; }
		Selection * getSelection2(void) { return _selection2; }

	private :

		Selection * _selection1;
		Selection * _selection2;
		Vector3f _force;
		float _forcemodule;
		float _distance;

	};

#endif
