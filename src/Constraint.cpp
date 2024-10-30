#include "Constraint.h"
#include "Selection.h"
#include "measure.hpp"


Constraint::Constraint(Selection * sel1, Selection * sel2, float module)
{
    _selection1 = sel1;
    _selection2 = sel2;
    _forcemodule = module;
    _force = Vector3f();
    _distance = biospring::measure::distance(sel1->getBarycentre(), sel2->getBarycentre());
}

void Constraint::updateDistance(void)
{
    _distance = biospring::measure::distance(_selection1->getBarycentre(), _selection2->getBarycentre());
}

void Constraint::apply()
{
    _force = _selection2->getBarycentre() - _selection1->getBarycentre();
    _distance = _force.norm();
    _force.normalize();
    _force = _force * _forcemodule;
    _selection1->addForce(_force);
    _selection2->addForce(-_force);
}
