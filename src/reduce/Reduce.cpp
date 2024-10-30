
#include "Reduce.h"

namespace biospring
{
namespace reduce
{
namespace legacy
{
void Reduce::addReduceRuleName(const std::string & reducerulename)
{
    if (_reducerulenametoreduceruleid.count(reducerulename) == 0)
    {
        size_t id = _reducerulenametoreduceruleid.size();
        _reducerulenametoreduceruleid[reducerulename] = id;
        _reduceruleidtoreducerulename[id] = reducerulename;
    }
}

void Reduce::addAtomName(const std::string & atomname)
{
    if (_atomnametoatomid.count(atomname) == 0)
    {
        size_t id = _atomnametoatomid.size();
        _atomnametoatomid[atomname] = id;
        _atomidtoatomname[id] = atomname;
    }
}

void Reduce::addResidueName(const std::string & resname)
{
    if (_resnametoresid.count(resname) == 0)
    {
        size_t id = _resnametoresid.size();
        _resnametoresid[resname] = id;
        _residtoresname[id] = resname;
    }
}

int Reduce::getReduceRuleIdFromName(const std::string & reducerulename) const
{
    if (_reducerulenametoreduceruleid.count(reducerulename) == 0)
        return -1;
    else
        return (int)_reducerulenametoreduceruleid.at(reducerulename);
}

int Reduce::getResidueIdFromName(const std::string & resname) const
{
    if (_resnametoresid.count(resname) == 0)
        return -1;
    else
        return (int)_resnametoresid.at(resname);
}

int Reduce::getAtomIdFromName(const std::string & atomname) const
{

    if (_atomnametoatomid.count(atomname) == 0)
        return -1;
    else
        return (int)_atomnametoatomid.at(atomname);
}

const std::string Reduce::getReduceRuleNameFromId(size_t id) const
{
    if (_reduceruleidtoreducerulename.count(id) == 0)
        return "";
    else
        return _reduceruleidtoreducerulename.at(id);
}

const std::string Reduce::getResidueNameFromId(size_t id) const
{
    if (_residtoresname.count(id) == 0)
        return "";
    else
        return _residtoresname.at(id);
}

const std::string Reduce::getAtomNameFromId(size_t id) const
{

    if (_atomidtoatomname.count(id) == 0)
        return "";
    else
        return _atomidtoatomname.at(id);
}

void Reduce::addGroup(ParticleGroup * grp) { _groups.push_back(grp); }

std::vector<ReduceRule *> Reduce::getReduceRulesByResName(const std::string & resname)
{
    std::vector<ReduceRule *> rulelist;
    for (size_t i = 0; i < _reducerules.size(); i++)
    {
        ReduceRule * rule = &_reducerules[i];
        if (rule->getResidueName() == resname)
            rulelist.push_back(rule);
    }
    return rulelist;
}

void Reduce::print()
{
    for (size_t i = 0; i < _reducerules.size(); i++)
        _reducerules[i].print();

    for (size_t i = 0; i < _groups.size(); i++)
        _groups[i]->print();
}

void Reduce::checkOneParticlePerRule()
{
    _oneparticleperrule = true;
    for (size_t i = 0; i < _reducerules.size(); i++)
        if (_reducerules[i].isAllAtom() or _reducerules[i].getNumberOfAtoms() != 1)
        {
            _oneparticleperrule = false;
        }
}

bool Reduce::OneParticlePerRule() { return _oneparticleperrule; }

} // namespace legacy
} // namespace reduce
} // namespace biospring
