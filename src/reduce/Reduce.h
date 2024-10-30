#ifndef _REDUCE_H_
#define _REDUCE_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "ParticleGroup.h"
#include "ReduceRule.h"

class ForceFieldData;

namespace biospring
{
namespace reduce
{
namespace legacy
{

class Reduce
{
  public:
    Reduce()
        : _reducerules(), _groups(), _reducerulenametoreduceruleid(), _atomnametoatomid(), _resnametoresid(),
          _reduceruleidtoreducerulename(), _atomidtoatomname(), _residtoresname(), _oneparticleperrule(false)
    {
    }

    void addReduceRule(const ReduceRule & rule) { _reducerules.push_back(rule); }

    size_t getNumberOfRules() const { return _reducerules.size(); }

    ReduceRule getReduceRule(const std::string & name) const { return _reducerules.at(getReduceRuleIdFromName(name)); }

    bool hasRuleNamed(const std::string & name) const
    {
        return _reducerulenametoreduceruleid.find(name) != _reducerulenametoreduceruleid.end();
    }

    void addReduceRuleName(const std::string & reducerulename);
    void addResidueName(const std::string & resname);
    void addAtomName(const std::string & atomname);

    int getReduceRuleIdFromName(const std::string & reducerulename) const;
    int getAtomIdFromName(const std::string & atomname) const;
    int getResidueIdFromName(const std::string & resname) const;
    const std::string getReduceRuleNameFromId(size_t id) const;
    const std::string getResidueNameFromId(size_t id) const;
    const std::string getAtomNameFromId(size_t id) const;

    std::vector<ReduceRule *> getReduceRulesByResName(const std::string & resname);

    void addGroup(ParticleGroup *);

    void reduce(ForceFieldData * ff);
    void checkOneParticlePerRule();
    bool OneParticlePerRule();
    void print();

  private:
    std::vector<ReduceRule> _reducerules;
    std::vector<ParticleGroup *> _groups;

    std::unordered_map<std::string, size_t> _reducerulenametoreduceruleid;
    std::unordered_map<std::string, size_t> _atomnametoatomid;
    std::unordered_map<std::string, size_t> _resnametoresid;

    std::unordered_map<size_t, std::string> _reduceruleidtoreducerulename;
    std::unordered_map<size_t, std::string> _atomidtoatomname;
    std::unordered_map<size_t, std::string> _residtoresname;
    bool _oneparticleperrule;
};

} // namespace legacy
} // namespace reduce
} // namespace biospring

#endif
