#ifndef __RIGIDBODYGROUP_RIGIDBODYRULE_H__
#define __RIGIDBODYGROUP_RIGIDBODYRULE_H__

#include <string>
#include <unordered_set>

namespace biospring
{
namespace rigidbodygroup
{

// Describes one named rigid-body group: a set of atom names, within a given
// residue type, that should all be pairwise spring-connected (a locally rigid
// substructure -- a bond/valence-angle tetrahedron, a peptide plane, a rigid
// ring...). Unrelated to biospring::rigidbody::RigidBody (classical rigid-body
// dynamics) -- this only describes which atoms belong together, nothing more.
//
// An atom name may carry a CHARMM-style "+"/"-" prefix, meaning "this atom in
// the next/previous residue of the chain" instead of the current one. The
// residue name always refers to the rule's own (unprefixed) residue.
class RigidBodyRule
{
  public:
    RigidBodyRule(const std::string & name, const std::string & residue_name)
        : _name(name), _residue_name(residue_name), _atomnames()
    {
    }

    RigidBodyRule() : RigidBodyRule("", "") {}

    void setName(const std::string & name) { _name = name; }
    const std::string & getName() const { return _name; }
    const std::string & name() const { return _name; }

    void setResidueName(const std::string & resname) { _residue_name = resname; }
    const std::string & getResidueName() const { return _residue_name; }
    const std::string & residue_name() const { return _residue_name; }

    const std::unordered_set<std::string> & getAtomNames() const { return _atomnames; }

    size_t getNumberOfAtoms() const { return _atomnames.size(); }
    size_t number_of_atoms() const { return _atomnames.size(); }

    void addAtom(const std::string & atomname) { _atomnames.insert(atomname); }

    bool hasAtomNamed(const std::string & name) const { return _atomnames.find(name) != _atomnames.end(); }

    void print() const;

  private:
    std::string _name;
    std::string _residue_name;
    std::unordered_set<std::string> _atomnames;
};

} // namespace rigidbodygroup
} // namespace biospring

#endif
