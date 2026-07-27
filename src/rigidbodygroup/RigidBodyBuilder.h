#ifndef __RIGIDBODYGROUP_RIGIDBODYBUILDER_H__
#define __RIGIDBODYGROUP_RIGIDBODYBUILDER_H__

#include <cstddef>
#include <string>
#include <vector>

#include "reduce/ReduceRule.h"
#include "reduce/ReduceRuleContainer.hpp"
#include "rigidbodygroup/RigidBodyRule.h"
#include "rigidbodygroup/RigidBodyRuleContainer.hpp"
#include "topology.hpp"

namespace biospring
{
namespace rigidbodygroup
{

// Builds springs directly from a set of RigidBodyRule groups: for every
// residue in the topology and every matching rule, resolves the rule's atom
// names (each possibly "+"/"-" prefixed for the next/previous residue) to
// actual particles and creates every pairwise spring between them.
//
// A referenced atom that does not exist (a missing chain neighbour at a
// terminus, proline's absent amide H, ...) is silently skipped rather than
// treated as an error. A spring that already exists (e.g. from the
// deliberate _PHI/_PSI overlap between consecutive residues) is silently
// left as-is rather than duplicated.
//
// Unrelated to biospring::rigidbody::RigidBody (classical rigid-body
// dynamics): this only ever adds springs to an existing all-atom topology.
class RigidBodyBuilder
{
  protected:
    using ResidueParticleIndices = std::vector<size_t>;

    topology::Topology & _topology;
    const RigidBodyRuleContainer & _rules;
    double _stiffness;

    // Optional naming-translation table: when a reduction (e.g. --grp
    // amber.grp --ff amber.ff) has renamed particles (CA -> ACA for Ala),
    // this lets atoms still be found by their original name. Must be an
    // all-atom identity mapping (exactly one atom per rule, like amber.grp)
    // -- never a real coarse-grain reduction; checked at construction time.
    const reduce::ReduceRuleContainer * _translation;

  public:
    RigidBodyBuilder(topology::Topology & topology, const RigidBodyRuleContainer & rules, double stiffness,
                      const reduce::ReduceRuleContainer * translation = nullptr)
        : _topology(topology), _rules(rules), _stiffness(stiffness), _translation(translation)
    {
        if (_translation != nullptr)
            _check_translation_is_one_atom_per_rule();
    }

    // Creates every rigid-body spring described by `_rules` on `_topology`.
    // Modifies the topology in place.
    void build();

  protected:
    // Groups the topology's particles into consecutive per-residue index
    // lists, in the order they appear in the topology (the same order a PDB
    // or NetCDF-derived topology naturally has: all atoms of a residue
    // together, residues in chain order).
    std::vector<ResidueParticleIndices> _group_particles_by_residue() const;

    // Applies one rule to the residue at `index`: resolves its atom names and
    // creates all pairwise springs among the ones that were found.
    void _apply_rule(const RigidBodyRule & rule, const std::vector<ResidueParticleIndices> & residues, size_t index);

    // Resolves a single (possibly +/- prefixed) atom name relative to the
    // residue at `index`. Returns nullptr if the neighbour residue does not
    // exist, belongs to a different chain, or does not contain that atom
    // (directly, or once translated through `_translation`).
    topology::Particle * _resolve_atom(const std::string & atomname, const std::vector<ResidueParticleIndices> & residues,
                                        size_t index);

    // Translates `atomname` for a residue named `resname` through
    // `_translation` (e.g. "CA" -> "ACA" for ALA). Returns an empty string if
    // no translation is found.
    std::string _translate(const std::string & resname, const std::string & atomname) const;

    // Dies if `_translation` contains a rule with more than one atom: this
    // table must be an all-atom identity mapping, never a real coarse-grain
    // reduction (which would make the translation ambiguous/meaningless).
    void _check_translation_is_one_atom_per_rule() const;
};

} // namespace rigidbodygroup
} // namespace biospring

#endif
