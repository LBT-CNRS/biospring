
#include "rigidbodygroup/RigidBodyBuilder.h"

#include "logging.h"

namespace biospring
{
namespace rigidbodygroup
{

void RigidBodyBuilder::build()
{
    std::vector<ResidueParticleIndices> residues = _group_particles_by_residue();

    for (size_t index = 0; index < residues.size(); index++)
    {
        const std::string & resname = _topology.get_particle(residues[index][0]).properties().residue_name();
        RigidBodyRuleContainer rules_for_residue = _rules.get_rules_for_residue(resname);

        for (const RigidBodyRule & rule : rules_for_residue)
            _apply_rule(rule, residues, index);
    }

    logging::info("RigidBodyBuilder: created %zu spring(s) from rigid-body rules.", _topology.number_of_springs());
}

std::vector<RigidBodyBuilder::ResidueParticleIndices> RigidBodyBuilder::_group_particles_by_residue() const
{
    std::vector<ResidueParticleIndices> residues;
    if (_topology.number_of_particles() == 0)
        return residues;

    ResidueParticleIndices current = {0};

    for (size_t i = 1; i < _topology.number_of_particles(); i++)
    {
        const auto & previous_properties = _topology.get_particle(current.back()).properties();
        const auto & properties = _topology.get_particle(i).properties();

        bool same_residue = properties.residue_id() == previous_properties.residue_id() &&
                             properties.chain_name() == previous_properties.chain_name();

        if (same_residue)
            current.push_back(i);
        else
        {
            residues.push_back(current);
            current = {i};
        }
    }
    residues.push_back(current);

    return residues;
}

topology::Particle * RigidBodyBuilder::_resolve_atom(const std::string & atomname,
                                                      const std::vector<ResidueParticleIndices> & residues,
                                                      size_t index)
{
    std::string name = atomname;
    size_t target_index = index;

    if (!name.empty() && (name[0] == '+' || name[0] == '-'))
    {
        bool is_next = name[0] == '+';
        name = name.substr(1);

        if (is_next && index + 1 >= residues.size())
            return nullptr;
        if (!is_next && index == 0)
            return nullptr;

        target_index = is_next ? index + 1 : index - 1;

        // The neighbour must be in the same chain: a chain break behaves the
        // same as the actual first/last residue (no neighbour there either).
        const std::string & this_chain = _topology.get_particle(residues[index][0]).properties().chain_name();
        const std::string & other_chain = _topology.get_particle(residues[target_index][0]).properties().chain_name();
        if (this_chain != other_chain)
            return nullptr;
    }

    for (size_t particle_index : residues[target_index])
    {
        topology::Particle & p = _topology.get_particle(particle_index);
        if (p.properties().name() == name)
            return &p;
    }

    // Not found under its original name: if a reduction renamed the
    // particles (e.g. amber.grp: CA -> ACA for ALA), try the translated name.
    if (_translation != nullptr)
    {
        const std::string & target_resname = _topology.get_particle(residues[target_index][0]).properties().residue_name();
        std::string translated = _translate(target_resname, name);

        if (!translated.empty())
            for (size_t particle_index : residues[target_index])
            {
                topology::Particle & p = _topology.get_particle(particle_index);
                if (p.properties().name() == translated)
                    return &p;
            }
    }

    return nullptr;
}

std::string RigidBodyBuilder::_translate(const std::string & resname, const std::string & atomname) const
{
    reduce::ReduceRuleContainer rules = _translation->get_rules_for_residue(resname);

    for (const auto & rule : rules)
        if (rule.hasAtomNamed(atomname))
            return rule.getName();

    return "";
}

void RigidBodyBuilder::_check_translation_is_one_atom_per_rule() const
{
    for (const auto & rule : *_translation)
        if (rule.number_of_atoms() != 1)
            logging::die("RigidBodyBuilder: naming-translation file must have exactly one atom per rule (rule "
                         "'%s' for residue '%s' has %zu) -- it must be an all-atom identity mapping (like "
                         "amber.grp), not a coarse-grain reduction.",
                         rule.getName().c_str(), rule.getResidueName().c_str(), rule.number_of_atoms());
}

void RigidBodyBuilder::_apply_rule(const RigidBodyRule & rule, const std::vector<ResidueParticleIndices> & residues,
                                    size_t index)
{
    std::vector<topology::Particle *> resolved;

    for (const std::string & atomname : rule.getAtomNames())
    {
        topology::Particle * p = _resolve_atom(atomname, residues, index);
        if (p != nullptr)
            resolved.push_back(p);
    }

    for (size_t i = 0; i < resolved.size(); i++)
    {
        for (size_t j = i + 1; j < resolved.size(); j++)
        {
            try
            {
                _topology.add_spring(*resolved[i], *resolved[j], -1.0, _stiffness);
            }
            catch (const topology::SpringAlreadyExistsException &)
            {
            }
            catch (const topology::SelfSpringException &)
            {
            }
        }
    }
}

} // namespace rigidbodygroup
} // namespace biospring
