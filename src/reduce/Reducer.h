
//
// Defines the reducer class which reduces an all atom topology
// to a coarse grain topology.
//

#ifndef __REDUCER_H__
#define __REDUCER_H__

#include <string>
#include <vector>

#include "Particle.h"
#include "Reduce.h"
#include "ReduceRuleContainer.hpp"
#include "SpringNetwork.h"
#include "forcefield/ForceField.h"
#include "logging.h"
#include "measure.hpp"
#include "topology.hpp"

namespace biospring
{
namespace reduce
{

using ParticleContainer = std::vector<topology::Particle>;

//
// Storage of reduction to coarse grain parameters.
//
// This is used by command-line parsers to store parameters and pass it to
// actual Reducer instance.
//
struct ReductionParameters
{
    std::string pathGroup;
    std::string pathForceField;
    bool ignoreMissing;   // ignore missing particles?
    bool ignoreDuplicate; // ignore duplicate particles?

    ReductionParameters() : pathGroup(""), pathForceField(""), ignoreMissing(false), ignoreDuplicate(false) {}
};

// Used to build a specific grain from a set of particles and a rule.
class GrainBuilder
{
  protected:
    using Grain = topology::Particle;

    const ParticleContainer & _input_particles;
    const ReduceRule & _rule;
    const forcefield::ForceField & _forcefield;
    bool _ignore_duplicate_particles;
    bool _ignore_missing_particles;
    Grain _grain; // output grain
    bool _grain_built;

  public:
    GrainBuilder(const ParticleContainer & input_particles, const ReduceRule & rule,
                 const forcefield::ForceField & forcefield)
        : _input_particles(input_particles), _rule(rule), _forcefield(forcefield), _ignore_duplicate_particles(false),
          _ignore_missing_particles(false), _grain(), _grain_built(false)
    {
    }

    // =============================================================================
    // Accessors.
    // =============================================================================
    void set_ignore_duplicate_particles(bool value) { _ignore_duplicate_particles = value; }
    bool get_ignore_duplicate_particles(void) const { return _ignore_duplicate_particles; }

    void set_ignore_missing_particles(bool value) { _ignore_missing_particles = value; }
    bool get_ignore_missing_particles(void) const { return _ignore_missing_particles; }

    const Grain & grain(void)
    {
        if (!_grain_built)
            throw std::runtime_error("Grain not built yet!");
        return _grain;
    }

    // =============================================================================
    // Main methods.
    // =============================================================================

    // Returns a grain built from the input particles and the rule.
    bool build(void)
    {
        ParticleContainer grain_particles;

        // Adds particles from `_input_particles` that are defined in `_rule`.
        for (const topology::Particle & p : _input_particles)
            if (_rule.hasAtomNamed(p.properties().name()))
                grain_particles.push_back(p);

        // Sanity checks: grain is not empty, has enough particles, and has not too many particles.
        bool success = _check_not_enough_particles(grain_particles) && _check_too_many_particles(grain_particles);
        if (!success)
            return false;

        // Creates the grain.
        _grain = _create_grain(grain_particles);
        _grain_built = true;

        return true;
    }

  protected:
    std::string _grain_name() const { return _rule.name(); }
    std::string _residue_name() const { return _rule.residue_name(); }
    int _residue_id() const { return _input_particles[0].properties().residue_id(); }

    // Returns true if a grain contains a particle with a given name.
    static bool contains_particle(const ParticleContainer & particles, const std::string & name)
    {
        for (const auto & p : particles)
            if (p.properties().name() == name)
                return true;
        return false;
    }

    // Creates a grain from a set of particles.
    Grain _create_grain(const ParticleContainer & grain_particles) const
    {
        Grain grain;
        grain.properties().set_name(_grain_name());
        grain.properties().set_residue_name(_residue_name());
        grain.properties().set_residue_id(_residue_id());
        grain.properties().set_chain_name(_input_particles[0].properties().chain_name());
        grain.properties().set_position(measure::centroid(grain_particles));

        // Set properties from forcefield.
        if (_forcefield.hasProperty(_grain_name()))
        {
            const spn::ParticleProperty & pp = _forcefield.getPropertiesFromName(_grain_name());
            grain.properties().set_charge(pp.getCharge());
            grain.properties().set_radius(pp.getRadius());
            grain.properties().set_mass(pp.getMass());
            grain.properties().set_epsilon(pp.getEpsilon());

            grain.properties().set_hydrophobicity(pp.getHydrophobicity());
            // grain.properties().set_imp(topology::IMPProperties::build().transfert_energy_by_accessible_surface(pp.getTransferEnergyByAccessibleSurface()));
            topology::IMPProperties impProperties = topology::IMPProperties::build();
            impProperties.set_solvent_accessible_surface(pp.getSolventAccessibilitySurface());
            impProperties.set_transfert_energy_by_accessible_surface(pp.getTransferEnergyByAccessibleSurface());
            grain.properties().set_imp(impProperties);
            // grain.properties().imp().set_transfert_energy_by_accessible_surface(
            //     pp.getTransferEnergyByAccessibleSurface());
        }
        else
        {
            logging::warning("No Forcefield found for '%s'!", _grain_name().c_str());
        }

        return grain;
    }

    // Checks if a grain is empty.
    // If so, displays a warning message.
    // Returns the success status, which is true if grain is not empty, false otherwise.
    bool _check_not_empty(const ParticleContainer & grain) const
    {
        bool success = true;
        if (grain.empty())
        {
            logging::warning("Residue %s:%d: No particle found to create grain %s...skipping it",
                             _residue_name().c_str(), _residue_id(), _grain_name().c_str());
            success = false;
        }
        return success;
    }

    // Checks if a grain has not enough particles.
    // If so, displays a warning message.
    // Returns the success status, which is true if grain has enough particles, or
    // false if grain has not enough particles and `_ignore_missing_particles` is false.
    bool _check_not_enough_particles(const ParticleContainer & grain) const
    {
        if (!_check_not_empty(grain))
            return false;

        bool success = true;
        if (grain.size() < _rule.number_of_atoms())
        {
            // Grain has missing particles.
            // Finds which one are missing to display warning message.
            for (const std::string & name : _rule.getAtomNames())
            {
                if (!contains_particle(grain, name))
                {
                    if (_ignore_missing_particles)
                    {
                        logging::warning("Residue %s:%d: Particle '%s' required for grain '%s' not found in "
                                         "topology...still adding grain to model",
                                         _residue_name().c_str(), _residue_id(), name.c_str(), _grain_name().c_str());
                    }
                    else
                    {
                        logging::warning("Residue %s:%d: Particle '%s' required for grain '%s' not found in "
                                         "topology...skipping this grain",
                                         _residue_name().c_str(), _residue_id(), name.c_str(), _grain_name().c_str());
                        success = false;
                    }
                }
            }
        }

        return success;
    }

    // Checks if a grain has too many particles.
    // If so, displays a warning message.
    // Returns the success status, which is true if grain has not too many particles, or
    // false if grain has too many particles and `_ignore_duplicate_particles` is false.
    bool _check_too_many_particles(const ParticleContainer & grain) const
    {
        bool success = true;
        if (grain.size() > _rule.number_of_atoms())
        {
            if (_ignore_duplicate_particles)
            {
                logging::warning("Residue %s:%d: Duplicate particles in grain %s...still adding grain to model",
                                 _residue_name().c_str(), _residue_id(), _grain_name().c_str());
            }
            else
            {
                logging::warning("Residue %s:%d: Duplicate particles in grain %s...skipping this grain",
                                 _residue_name().c_str(), _residue_id(), _grain_name().c_str());
                success = false;
            }
        }
        return success;
    }
};

class ResidueReducer
{
  protected:
    const ParticleContainer & _particles;
    const ReduceRuleContainer & _rules;
    const forcefield::ForceField & _forcefield;
    bool _ignore_missing_particles;
    bool _ignore_duplicate_particles;

  public:
    ResidueReducer(const ParticleContainer & particles, const ReduceRuleContainer & rules,
                   const forcefield::ForceField & forcefield)
        : _particles(particles), _rules(rules), _forcefield(forcefield), _ignore_missing_particles(false),
          _ignore_duplicate_particles(false)
    {
    }

    // =============================================================================
    // Accessors.
    // =============================================================================
    void set_ignore_missing_particles(bool value) { _ignore_missing_particles = value; }
    bool get_ignore_missing_particles(void) const { return _ignore_missing_particles; }

    void set_ignore_duplicate_particles(bool value) { _ignore_duplicate_particles = value; }
    bool get_ignore_duplicate_particles(void) const { return _ignore_duplicate_particles; }

    // =============================================================================
    // Main methods.
    // =============================================================================

    // Returns a vector of grains representing the residue.
    std::vector<topology::Particle> build(void) const
    {
        // Displays warnings if some particles are present in the residue but not in the rule.
        _check_unknown_particles();

        std::vector<topology::Particle> grains;

        for (const auto & rule : _rules)
        {
            GrainBuilder grain_builder(_particles, rule, _forcefield);
            grain_builder.set_ignore_duplicate_particles(_ignore_duplicate_particles);
            grain_builder.set_ignore_missing_particles(_ignore_missing_particles);

            bool success = grain_builder.build();
            if (success)
                grains.push_back(grain_builder.grain());
        }

        return grains;
    }

  protected:
    // Displays a warning message if some particles are present in the residue but not in the rule.
    void _check_unknown_particles() const
    {
        for (const topology::Particle & p : _particles)
            if (!_particle_belongs_to_a_grain(p))
                logging::warning("Particle '%s' of residue '%s':'%d' not found in any grain", p.properties().name().c_str(),
                                                                                              p.properties().residue_name().c_str(),
                                                                                              p.properties().residue_id());
    }

    // Returns true if a particle belongs to a grain, i.e. if particle's names is found in at least one rule.
    bool _particle_belongs_to_a_grain(const topology::Particle & p) const
    {
        for (const ReduceRule & rule : _rules)
            if (rule.hasAtomNamed(p.properties().name()))
                return true;

        return false;
    }
};

class Reducer
{
  protected:
    using ParticleGroups = std::vector<ParticleContainer>;

    const topology::Topology & _source_topology;

    forcefield::ForceField _forcefield;
    bool _forcefield_initialized;

    ReduceRuleContainer _rules;
    bool _rules_initialized;

    topology::Topology _target_topology;
    bool _ignore_duplicate_particles;
    bool _ignore_missing_particles;

  public:
    Reducer(const topology::Topology & source_topology)
        : _source_topology(source_topology), _forcefield(), _forcefield_initialized(false), _rules(),
          _rules_initialized(false), _target_topology(), _ignore_duplicate_particles(false),
          _ignore_missing_particles(false)
    {
    }

    // =============================================================================
    // Accessors.
    // =============================================================================
    const topology::Topology & source_topology(void) const { return _source_topology; }
    const biospring::forcefield::ForceField & forcefield(void) const { return _forcefield; }
    const ReduceRuleContainer & rules(void) const { return _rules; }
    const topology::Topology & target_topology(void) const { return _target_topology; }

    void set_ignore_duplicate_particles(bool value) { _ignore_duplicate_particles = value; }
    bool get_ignore_duplicate_particles(void) const { return _ignore_duplicate_particles; }

    void set_ignore_missing_particles(bool value) { _ignore_missing_particles = value; }
    bool get_ignore_missing_particles(void) const { return _ignore_missing_particles; }

    // =============================================================================
    // Main methods.
    // =============================================================================

    void initialize_forcefield(const std::string & path);
    void initialize_rules(const std::string & path);

    void reduce(void)
    {
        if (!_forcefield_initialized)
            throw std::runtime_error("Reducer:: Forcefield not initialized!");

        if (!_rules_initialized)
            throw std::runtime_error("Reducer:: Rules not initialized!");

        logging::status("Reducing model using group and forcefield");

        ParticleGroups residues = _group_particles_by_residue(_source_topology.particles().data());

        for (const ParticleContainer & residue : residues)
        {
            ParticleContainer grains = _reduce_residue(residue);
            _target_topology.particles().push_back(grains);
        }
    }

    void reduce(const ReductionParameters & parameters)
    {
        // Dies if reduce parameter file is not set.
        if (parameters.pathForceField.empty())
            logging::die("Reducer: force field parameter file not provided");

        // Dies if force field parameter file is not set.
        if (parameters.pathGroup.empty())
            logging::die("Reducer: coarse grain parameter file not provided");

        if (not parameters.pathGroup.empty())
        {
            initialize_forcefield(parameters.pathForceField);
            initialize_rules(parameters.pathGroup);
            set_ignore_duplicate_particles(parameters.ignoreDuplicate);
            set_ignore_missing_particles(parameters.ignoreMissing);
            reduce();
        }
    }

  protected:
    std::vector<topology::Particle> _reduce_residue(const ParticleContainer & residue) const
    {
        const auto & rules = _rules.get_rules_for_residue(residue[0].properties().residue_name());

        // Returns an empty vector if found no rule for this residue.
        if (rules.empty())
        {
            std::string resname = residue[0].properties().residue_name();
            logging::warning("Residue '%s' not found in reduce rules. Skipping it.", resname.c_str());
            return {};
        }

        ResidueReducer residue_builder(residue, rules, _forcefield);
        residue_builder.set_ignore_duplicate_particles(_ignore_duplicate_particles);
        residue_builder.set_ignore_missing_particles(_ignore_missing_particles);

        return residue_builder.build();
    }

    static bool __still_same_residue(const ParticleContainer & current_residue, const topology::Particle & p)
    {
        return p.properties().residue_id() == current_residue[0].properties().residue_id() and
               p.properties().chain_name() == current_residue[0].properties().chain_name();
    }

    ParticleGroups _group_particles_by_residue(const ParticleContainer & particles) const
    {
        ParticleGroups residues;
        ParticleContainer current_residue = {particles[0]};

        for (size_t i = 1; i < particles.size(); i++)
        {
            const topology::Particle & p = particles[i];

            if (__still_same_residue(current_residue, p))
                current_residue.push_back(p);

            else
            {
                residues.push_back(current_residue);
                current_residue.clear();
                current_residue.push_back(p);
            }
        }
        residues.push_back(current_residue);
        return residues;
    }
};

namespace legacy
{

class Reducer
{
  public:
    Reducer() : _spn(0), _reduce(0), _ff(), _ignoreDuplicateParticles(false), _ignoreMissingParticle(false) {}

    void setIgnoreDuplicateParticles(bool value) { _ignoreDuplicateParticles = value; }
    void ignoreDuplicateParticles() { _ignoreDuplicateParticles = true; }
    bool getIgnoreDuplicateParticles() const { return _ignoreDuplicateParticles; }

    void setIgnoreMissingParticles(bool value) { _ignoreMissingParticle = value; }
    void ignoreMissingParticles() { _ignoreMissingParticle = true; }
    bool getIgnoreMissingParticles() const { return _ignoreMissingParticle; }

    spn::SpringNetwork * getSpringNetwork(void) const { return _spn; }
    void setSpringNetwork(spn::SpringNetwork * const spn) { _spn = spn; }

    Reduce * getReduce(void) const { return _reduce; }
    void setReduce(Reduce * const red) { _reduce = red; }
    void setReduce(const std::string & path);

    const biospring::forcefield::ForceField & getForceField(void) const { return _ff; }
    void setForceField(const std::string & path);

    void reduce(void) const;
    void reduceToCoarseGrain(spn::SpringNetwork & spn, const ReductionParameters & parameters);

    bool isSpringNetworkSet(void) const { return _spn->getNumberOfParticles() > 0; }
    bool isReduceSet(void) const { return _reduce->getNumberOfRules() > 0; }
    bool isForceFieldSet(void) const { return _ff.getNumberOfProperties() > 0; }

  protected:
    spn::SpringNetwork * _spn;
    Reduce * _reduce;
    biospring::forcefield::ForceField _ff;

    bool _ignoreDuplicateParticles;
    bool _ignoreMissingParticle;

    void _reduceAA(std::vector<const biospring::spn::Particle *> & particles,
                   std::vector<biospring::spn::Particle *> & cgParticles, const std::string & resname,
                   const size_t resid) const;
};

} // namespace legacy
} // namespace reduce
} // namespace biospring

#endif // __REDUCER_H__
