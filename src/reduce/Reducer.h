
//
// Defines the reducer class which reduces an all atom topology
// to a coarse grain topology.
//

#ifndef __REDUCER_H__
#define __REDUCER_H__

#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Particle.h"
#include "ReduceRuleContainer.hpp"
#include "SpringNetwork.h"
#include "forcefield/ForceField.h"
#include "logging.h"
#include "measure.hpp"
#include "topology.hpp"
#include "topology/TerminalPosition.hpp"
#include "utils/string.hpp"

namespace biospring
{
namespace reduce
{

using ParticleContainer = std::vector<topology::Particle>;

// A grain, together with the unique ids of the source particles it was built
// from. The provenance is what lets springs survive the reduction: a spring
// joins two source ATOMS, and after reduction it must join the two GRAINS
// those atoms ended up in (see Reducer::_carry_springs_over).
//
// The uids are usable as a key because a plain Particle copy preserves
// unique_id() -- only ParticleCollection::push_back, which goes through
// Particle::copy(), issues a new one. So they still identify particles of
// `_source_topology` here, and can never collide with the target's own.
struct ReducedGrain
{
    topology::Particle grain;
    std::vector<topology::pid_t> member_uids;
};

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
    Grain _grain;                                // output grain
    std::vector<topology::pid_t> _member_uids;   // source particles it was built from
    bool _grain_built;
    bool _warn_when_empty;
    // Set from _check_not_empty, which is const.
    mutable bool _was_empty;
    // Set from _create_grain, which is const.
    mutable std::string _forcefield_fallback;

  public:
    GrainBuilder(const ParticleContainer & input_particles, const ReduceRule & rule,
                 const forcefield::ForceField & forcefield)
        : _input_particles(input_particles), _rule(rule), _forcefield(forcefield), _ignore_duplicate_particles(false),
          _ignore_missing_particles(false), _grain(), _member_uids(), _grain_built(false), _warn_when_empty(true),
          _was_empty(false), _forcefield_fallback()
    {
    }

    // =============================================================================
    // Accessors.
    // =============================================================================
    void set_ignore_duplicate_particles(bool value) { _ignore_duplicate_particles = value; }
    bool get_ignore_duplicate_particles(void) const { return _ignore_duplicate_particles; }

    void set_ignore_missing_particles(bool value) { _ignore_missing_particles = value; }
    bool get_ignore_missing_particles(void) const { return _ignore_missing_particles; }

    // Several rules may build the SAME grain from different spellings of the
    // same atom -- amber.dna.grp accepts both OP1 and O1P for dOP1, both H5'
    // and H5'1 for dH5. Exactly one spelling can match, so the others find
    // nothing, and warning from here announced as "skipped" a grain the
    // sibling rule had just built. ResidueReducer::build silences this one and
    // warns once per residue, only for grains that no rule produced at all.
    void set_warn_when_empty(bool value) { _warn_when_empty = value; }

    // Whether build() failed because no atom matched, as opposed to matching
    // the wrong number of them.
    bool was_empty(void) const { return _was_empty; }

    // The force-field key this grain asked for and did not get, having been
    // parameterised from the unqualified one instead; empty when nothing fell
    // back. Reported by the caller across the whole model -- see
    // Reducer::_report_forcefield_fallbacks.
    const std::string & forcefield_fallback(void) const { return _forcefield_fallback; }

    const Grain & grain(void)
    {
        if (!_grain_built)
            throw std::runtime_error("Grain not built yet!");
        return _grain;
    }

    // Unique ids of the source particles this grain aggregates.
    const std::vector<topology::pid_t> & member_uids(void) const { return _member_uids; }

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

        _member_uids.clear();
        _member_uids.reserve(grain_particles.size());
        for (const topology::Particle & p : grain_particles)
            _member_uids.push_back(p.unique_id());

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
        // A grain built from ONE source particle is a renaming, and this is the
        // name it replaced -- taken from the RULE, not from the structure. A
        // type may accept several spellings (amber.grp has "XH1 ACE H1" and
        // "XH1 ACE HH31", "ZC NME C" and "ZC NME CH3"), and a parameter table
        // is written against one of them, so recording whichever spelling this
        // particular PDB happened to use would match that table only by luck.
        //
        // A grain built from several particles is a real average and has no
        // single original, so it keeps none.
        //
        // This is also why no heuristic is needed anywhere: the question is
        // answered HERE, with the .grp in hand, and the answer is written into
        // the .nc. Asking it again at run time is what forced a guess.
        if (grain_particles.size() == 1)
            grain.properties().set_original_name(_rule.getCanonicalAtomName());
        grain.properties().set_residue_name(_residue_name());
        grain.properties().set_residue_id(_residue_id());
        grain.properties().set_chain_name(_input_particles[0].properties().chain_name());
        grain.properties().set_position(measure::centroid(grain_particles));

        // Set properties from forcefield.
        const std::string key = _forcefield_key();
        if (_forcefield.hasProperty(key))
        {
            const spn::ParticleProperty & pp = _forcefield.getPropertiesFromName(key);
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
            logging::warning("No Forcefield found for '%s'!", key.c_str());
        }

        return grain;
    }

    // Which force-field type this grain is parameterised from.
    //
    // A terminal rule asks for the qualified type (see
    // ReduceRule::forcefield_key), and a force field that does not carry it
    // falls back on the unqualified one. That keeps a .grp whose terminal rules
    // are newer than the .ff it is paired with at today's behaviour -- slightly
    // wrong charges -- rather than at a particle with no mass, which is what a
    // missing type gives. It is reported, once per key, by the caller.
    std::string _forcefield_key() const
    {
        const std::string key = _rule.forcefield_key();
        if (key == _rule.name() || _forcefield.hasProperty(key))
            return key;

        _forcefield_fallback = key;
        return _rule.name();
    }

    // Checks if a grain is empty.
    // If so, displays a warning message.
    // Returns the success status, which is true if grain is not empty, false otherwise.
    bool _check_not_empty(const ParticleContainer & grain) const
    {
        bool success = true;
        if (grain.empty())
        {
            _was_empty = true;
            if (_warn_when_empty)
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

    // Which grain names this residue produced and which its rules left empty.
    // Reported by the caller across the WHOLE model, not here -- see build().
    std::set<std::string> _built;
    std::set<std::string> _empty;

    // Force-field keys the rules asked for and the force field did not carry,
    // also reported by the caller -- see GrainBuilder::_forcefield_key.
    std::set<std::string> _forcefield_fallbacks;

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

    // Valid after build(): the grain names this residue produced, and those its
    // rules found no atom for.
    const std::set<std::string> & grains_built(void) const { return _built; }
    const std::set<std::string> & grains_left_empty(void) const { return _empty; }

    // Valid after build(): the qualified force-field keys this residue's rules
    // asked for in vain.
    const std::set<std::string> & forcefield_fallbacks(void) const { return _forcefield_fallbacks; }

    // =============================================================================
    // Main methods.
    // =============================================================================

    // Returns a vector of grains representing the residue, each carrying the
    // source particles it was built from.
    std::vector<ReducedGrain> build(void)
    {
        // Displays warnings if some particles are present in the residue but not in the rule.
        _check_unknown_particles();

        std::vector<ReducedGrain> grains;

        // A grain name can appear in several rules, one per accepted spelling
        // of its atoms (amber.dna.grp takes OP1 or O1P for dOP1, H5' or H5'1
        // for dH5). Only one spelling can be present, so the others necessarily
        // find nothing -- which is not a missing grain and must not be reported
        // as one. Decide per grain NAME, after every rule has had its turn.
        std::set<std::string> & built = _built;
        std::set<std::string> & empty = _empty;
        built.clear();
        empty.clear();
        _forcefield_fallbacks.clear();

        for (const auto & rule : _rules)
        {
            GrainBuilder grain_builder(_particles, rule, _forcefield);
            grain_builder.set_ignore_duplicate_particles(_ignore_duplicate_particles);
            grain_builder.set_ignore_missing_particles(_ignore_missing_particles);
            grain_builder.set_warn_when_empty(false);

            bool success = grain_builder.build();
            if (success)
            {
                grains.push_back({grain_builder.grain(), grain_builder.member_uids()});
                built.insert(rule.name());
            }
            else if (grain_builder.was_empty())
            {
                empty.insert(rule.name());
            }

            if (!grain_builder.forcefield_fallback().empty())
                _forcefield_fallbacks.insert(grain_builder.forcefield_fallback());
        }

        // Nothing is reported here. A grain whose rules found no atom in THIS
        // residue is the normal case for a rule that only applies at one end of
        // a chain: amber.grp declares the N-terminal ammonium hydrogens and the
        // C-terminal second oxygen for all 20 residues, so every interior
        // residue leaves four of them empty. Reported per residue, that was 42%
        // of the whole example suite's build output -- 10334 of 076's 10434
        // lines -- which buried the four warnings that mattered there.
        //
        // What a per-residue report cannot distinguish, and a model-wide one
        // can, is a rule that matches NOWHERE: a typo in a .grp entry. So the
        // two sets are handed to the caller, which reports the names left empty
        // across the whole model and built in none of it. See
        // Reducer::_report_grains_never_built.
        for (const std::string & name : built)
            empty.erase(name);

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

    // Accumulated over every residue of the model, so that a grain name can be
    // reported once if no residue anywhere produced it -- see
    // _report_grains_never_built.
    std::set<std::string> _grains_built;
    std::set<std::string> _grains_empty;

    // Likewise for the qualified force-field types a terminal rule asked for
    // and the force field did not carry -- see _report_forcefield_fallbacks.
    std::set<std::string> _forcefield_fallbacks;

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

    // Valid after reduce(): how many (residue, grain) rules found no atom in
    // any residue of the whole model. Expected to be nonzero as soon as the
    // rules carry position-dependent grains, and not a defect -- see
    // _report_grains_never_built for why, and for what is one.
    size_t grains_never_built(void) const
    {
        size_t never = 0;
        for (const std::string & key : _grains_empty)
            if (_grains_built.find(key) == _grains_built.end())
                never++;
        return never;
    }

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

        // Which (residue, position) pairs the .grp carries, so a position it
        // says nothing about falls back to the interior rules -- which is every
        // position of every .grp but amber.grp, hence unchanged behaviour
        // there. See topology::wanted_position.
        const auto declared = _rules.declared_positions();

        // Where each source particle ended up. Built as the grains are added
        // so the target index is simply the collection's size at that point.
        std::unordered_map<topology::pid_t, size_t> source_to_target;

        for (size_t index = 0; index < residues.size(); index++)
        {
            const ParticleContainer & residue = residues[index];

            // Read off the chain the same way RigidBodyBuilder and
            // BondedForceFieldReader read "+"/"-": a residue is a terminus when
            // the neighbouring group in particle order belongs to another
            // chain, or does not exist. So a chain BROKEN in the middle -- a
            // gap in the numbering, a missing loop -- is not seen as a
            // terminus, exactly as a "+N" there still resolves to the next
            // residue present.
            const std::string & chain = residue[0].properties().chain_name();
            const bool is_nterm = index == 0 || residues[index - 1][0].properties().chain_name() != chain;
            const bool is_cterm =
                index + 1 == residues.size() || residues[index + 1][0].properties().chain_name() != chain;

            // And being at the end of a chain is not the same as being a free
            // terminus: a carved-out fragment has neither the ammonium nor the
            // carboxylate. See topology::is_nterminal_only_atom.
            bool has_nterminal_atom = false;
            bool has_cterminal_atom = false;
            for (const topology::Particle & p : residue)
            {
                has_nterminal_atom |= topology::is_nterminal_only_atom(p.properties().name());
                has_cterminal_atom |= topology::is_cterminal_only_atom(p.properties().name());
            }

            for (const ReducedGrain & reduced :
                 _reduce_residue(residue, is_nterm, is_cterm, has_nterminal_atom, has_cterminal_atom, declared))
            {
                const size_t target_index = _target_topology.number_of_particles();
                _target_topology.add_particle(reduced.grain);

                // An atom claimed by two rules producing different grains is
                // ambiguous; the first grain wins, which is arbitrary but
                // stable. emplace() keeps the first, insert_or_assign would not.
                for (topology::pid_t uid : reduced.member_uids)
                    source_to_target.emplace(uid, target_index);
            }
        }

        _report_grains_never_built();
        _report_forcefield_fallbacks();
        _carry_springs_over(source_to_target);
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
    // Re-draws the source topology's springs between the grains its atoms
    // ended up in.
    //
    // The source's springs are the ones parsed from the input file -- in
    // practice the CONECT records of a coarse-grained PDB, which Topology's
    // assignment operator goes out of its way to preserve. They are stored as
    // references to source particles, and the reduction replaces the particle
    // list wholesale: grains come out in rule order, which has no relation to
    // the input's atom order. Carrying them over by position, as used to
    // happen, left every one of them on an unrelated pair of grains --
    // silently, because a spring's stored equilibrium stays whatever it was.
    //
    // The equilibrium is recomputed (-1.0 = "use the current distance")
    // rather than carried over: a grain sits at its members' centroid, so the
    // atom-atom distance the source spring stored is not the grain-grain
    // distance. When the reduction is a 1:1 renaming -- the martinize case,
    // and every all-atom .grp here -- the centroid of a single atom is that
    // atom, and the recomputed value is the stored one.
    void _carry_springs_over(const std::unordered_map<topology::pid_t, size_t> & source_to_target)
    {
        if (_source_topology.number_of_springs() == 0)
            return;

        for (const topology::Spring & source : _source_topology.springs())
        {
            const auto first = source_to_target.find(source.first().unique_id());
            const auto second = source_to_target.find(source.second().unique_id());

            // An endpoint no rule claimed, or both endpoints swallowed by the
            // same grain: there is no spring left to draw.
            if (first == source_to_target.end() || second == source_to_target.end() ||
                first->second == second->second)
                continue;

            try
            {
                _target_topology.add_spring(first->second, second->second, -1.0, source.stiffness());
            }
            // Two source springs can collapse onto the same pair of grains.
            catch (const topology::SpringAlreadyExistsException &)
            {
            }
        }

        logging::info("Kept %zu of the input topology's %zu springs through reduction",
                      _target_topology.number_of_springs(), _source_topology.number_of_springs());
    }

    std::vector<ReducedGrain>
    _reduce_residue(const ParticleContainer & residue, bool is_nterm, bool is_cterm, bool has_nterminal_atom,
                    bool has_cterminal_atom,
                    const std::set<std::pair<std::string, topology::TerminalPosition>> & declared)
    {
        const std::string & residue_name = residue[0].properties().residue_name();

        // A free chain end is AMBER's own residue: NALA, not ALA with three
        // hydrogens bolted on. The rules qualified for this position replace
        // the unqualified ones rather than adding to them, and a position the
        // .grp says nothing about falls back to them -- see
        // topology::wanted_position.
        const topology::TerminalPosition position = topology::wanted_position(
            residue_name, is_nterm, is_cterm, has_nterminal_atom, has_cterminal_atom, declared);

        const auto & rules = _rules.get_rules_for_residue(residue_name, position);

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

        std::vector<ReducedGrain> grains = residue_builder.build();

        // A grain name is keyed by residue as well: ALA's AOXT and GLY's GOXT
        // are different rules, and one being built says nothing about the
        // other. Without the residue name, one terminal alanine in a model
        // would silence every other residue's terminal grains.
        const std::string & resname = residue[0].properties().residue_name();
        for (const std::string & name : residue_builder.grains_built())
            _grains_built.insert(resname + ":" + name);
        for (const std::string & name : residue_builder.grains_left_empty())
            _grains_empty.insert(resname + ":" + name);

        for (const std::string & key : residue_builder.forcefield_fallbacks())
            _forcefield_fallbacks.insert(key);

        return grains;
    }

    // How many rule grains the model never produced -- one line, not a warning.
    //
    // It is not a defect, and that is the whole point of counting it here
    // instead of reporting it per residue. A positional rule only matches if
    // the model happens to put that residue type at that position: amber.grp
    // declares the N-terminal ammonium and the C-terminal OXT for all 20
    // residues, so a four-chain protein can satisfy at most four of each and
    // the other sixteen types necessarily match nothing. There is nothing to
    // check and nothing is lost.
    //
    // What IS a defect is the dual, and it is already reported, per atom, by
    // _check_unknown_particles: an atom no rule claims is silently dropped from
    // the network (072 lost a third of its protein that way). That check also
    // covers the typo this one looks like it covers -- misspell ALA's AH1 rule
    // and the real H1 of a terminal alanine goes unclaimed, which it reports.
    void _report_grains_never_built(void) const
    {
        const size_t never = grains_never_built();
        if (never)
            logging::info("%zu reduce-rule grain(s) matched no atom anywhere, which is expected of a rule "
                          "that only applies at a chain end; an atom no rule claims is reported separately",
                          never);
    }

    // Which qualified force-field types the .grp asked for in vain -- one line,
    // naming them, because this one IS a mismatch worth seeing.
    //
    // It means the .grp's terminal rules are newer than the .ff paired with it:
    // the grains were parameterised from the interior types instead, so a
    // terminal residue carries the interior net charge (-1.2719 e at an
    // N-terminus, +1.0000 e at a C-terminus, where AMBER has +1 and -1) rather
    // than no mass at all, which is what a missing type would give.
    void _report_forcefield_fallbacks(void) const
    {
        if (_forcefield_fallbacks.empty())
            return;

        logging::warning("%zu terminal force-field type(s) are missing from the force field, so those grains kept "
                         "the interior parameters -- the terminal residues' net charge is AMBER's interior one. "
                         "Pair the .grp with a force field that declares them (amber.nbi.ff does).",
                         _forcefield_fallbacks.size());
        logging::info("Missing terminal force-field type(s): %s",
                      utils::string::join(std::vector<std::string>(_forcefield_fallbacks.begin(),
                                                                   _forcefield_fallbacks.end()),
                                          " ")
                          .c_str());
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

} // namespace reduce
} // namespace biospring

#endif // __REDUCER_H__
