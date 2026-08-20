#ifndef __TOPOLOGY_TOPOLOGY_HPP__
#define __TOPOLOGY_TOPOLOGY_HPP__

#include <array>
#include <unordered_map>
#include <utility>

#include "Particle.hpp"
#include "ParticleCollection.hpp"
#include "Spring.hpp"
#include "SpringCollection.hpp"
#include "nsearch.hpp"

#include "SpringNetwork.h"
namespace biospring
{
namespace topology
{

// Topology class.
//
// A `Topology` is a collection of `Particle`s and `Spring`s.

class Topology
{
  protected:
    // The particles in the topology.
    ParticleCollection _particles;

    // The springs in the topology.
    SpringCollection _springs;

    // Ghost-spring collections for dihedral torsions, kept separate from
    // _springs (rather than tagging individual Spring objects) so each
    // family stays identifiable through serialization (.nc I/O) without
    // touching Spring/SpringCollection at all -- see
    // doc/BondedForceFieldSprings.md, Section 3.2. Which families are
    // actually built is decided once, at build time, by
    // BondedForceFieldReader::buildSprings's enableDihedral* parameters
    // (see -dihedralbackbone/-dihedralsidechain in pdb2spn-cli.cpp) -- PHI/
    // PSI/OMEGA are always built together under -dihedralbackbone (kept as
    // three separate collections here purely so SpringNetwork can enable/
    // disable each independently at runtime, see its dihedral.phi/psi/
    // omega .msp settings). A dihedral ghost spring is always a new
    // addition (never a retune of a real 1-2/1-3 spring), so unlike
    // _springs these are never populated by RigidBodyBuilder.
    // One collection per family, indexed by the family itself rather than
    // held as parallel members. Adding a family used to mean repeating the
    // same declaration/clear/copy/convert lines here and in four other
    // layers; indexed, everything below loops instead, and a new family is
    // one enum value.
    //
    // The family list itself is spn::SpringNetwork's -- these collections
    // are copied straight into it, family by family, by to_spring_network
    // below -- and is aliased here rather than restated: two independent
    // enums that have to agree in *order* behind a bare `unsigned family`
    // argument is exactly the kind of drift indexing them was meant to
    // remove.
    using DihedralFamilyIndex = spn::SpringNetwork::DihedralFamilyIndex;
    static constexpr unsigned DIHEDRAL_FAMILY_COUNT = spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT;
    std::array<SpringCollection, DIHEDRAL_FAMILY_COUNT> _dihedral_springs;

    // std::array has no "fill with N copies of this value" constructor, and
    // SpringCollection isn't default-constructible (it holds the particle
    // collection its springs index into), so the array can't just be
    // default-initialized in the member-init list. Built element by element
    // from the family count instead of restating one SpringCollection
    // (_particles) per family in each constructor.
    template <std::size_t... Is>
    static std::array<SpringCollection, DIHEDRAL_FAMILY_COUNT>
    _make_dihedral_springs(ParticleCollection & particles, std::index_sequence<Is...>)
    {
        return {(static_cast<void>(Is), SpringCollection(particles))...};
    }
    static std::array<SpringCollection, DIHEDRAL_FAMILY_COUNT> _make_dihedral_springs(ParticleCollection & particles)
    {
        return _make_dihedral_springs(particles, std::make_index_sequence<DIHEDRAL_FAMILY_COUNT>{});
    }

  public:
    // Binds a ghost (massless virtual-site) Particle -- created by
    // add_ghost_particle below, never 1:1 with a PDB atom -- to its 3 real
    // anchor particles (by unique id, not by index: unlike SpringNetwork's
    // particle vector, Topology's ParticleCollection can be reordered by
    // remove_particle, so a uid-keyed lookup via ParticleCollection::at_uid
    // stays correct even then) and its 3 calibrated placement parameters.
    // See spn::GhostParticle for the placement/force-redistribution
    // formulas this drives, applied once the topology is converted to a
    // SpringNetwork (see to_spring_network below).
    struct GhostParticleInfo
    {
        pid_t anchor_B_uid;
        pid_t anchor_C_uid;
        pid_t anchor_ref_uid;
        double r;
        double theta_deg;
        double delta_deg;
        // Which construction places this ghost -- see spn::GhostPlacement.
        // The two use disjoint parameters (axial reads r and ignores the
        // reference atom; rotation reads the reference atom and delta and
        // ignores r/theta), so this cannot be inferred from the values and
        // is carried explicitly, including through the .nc.
        unsigned placement;
    };

  protected:
    std::unordered_map<pid_t, GhostParticleInfo> _ghost_particles;

  public:
    // =============================================================================
    // Initialization methods.
    // =============================================================================

    Topology()
        : _springs(_particles), _dihedral_springs(_make_dihedral_springs(_particles))
    {
    }

    // Copy constructor.
    Topology(const Topology & other)
        : _springs(_particles), _dihedral_springs(_make_dihedral_springs(_particles))
    {
        _copy_particles(other);
        _copy_springs(other);
        _copy_ghost_particles(other);
    }

    // Assignment operator.
    Topology & operator=(const Topology & other)
    {
        // This preserves any springs that may already have been added when
        // parsing the input topology file (see pdb2spn and PDBPreader).
        // Useful, for example, for keeping the springs defined by the CONECT
        // lines of a coarse-grained topology generated by the martinize tool,
        // which is given as input to pdb2spn.
        // Reset the springs only if the number of particles in the new topology
        // after reduction is not equal to that in the input topology.
        // If the topology has not been reduced, for example an already
        // coarse-grained martini representation structure to which we simply
        // assign the names of the particles and their parameters via the --ff
        // and --grp options in pdb2spn (see the Martini ff and grp file in
        // /data/martini).
        if (_particles.size() != other._particles.size())
        {
            _springs.clear();
            for (auto & family : _dihedral_springs)
                family.clear();
        }

        _copy_particles(other);
        _copy_springs(other);
        _copy_ghost_particles(other);

        return *this;
    }

    // =============================================================================
    // Capacity.
    // =============================================================================

    size_t number_of_particles() const { return _particles.size(); }
    size_t number_of_springs() const { return _springs.size(); }

    // =============================================================================
    // Particle manipulation.
    // =============================================================================

    // Adds a particle to the topology.
    void add_particle(const Particle & particle) { _particles.push_back(particle); }

    // Adds a collection of particles to the topology.
    template <typename container> void add_particles(const container & particles) { _particles.push_back(particles); }
    void add_particles(const std::initializer_list<Particle> & particles) { _particles.push_back(particles); }

    // Reserves capacity for at least `n` particles. Growing the underlying
    // vector past its current capacity reallocates its buffer, which would
    // silently turn every Spring's cached Particle& (a real C++ reference,
    // not an index) into a dangling one (read back as garbage -- huge
    // unique_id, empty name -- the exact failure mode found and fixed by
    // reserving before --rigidbody/ghost-particle creation, see
    // BondedForceFieldReader and pdb2spn-cli.cpp's call site).
    //
    // A second, subtler instance of the same problem (found later, via a
    // real Fs-peptide PDB with CONECT records -- ubiquitin.pdb has none,
    // so it never exercised this path): a spring can already exist by the
    // time this is called even before --rigidbody ever runs, if the input
    // PDB had its own CONECT lines (parsed into real springs inside
    // PDBReader, long before pdb2spn-cli.cpp gets a chance to reserve
    // anything). A plain `_particles.data().reserve(n)` in that case would
    // still invalidate those already-existing springs' references.
    // Handled by detecting that case and rebuilding instead: snapshot the
    // current state into a temporary (a plain copy, always safe -- see the
    // copy constructor), clear this topology's own containers, reserve on
    // the now-empty (so reservation itself can never invalidate anything)
    // particle buffer, then replay the snapshot back in through the same
    // index-based re-resolution _copy_particles/_copy_springs/
    // _copy_ghost_particles already use for copying -- never lets
    // _particles grow again once anything actually references it.
    void reserve_particles(size_t n)
    {
        bool no_dihedral = true;
        for (const auto & family : _dihedral_springs)
            no_dihedral = no_dihedral && family.size() == 0;
        if (_springs.size() == 0 && no_dihedral && _ghost_particles.empty())
        {
            _particles.data().reserve(n);
            return;
        }

        Topology snapshot(*this);
        _particles.clear();
        _springs.clear();
        for (auto & family : _dihedral_springs)
            family.clear();
        _ghost_particles.clear();
        _particles.data().reserve(n);
        _copy_particles(snapshot);
        _copy_springs(snapshot);
        _copy_ghost_particles(snapshot);
    }

    // Adds a ghost (massless virtual-site) particle -- the first case
    // where Topology creates a particle with no 1:1 PDB atom behind it.
    // `particle` carries this ghost's own name/resname/residue_id/
    // chain_name (typically copied from `anchor_B`'s, since a ghost
    // conceptually belongs to the same residue as its anchors); its
    // position is not set here (to_spring_network/BondedForceFieldReader
    // compute it from the anchors via spn::GhostParticle::computePosition
    // once all 3 anchors are known to actually be resolved). Returns the
    // newly-added particle (with its own freshly-minted unique id, see
    // Particle::copy()/ParticleCollection::push_back).
    Particle & add_ghost_particle(const Particle & particle, const Particle & anchor_B, const Particle & anchor_C,
                                  const Particle & anchor_ref, double r, double theta_deg, double delta_deg,
                                  unsigned placement)
    {
        _particles.push_back(particle);
        Particle & added = _particles[_particles.size() - 1];
        _ghost_particles[added.unique_id()] =
            GhostParticleInfo{anchor_B.unique_id(), anchor_C.unique_id(), anchor_ref.unique_id(), r, theta_deg,
                              delta_deg, placement};
        return added;
    }

    bool is_ghost_particle(pid_t uid) const { return _ghost_particles.find(uid) != _ghost_particles.end(); }

    const GhostParticleInfo & get_ghost_particle_info(pid_t uid) const { return _ghost_particles.at(uid); }

    // Registers ghost-particle info for a particle that already exists (by index)
    // in the collection, without creating a new one. Used by NetCDFReader, which
    // reads ghost particles as regular particles first (see addParticlesToSpn())
    // and only learns which ones are ghosts, and their anchors, from a separate
    // buffer read afterwards.
    void register_ghost_particle(size_t particle_index, size_t anchor_B_index, size_t anchor_C_index,
                                  size_t anchor_ref_index, double r, double theta_deg, double delta_deg,
                                  unsigned placement)
    {
        pid_t particle_uid = _particles[particle_index].unique_id();
        _ghost_particles[particle_uid] = GhostParticleInfo{
            _particles[anchor_B_index].unique_id(), _particles[anchor_C_index].unique_id(),
            _particles[anchor_ref_index].unique_id(), r, theta_deg, delta_deg, placement};
    }

    // =============================================================================
    // Spring manipulation.
    // =============================================================================

    // Creates a spring between two particles.
    auto & add_spring(Particle & p1, Particle & p2, double equilibrium = -1.0, double stiffness = 1.0)
    {
        return _springs.add_spring(p1, p2, equilibrium, stiffness);
    }

    // Creates a spring between two particles using their position in the list of particles.
    auto & add_spring(size_t p1, size_t p2, double equilibrium = -1.0, double stiffness = 1.0)
    {
        return _springs.add_spring(_particles[p1], _particles[p2], equilibrium, stiffness);
    }

    // Creates a dihedral ghost spring between two real substituent atoms
    // (never a real 1-2 bond) -- one per DIHEDRAL entry, in the family the
    // entry names (a spn::SpringNetwork::DihedralFamilyIndex, resolved by
    // BondedForceFieldReader::buildSprings). See the _dihedral_springs
    // member comment above.
    auto & add_dihedral_spring(unsigned family, Particle & p1, Particle & p2, double equilibrium = -1.0,
                               double stiffness = 1.0)
    {
        // .at(), not [] -- see dihedral_springs() below for why the index is
        // range-checked.
        return _dihedral_springs.at(family).add_spring(p1, p2, equilibrium, stiffness);
    }

    // Adds springs between all particles within a cutoff distance.
    // Equilibrium distance is set to the distance bewteen each pair.
    void add_springs_from_cutoff(double cutoff, double stiffness = 1.0)
    {
        // Initializes the neighbor search object.
        nsearch::NeighborSearch nsearch(_particles, cutoff);

        // Iterates over all particles.
        for (Particle & p : _particles)
        {
            // Gets the neighbors of the current particle.
            auto neighbor_indexes = nsearch.get_neighbors(p);

            // Iterates over the neighbors.
            for (auto & neighbor_index : neighbor_indexes)
            {
                Particle & neighbor = _particles[neighbor_index];

                // If both particles are static, skip them.
                if (p.properties().is_static() && neighbor.properties().is_static())
                    continue;

                // Adds a spring between the current particle and the neighbor.
                // In this context, it's possible that we're trying to add a spring that
                // already exists.
                try
                {
                    auto & spring = add_spring(p, neighbor);
                    spring.set_stiffness(stiffness);
                }
                catch (const SpringAlreadyExistsException & e)
                {
                }
            }
        }
    }

    // Adds springs betweens particles that originate from different topologies.
    // This function's main purpose is to add springs between particles during
    // topology merging.
    void add_springs_between_topologies_from_cutoff(double cutoff, double stiffness = 1.0)
    {
        // Initializes the neighbor search object.
        nsearch::NeighborSearch nsearch(_particles, cutoff);

        // Iterates over all particles.
        for (Particle & p : _particles)
        {
            // Gets the neighbors of the current particle.
            auto neighbor_indexes = nsearch.get_neighbors(p);

            // Iterates over the neighbors.
            for (auto & neighbor_index : neighbor_indexes)
            {
                Particle & neighbor = _particles[neighbor_index];

                // If the neighbor is in the same topology, skip it.
                if (p.properties().topology_id() == neighbor.properties().topology_id())
                    continue;

                // If both particles are static, skip them.
                if (p.properties().is_static() && neighbor.properties().is_static())
                    continue;

                // Adds a spring between the current particle and the neighbor.
                // In this context, it's possible that we're trying to add a spring that
                // already exists.
                try
                {
                    auto & spring = add_spring(p, neighbor);
                    spring.set_stiffness(stiffness);
                }
                catch (const SpringAlreadyExistsException & e)
                {
                }
            }
        }
    }

    // =============================================================================
    // Particle Getters.
    // =============================================================================

    auto & particles() { return _particles; }
    const auto & particles() const { return _particles; }

    auto & get_particle(size_t position) { return _particles[position]; }
    const auto & get_particle(size_t position) const { return _particles[position]; }

    auto & get_particle_at_uid(pid_t uid) { return _particles.at_uid(uid); }
    const auto & get_particle_at_uid(pid_t uid) const { return _particles.at_uid(uid); }

    // =============================================================================
    // Spring Getters.
    // =============================================================================

    auto & springs() { return _springs; }
    const auto & springs() const { return _springs; }

    // One dihedral family's collection, by
    // spn::SpringNetwork::DihedralFamilyIndex. Range-checked: a family
    // index is routinely computed (parsed from a .bi.ff, walked over in a
    // loop) rather than written out literally, so an out-of-range one is a
    // real possibility and reads as garbage rather than failing.
    auto & dihedral_springs(unsigned family) { return _dihedral_springs.at(family); }
    const auto & dihedral_springs(unsigned family) const { return _dihedral_springs.at(family); }

    auto & get_spring(size_t position) { return _springs[position]; }
    const auto & get_spring(size_t position) const { return _springs[position]; }

    auto & get_spring(pid_t uid1, pid_t uid2) { return _springs.at_uid({uid1, uid2}); }
    const auto & get_spring(pid_t uid1, pid_t uid2) const { return _springs.at_uid({uid1, uid2}); }

    auto & get_spring(const Particle & p1, const Particle & p2)
    {
        return _springs.at_uid({p1.unique_id(), p2.unique_id()});
    }

    auto & get_spring_at_uid(sid_t uid) { return _springs.at_uid(uid); }
    const auto & get_spring_at_uid(sid_t uid) const { return _springs.at_uid(uid); }

    // =============================================================================
    // Spring Lookup.
    // =============================================================================

    bool has_spring_between(const Particle & p1, const Particle & p2) const { return _springs.exists(p1, p2); }

    // =============================================================================
    // Merge.
    // =============================================================================

    // Returns a new Topology which is the result of merging two topologies.
    // The new topology will contain all particles and springs from both topologies.
    static Topology merge(const Topology & top1, const Topology & top2)
    {
        Topology result;
        result._particles += top1._particles + top2._particles;
        result._copy_springs(top1, 0);
        result._copy_springs(top2, top1.number_of_particles());
        return result;
    }

    Topology merge(const Topology & other) const { return merge(*this, other); }

    // =============================================================================
    // Conversion.
    // =============================================================================

    // Converts the topology to a SpringNetwork.
    // The result is passed as argument because SpringNetwork has no copy constructor.
    void to_spring_network(spn::SpringNetwork & spn) const
    {
        // Removes all springs and particles from the SpringNetwork.
        spn.clear();

        // Copies particles. Ghost particles (see add_ghost_particle) are
        // never copied via the regular spn::Particle path: their 3 anchors
        // must already exist on the SpringNetwork side (true here, since
        // they were always added earlier in `_particles` -- a ghost's
        // anchors are real PDB atoms, always created before any ghost
        // referencing them), so their spn-side indices are resolved via
        // the same by_uid() lookup already used for springs below.
        for (const topology::Particle & source : _particles)
        {
            if (is_ghost_particle(source.unique_id()))
            {
                const GhostParticleInfo & info = get_ghost_particle_info(source.unique_id());
                const unsigned anchor_b_index = static_cast<unsigned>(_particles.by_uid().at(info.anchor_B_uid));
                const unsigned anchor_c_index = static_cast<unsigned>(_particles.by_uid().at(info.anchor_C_uid));
                const unsigned anchor_ref_index = static_cast<unsigned>(_particles.by_uid().at(info.anchor_ref_uid));
                spn.addGhostParticle(info.placement, anchor_b_index, anchor_c_index, anchor_ref_index,
                                     static_cast<float>(info.r),
                                    static_cast<float>(info.theta_deg), static_cast<float>(info.delta_deg));
                continue;
            }

            spn::Particle target;
            target.setName(source.properties().name());
            target.setResName(source.properties().residue_name());
            // residue_id()/atom_id() stay signed ints (some PDB files use
            // negative residue numbering), while the simulation-side
            // Particle assumes a non-negative id; cast explicitly here.
            target.setResId(static_cast<unsigned>(source.properties().residue_id()));
            target.setChainName(source.properties().chain_name());
            target.setExtid(static_cast<unsigned>(source.properties().atom_id()));
            target.setMass(source.properties().mass());
            target.setCharge(source.properties().charge());
            target.setRadius(source.properties().radius());
            // epsilon was missing here, and it silently disabled the whole
            // steric term: it stayed at spn::ParticleProperty's 0.0 default,
            // and every Lennard-Jones form is proportional to epsilon_ij, so
            // both the steric energy AND force came out exactly zero on any
            // topology built this way (found 2026-08-13 while writing
            // example/072: "Steric energy: 0.00" with correct radii/charges
            // in the same .nc). Every other link in the chain was already
            // correct -- ForceFieldReader parses it, Reducer sets it on the
            // grain, SpnBuffer writes it, NetCDFReader restores it -- only
            // this topology-to-network hand-off dropped it.
            target.setEpsilon(source.properties().epsilon());
            target.setPosition(source.properties().position());
            target.setVelocity(source.properties().velocity());
            target.setForce(source.properties().force());
            target.setDynamic(source.properties().is_dynamic());
            target.setHydrophobicity(source.properties().hydrophobicity());
            target.setSolventAccessibilitySurface(source.properties().imp().solvent_accessible_surface());
            target.setTransferEnergyByAccessibleSurface(source.properties().imp().transfert_energy_by_accessible_surface());

            spn.addParticle(target);
        }

        // Copies springs.
        for (const topology::Spring & source : _springs)
        {
            // Retrieves the positions of the particles making up the spring.
            size_t i = _particles.by_uid().at(source.first().unique_id());
            size_t j = _particles.by_uid().at(source.second().unique_id());

            spn.addSpring(i, j, source.equilibrium(), source.stiffness());
        }

        // Copies dihedral ghost springs, one family at a time (see the
        // _dihedral_springs member comment above for why these stay
        // separate from _springs) -- into the same family index on the
        // network side, which is the very index these are stored under.
        for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
        {
            for (const topology::Spring & source : _dihedral_springs[family])
            {
                size_t i = _particles.by_uid().at(source.first().unique_id());
                size_t j = _particles.by_uid().at(source.second().unique_id());
                spn.addDihedralSpring(family, i, j, source.equilibrium(), source.stiffness(), source.dc_offset());
            }
        }
    }

  protected:
    // Copies `other`'s particles to this topology.
    void _copy_particles(const Topology & other)
    {
        // Copies particles.
        _particles = other._particles;
    }

    // Copies one spring collection (`src`, belonging to `src_particles`)
    // into `dst` (belonging to `this->_particles`), shifting particle
    // positions by `offset`. Shared by _copy_springs below for every spring
    // collection a Topology carries.
    void _copy_spring_collection(SpringCollection & dst, const SpringCollection & src,
                                  const ParticleCollection & src_particles, size_t offset)
    {
        for (const Spring & source : src)
        {
            size_t i = offset + src_particles.by_uid().at(source.first().unique_id());
            size_t j = offset + src_particles.by_uid().at(source.second().unique_id());

            // Two merged topologies may share overlapping particles (e.g. the
            // boundary residues between two structures reduced separately
            // before merging), so the same spring can already exist.
            try
            {
                dst.add_spring(_particles[i], _particles[j], source.equilibrium(), source.stiffness())
                    .set_dc_offset(source.dc_offset());
            }
            catch (const SpringAlreadyExistsException & e)
            {
            }
        }
    }

    // Copies `other`'s springs (all four collections) to this topology.
    // `offset` is the position at which `other`'s particles begin in
    // `this` collection: 0 when `this` is an exact copy of `other` (copy
    // constructor/assignment operator, same particles in the same order),
    // or the particle count already present in `this` when appending
    // `other` as part of a merge.
    //
    // Note this can't be resolved via unique_id() instead of a plain
    // offset: Particle::copy() (used by ParticleCollection::push_back(),
    // which is what `_particles +=` calls into) assigns each copy a *new*
    // unique id, so by the time this runs, `this->_particles` no longer
    // shares any unique_id with `other._particles` at all -- `other`'s own
    // unique_id lookup only ever finds `other`'s own *local* position.
    void _copy_springs(const Topology & other, size_t offset = 0)
    {
        _copy_spring_collection(_springs, other._springs, other._particles, offset);
        for (unsigned family = 0; family < DIHEDRAL_FAMILY_COUNT; ++family)
            _copy_spring_collection(_dihedral_springs[family], other._dihedral_springs[family], other._particles,
                                    offset);
    }

    // Copies `other`'s ghost-particle bindings to this topology. Same uid
    // problem as _copy_springs above (see its note): `other`'s ghost/anchor
    // unique ids only resolve to a *position* in `other._particles`, which
    // must be re-looked-up in `this->_particles` (already copied, same
    // order, by _copy_particles) to get the corresponding fresh unique id.
    void _copy_ghost_particles(const Topology & other)
    {
        _ghost_particles.clear();
        for (const auto & [old_ghost_uid, info] : other._ghost_particles)
        {
            size_t ghost_index = other._particles.by_uid().at(old_ghost_uid);
            size_t b_index = other._particles.by_uid().at(info.anchor_B_uid);
            size_t c_index = other._particles.by_uid().at(info.anchor_C_uid);
            size_t ref_index = other._particles.by_uid().at(info.anchor_ref_uid);

            pid_t new_ghost_uid = _particles[ghost_index].unique_id();
            _ghost_particles[new_ghost_uid] =
                GhostParticleInfo{_particles[b_index].unique_id(), _particles[c_index].unique_id(),
                                  _particles[ref_index].unique_id(), info.r, info.theta_deg, info.delta_deg,
                                  info.placement};
        }
    }
};

} // namespace topology
} // namespace biospring

#endif // __TOPOLOGY_TOPOLOGY_HPP__