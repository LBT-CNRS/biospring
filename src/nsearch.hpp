// Neighbor search.
//
// A cell list over the particles, queried one particle at a time. The grid is
// rebuilt by `update()`, which a positive skin lets skip until something has
// actually drifted; `for_each_neighbor` then answers from the current
// positions.
//
// A cell is NOT the size of the cutoff -- see forcefield/shared/cellgrid_shared.h
// for what that buys and what picks the width. Callers that pass no width get
// one cell per search radius and the 3x3x3 stencil this started out with.
//
// Example:
//
//     nsearch::NeighborSearch search(particles, cutoff);
//     for (size_t j : search.get_neighbors(particles[0]))
//         ...
//
//     search.update();  // after the particles have moved
//
// `NeighborSearchO2` is the same query answered by brute force. It exists to
// hold the cell list to account in the tests, and is too slow for anything else.

#ifndef __NSEARCH_HPP__
#define __NSEARCH_HPP__

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

#include "box.hpp"
#include "concepts.hpp"
#include "forcefield/shared/cellgrid_shared.h"
#include "measure.hpp"

namespace biospring
{
namespace nsearch
{

// Whether keeping a neighbour list is paying for itself, decided from what
// actually happens rather than from a setting.
//
// The criterion is REUSE: how many steps a list served before it had to be
// rebuilt. That is a complete criterion here and only here, which is worth
// spelling out, because the device uses no such policy and the asymmetry looks
// arbitrary until the three things a list buys are separated:
//
//   A. the search amortised over the steps the list survives;
//   B. a contiguous read where the cell walk chases a linked list;
//   C. only the particles a term acts on -- the charged ones for Coulomb.
//
// Only A depends on reuse. This class measures A, so it is the whole story
// exactly when B and C are nil, and on this path they are: C the CPU has had
// since before there was a device, through the per-term searchers, and B is
// small because an out-of-order core with a large cache absorbs a pointer chase.
// Measured on 023 with a 1 A skin -- a structure whose fastest beads force a
// rebuild at every step, so A = 0 by construction -- the list costs 7% here and
// saves 59% on the device. Same structure, same motion, same absent reuse: what
// differs is B and C, and there they are worth 40% and 58% on their own.
//
// So: a judgement about THIS backend's cost structure, not about the
// simulation. Do not lift it to the device without measuring there.
//
// Why it is needed at all: a list only earns its cost by being reused, and if
// the structure crosses half the skin every step it is rebuilt every step --
// two walks where the plain cell search does one. On 023 a 1 A skin ran 2.4x
// slower than no list at all before this existed. It is tried again later,
// because a structure changes regime: one that is equilibrating slows down.
//
// (On 023 the beads that force it are hydroxyl hydrogens moving at 84 km/s,
// thirty times thermal -- the known collapse of an AMBER hydroxyl H, which has
// no Lennard-Jones while BioSpring excludes only SPRUNG pairs from Coulomb. The
// list is not what is wrong there.)
class ListPolicy
{
    // How many steps a list has to survive to have been worth building. Two:
    // one rebuild serving one step is exactly break-even at best, since the
    // build walks the same cells the search would have walked.
    static constexpr unsigned MINIMUM_SURVIVAL = 2;
    // How many useless rebuilds to sit through before giving up. Not one: a
    // single fast step -- an interaction, a clash relaxing -- should not cost
    // the list for the rest of the run.
    static constexpr unsigned PATIENCE = 3;
    // And how long before trying again.
    static constexpr unsigned RETRY_AFTER = 1000;

    unsigned _sincerebuild = MINIMUM_SURVIVAL; // so the first build is not judged
    unsigned _wasted = 0;
    unsigned _sinceretry = 0;
    bool _worthit = true;

  public:
    bool worth_it() const { return _worthit; }
    unsigned wasted_rebuilds() const { return _wasted; }

    // Once per step, whether or not anything is rebuilt.
    void step()
    {
        _sincerebuild++;
        if (_worthit)
            return;
        if (++_sinceretry >= RETRY_AFTER)
        {
            _worthit = true;
            _wasted = 0;
            _sinceretry = 0;
        }
    }

    // When a rebuild actually happens.
    void rebuilt()
    {
        if (_sincerebuild < MINIMUM_SURVIVAL)
        {
            if (++_wasted >= PATIENCE)
                _worthit = false;
        }
        else
            _wasted = 0;
        _sincerebuild = 0;
    }
};

template <concepts::LocatableContainer ContainerType> class NeighborSearchBase
{
  protected:
    // The Locatable list.
    const ContainerType * _system;

    // The cutoff distance.
    float _cutoff;

  public:
    // Initializes the neighbor search object with the particles.
    NeighborSearchBase(const ContainerType & container, float cutoff) : _system(&container), _cutoff(cutoff)
    {
        if (container.empty())
            throw std::invalid_argument("the particle list is empty");

        if (_cutoff < 1e-6f)
            throw std::invalid_argument("the cutoff distance must be positive and different from zero");
    }
};

// The `NeighborSearchO2` class is a neighbor search implementation
// that uses the O(N^2) algorithm. It is used as a reference
// implementation to test the `NeighborSearch` class.
template <concepts::LocatableContainer ContainerType> class NeighborSearchO2 : public NeighborSearchBase<ContainerType>
{
  protected:
    using NeighborSearchBase<ContainerType>::_cutoff;
    using NeighborSearchBase<ContainerType>::_system;

  public:
    // Initializes the neighbor search object with the particles.
    NeighborSearchO2(const ContainerType & container, float cutoff)
        : NeighborSearchBase<ContainerType>(container, cutoff)
    {
    }

    // Finds and returns the neighbors of the given element.
    template <concepts::Locatable T> std::vector<size_t> get_neighbors(const T & element) const
    {
        // Loops over the particles.
        std::vector<size_t> neighbors;
        for (size_t i = 0; i < _system->size(); i++)
        {
            const T & candidate = _system->at(i);

            // If candidate if not `element` itself and is within the cutoff distance, add it to the list of neighbors.
            if (&candidate != &element && measure::distance(element, candidate) < _cutoff)
                neighbors.push_back(i);
        }
        return neighbors;
    }
};

template <concepts::LocatableContainer ContainerType> class NeighborSearch : public NeighborSearchBase<ContainerType>
{
  protected:
    using NeighborSearchBase<ContainerType>::_cutoff;
    using NeighborSearchBase<ContainerType>::_system;

    // The number of cells in each direction.
    size_t _ncells_x = 0;
    size_t _ncells_y = 0;
    size_t _ncells_z = 0;

    // The cell list, in the layout a counting sort produces: `_cellitems` holds
    // every binned particle index grouped by cell, and cell `c` owns the slice
    // [_cellstart[c], _cellstart[c + 1]). `_cellstart` therefore has one entry
    // more than there are cells.
    //
    // This used to be an unordered_map from cell id to a vector of indices,
    // which cost a hash lookup and a pointer chase per cell VISITED rather than
    // per cell occupied. That was tolerable while a walk visited 27 cells; it is
    // not once the cells are narrower than the cutoff and a walk visits several
    // hundred (see forcefield/shared/cellgrid_shared.h). The flat layout also
    // puts a cell's particles next to each other in memory, which is what the
    // walk reads them in.
    std::vector<size_t> _cellstart;
    std::vector<size_t> _cellitems;

    // Scratch for the second pass of the counting sort. A member rather than a
    // local so that a rebuild -- which happens every step without a skin --
    // allocates nothing.
    std::vector<size_t> _cellcursor;

    // The particle' bounding box.
    Box _box;

    // Optional particle index to exclude from the cell list. This is used by
    // SpringNetwork for the probe particle, whose interactions are handled
    // explicitly and must not be counted through the regular neighbor grids.
    std::optional<size_t> _excluded_index;

    // Optional list of particle indices to insert in the cell list. When empty,
    // every particle from the system is inserted, except the excluded index.
    // This lets specialised grids stay compact: steric uses all physical
    // particles, electrostatic uses charged particles, and hydrophobic uses
    // hydrophobic particles. Neighbor queries are still evaluated on demand from
    // the current cell list; no particle-to-particle neighbor list is cached.
    std::vector<size_t> _included_indices;

    // Extra distance added to `_cutoff` when sizing grid cells. A positive skin
    // lets `update()` skip rebuilding the grid until a tracked particle has
    // moved far enough that a neighbor could otherwise be missed, trading exact
    // per-call freshness for fewer O(N) rebuilds. Zero (the default) preserves
    // the original behaviour of rebuilding on every `update()` call.
    float _skin = 0.0f;

    // Position of each tracked particle at the last grid rebuild. Only
    // populated and consulted when `_skin > 0`.
    std::vector<std::array<double, 3>> _referencePositions;

    // The cached pair list, in the same layout as the cells: particle `i` owns
    // the slice [_listoffsets[i], _listoffsets[i + 1]) of _listitems.
    //
    // WHY IT EXISTS. Without it every query walks the cells again, so the
    // neighbour search is redone at every step of every term -- measured at 87%
    // of the non-bonded cost, against 13% for the force laws themselves. The
    // list is built once at `cutoff + skin` and reused until something has moved
    // far enough to invalidate it, which is what `_exceeds_skin` decides.
    //
    // It is not an approximation. A pair missing from the list was farther than
    // cutoff + skin when the list was built, so it cannot be inside `cutoff`
    // until the two have closed by `skin`; rebuilding before that can happen
    // makes the set of pairs identical to a fresh search, every step.
    //
    // Only populated when a skin was asked for: with no skin the list would be
    // rebuilt every step and would only add a pass.
    std::vector<size_t> _listoffsets;
    std::vector<uint32_t> _listitems;
    std::vector<bool> _selected;
    bool _listisvalid = false;
    ListPolicy _listpolicy;

    // Width of a cell, in the same unit as the cutoff. NOT the search radius:
    // see forcefield/shared/cellgrid_shared.h for why the two are different
    // things and what narrower cells buy. Zero asks for the historical
    // behaviour, one cell per search radius and a 3x3x3 stencil.
    float _cellwidth = 0.0f;

    // What `_size_grid` settled on, which is `_cellwidth` unless the box needed
    // more than `MAX_CELLS` of it, and the stencil radius that goes with it.
    float _gridwidth = 0.0f;
    int _stencilradius = 1;

    // Ceiling on the number of cells, and so on what `_cellstart` allocates:
    // 8M cells is 64 MB here, and the same cap the OpenCL grid uses.
    static constexpr size_t MAX_CELLS = 8u << 20;

  public:
    // Initializes the neighbor search object with the particles.

    NeighborSearch(const ContainerType & container, float cutoff, float skin = 0.0f, float cellwidth = 0.0f)
        : NeighborSearchBase<ContainerType>(container, cutoff), _skin(skin), _cellwidth(cellwidth)
    {
        _build_grid();
        _build_list();
    }

    NeighborSearch(const ContainerType & container, float cutoff, std::vector<size_t> included_indices,
                    float skin = 0.0f, float cellwidth = 0.0f)
        : NeighborSearchBase<ContainerType>(container, cutoff), _included_indices(std::move(included_indices)),
          _skin(skin), _cellwidth(cellwidth)
    {
        _build_grid();
        _build_list();
    }

    // The cell width the grid actually uses, and the stencil radius that goes
    // with it: how many cells out of its own a particle has to look. Reported
    // for tests and for logging; both are decided by `_size_grid`.
    float cell_width() const { return _gridwidth; }
    int stencil_radius() const { return _stencilradius; }

    // Applies a callback to each neighbor of the given element.
    template <concepts::Locatable T, typename Callback> void for_each_neighbor(const T & element, Callback && callback) const
    {
        // Loops over the particles in the cell of the given particle and in the neighboring cells.
        _for_each_neighbor_cell(concepts::locatable::get_position(element), [&](size_t neighbor_cell_id) {
            const size_t first = _cellstart[neighbor_cell_id];
            const size_t last = _cellstart[neighbor_cell_id + 1];

            for (size_t slot = first; slot < last; slot++)
            {
                const size_t particle_index = _cellitems[slot];
                const T & candidate = _system->at(particle_index);

                // Visit the particle if:
                //   - it is not `element` itself, and
                //   - it is within the cutoff distance.
                if (&candidate != &element && measure::distance(element, candidate) < _cutoff)
                    callback(particle_index);
            }
        });
    }

    // Returns the neighbors of the given element.
    template <concepts::Locatable T> std::vector<size_t> get_neighbors(const T & element) const
    {
        std::vector<size_t> neighbors;
        for_each_neighbor(element, [&](size_t particle_index) { neighbors.push_back(particle_index); });
        return neighbors;
    }

    // Applies a callback to each neighbor of the particle at `index`.
    //
    // This is the query the force loops make, and the only one the cached list
    // can answer: a list is per particle, so the caller has to say which one.
    // The element-based overload above stays for callers that hold a Locatable
    // and no index (Topology, the tests); it always walks the cells.
    template <typename Callback> void for_each_neighbor(size_t index, Callback && callback) const
    {
        if (!_listisvalid)
        {
            for_each_neighbor(_system->at(index), std::forward<Callback>(callback));
            return;
        }

        const auto & element = _system->at(index);
        const size_t first = _listoffsets[index];
        const size_t last = _listoffsets[index + 1];
        for (size_t slot = first; slot < last; slot++)
        {
            const size_t candidate = _listitems[slot];
            // The list was built at cutoff + skin, so it holds pairs that are
            // not neighbours yet. This is the test that makes the answer the
            // same one a fresh search would give.
            if (measure::distance(element, _system->at(candidate)) < _cutoff)
                callback(candidate);
        }
    }

    // Returns the neighbors of the element located at `index` in `_system`.
    std::vector<size_t> get_neighbors(size_t i) const
    {
        std::vector<size_t> neighbors;
        for_each_neighbor(i, [&](size_t particle_index) { neighbors.push_back(particle_index); });
        return neighbors;
    }

    // Rebuilds the cell list based on the system coordinates. When a positive
    // skin was requested, the rebuild is skipped until a tracked particle has
    // drifted far enough since the last rebuild that a neighbor could
    // otherwise be missed (see `_exceeds_skin`).
    void update()
    {
        _listpolicy.step();
        if (!_exceeds_skin())
            return;

        _build_grid();
        if (_listpolicy.worth_it())
            _build_list();
        else
        {
            // Dropped: the queries fall back to walking the cells, which is
            // what they did before there was a list.
            _listisvalid = false;
            _listoffsets.clear();
            _listitems.clear();
        }
        _listpolicy.rebuilt();
    }

    // Whether the list is still judged worth keeping -- see ListPolicy.
    bool list_is_paying() const { return _listpolicy.worth_it(); }

    // How many pairs the cached list holds, and whether there is one at all.
    // For the tests and for reporting; a searcher with no skin has no list.
    bool has_list() const { return _listisvalid; }
    size_t list_size() const { return _listitems.size(); }

    // Excludes one particle index from the grid. The grid is rebuilt
    // unconditionally, bypassing the skin check, so future neighbor queries
    // immediately reflect the exclusion.
    void exclude_index(size_t index)
    {
        _excluded_index = index;
        _build_grid();
        _build_list();
    }

  protected:
    // Returns the search radius used to size grid cells: the physical cutoff
    // used to filter pairs, plus the skin margin. Equal to `_cutoff` when no
    // skin was requested.
    float _search_radius() const { return _cutoff + _skin; }

    // The cell width asked for, which defaults to the search radius so that a
    // caller that says nothing gets the 3x3x3 stencil it always got.
    float _requested_cell_width() const { return _cellwidth > 0.0f ? _cellwidth : _search_radius(); }

    // Returns the total number of cells.
    size_t _number_of_cells() const { return _ncells_x * _ncells_y * _ncells_z; }

    // Returns the cell id of the given position.
    //
    // Clamped at BOTH ends. The grid is sized to the box at the last rebuild
    // and a query is answered from the current position, so with a skin a
    // particle can have drifted below the box minimum -- and the floor of a
    // negative quotient, cast to size_t, is not a small number. It used to be
    // clamped from above only, which sent such a particle to the cell at the
    // far end of the axis instead of the near one.
    size_t _compute_cell(const std::array<double, 3> & position) const
    {
        const auto along = [&](double coordinate, double minimum, size_t ncells) -> size_t {
            const double index = std::floor((coordinate - minimum) / _gridwidth);
            if (index <= 0.0)
                return 0;
            if (index >= static_cast<double>(ncells - 1))
                return ncells - 1;
            return static_cast<size_t>(index);
        };

        const size_t cell_x = along(position[0], _box.min_x(), _ncells_x);
        const size_t cell_y = along(position[1], _box.min_y(), _ncells_y);
        const size_t cell_z = along(position[2], _box.min_z(), _ncells_z);

        // Calculate a unique cell ID for the position
        return cell_x + cell_y * _ncells_x + cell_z * _ncells_x * _ncells_y;
    }

    // Applies a callback to every cell that can hold a neighbour of a particle
    // in `cell_id`: the cube of `_stencilradius` cells around it, minus the
    // corners that the search radius does not reach and minus whatever falls
    // outside the grid.
    //
    // The corner test is what pays for a narrow cell. A stencil is a cube and a
    // cutoff is a sphere, so the wider the stencil the larger the share of it
    // that cannot hold anything -- at 11x11x11 it is a quarter of the cells,
    // dropped for three multiplies each, before a single position is read.
    template <typename Callback> void _for_each_neighbor_cell(size_t cell_id, Callback && callback) const
    {
        const int k = _stencilradius;
        const float radius = _search_radius();
        const float radiussquared = radius * radius;

        const int cx = static_cast<int>(cell_id % _ncells_x);
        const int cy = static_cast<int>((cell_id / _ncells_x) % _ncells_y);
        const int cz = static_cast<int>(cell_id / (_ncells_x * _ncells_y));

        for (int dx = -k; dx <= k; dx++)
        {
            const int x = cx + dx;
            if (x < 0 || x >= static_cast<int>(_ncells_x))
                continue;

            for (int dy = -k; dy <= k; dy++)
            {
                const int y = cy + dy;
                if (y < 0 || y >= static_cast<int>(_ncells_y))
                    continue;

                for (int dz = -k; dz <= k; dz++)
                {
                    const int z = cz + dz;
                    if (z < 0 || z >= static_cast<int>(_ncells_z))
                        continue;

                    if (!biospring_cell_in_range(dx, dy, dz, _gridwidth, radiussquared))
                        continue;

                    callback(static_cast<size_t>(x) + static_cast<size_t>(y) * _ncells_x +
                             static_cast<size_t>(z) * _ncells_x * _ncells_y);
                }
            }
        }
    }

    template <typename Callback>
    void _for_each_neighbor_cell(const std::array<double, 3> & position, Callback && callback) const
    {
        _for_each_neighbor_cell(_compute_cell(position), std::forward<Callback>(callback));
    }

    // Returns the cell ids of the neighboring cells.
    std::vector<size_t> _compute_neighbor_cells(size_t cell_id) const
    {
        std::vector<size_t> neighbor_cell_ids;
        _for_each_neighbor_cell(cell_id, [&](size_t neighbor_cell_id) { neighbor_cell_ids.push_back(neighbor_cell_id); });
        return neighbor_cell_ids;
    }

    // Returns the cells ids of the neighboring cells, given a position.
    std::vector<size_t> _compute_neighbor_cells(const std::array<double, 3> & position) const
    {
        return _compute_neighbor_cells(_compute_cell(position));
    }

    // Applies a callback to the index of every particle this grid bins: all of
    // them, or `_included_indices` when the caller restricted the grid to a
    // subset, minus the excluded one in both cases.
    template <typename Callback> void _for_each_selected(Callback && callback) const
    {
        const auto visit = [&](size_t i) {
            if (_excluded_index && i == *_excluded_index)
                return;
            callback(i);
        };

        if (_included_indices.empty())
        {
            for (size_t i = 0; i < _system->size(); i++)
                visit(i);
        }
        else
        {
            for (size_t i : _included_indices)
                if (i < _system->size())
                    visit(i);
        }
    }

    // Sizes the grid to the box the particles currently occupy, widening the
    // cells if that would need more of them than `MAX_CELLS`.
    //
    // The cap is what makes a flat cell array safe. A hash map only ever held
    // the occupied cells, so an absurd bounding box cost nothing; an array is
    // allocated for the whole box, and a diverging structure produces boxes of
    // any size. Widening rather than refusing keeps the answer exact -- wider
    // cells search more than the cutoff asks for and the distance test drops the
    // surplus -- so the grid degrades towards brute force instead of vanishing.
    void _size_grid()
    {
        _box = measure::box(*_system);

        _gridwidth = _requested_cell_width();
        for (int attempt = 0; attempt < 64; attempt++)
        {
            const auto span = [&](size_t d) { return std::ceil(_box.length()[d] / _gridwidth) + 1.0; };
            const double cells = span(0) * span(1) * span(2);
            if (cells <= static_cast<double>(MAX_CELLS))
            {
                _ncells_x = static_cast<size_t>(span(0));
                _ncells_y = static_cast<size_t>(span(1));
                _ncells_z = static_cast<size_t>(span(2));
                _stencilradius = biospring_stencil_radius(_search_radius(), _gridwidth);
                return;
            }
            _gridwidth *= 2.0f;
        }

        // 64 doublings and still too big means the coordinates are not a
        // structure any more. One cell is degenerate but well defined: every
        // particle lands in it and the search is brute force.
        _ncells_x = _ncells_y = _ncells_z = 1;
        _stencilradius = 1;
    }

    // Builds the cell list, as a counting sort over the cells.
    void _build_grid()
    {
        _size_grid();

        const size_t ncells = _number_of_cells();
        _cellstart.assign(ncells + 1, 0);

        // First pass: how many particles each cell holds, written one slot to
        // the right so that the prefix sum below turns counts into starts.
        _for_each_selected([&](size_t i) { _cellstart[_compute_cell_of(i) + 1]++; });

        for (size_t c = 0; c < ncells; c++)
            _cellstart[c + 1] += _cellstart[c];

        // Second pass: place each particle in its cell's slice.
        _cellitems.resize(_cellstart[ncells]);
        _cellcursor.assign(_cellstart.begin(), _cellstart.end() - 1);
        if (_skin > 0.0f)
            _referencePositions.assign(_system->size(), std::array<double, 3>{});

        _for_each_selected([&](size_t i) {
            const auto & position = concepts::locatable::get_position(_system->at(i));
            _cellitems[_cellcursor[_compute_cell(position)]++] = i;
            if (_skin > 0.0f)
                _referencePositions[i] = position;
        });
    }

    // Collects, for every particle this grid serves, the particles within
    // `cutoff + skin` of it.
    //
    // Skipped entirely without a skin: the list would then be exactly the
    // neighbours of this step, rebuilt next step, so it would cost a pass and
    // a copy of itself to save nothing.
    void _build_list()
    {
        _listisvalid = false;
        if (_skin <= 0.0f)
        {
            _listoffsets.clear();
            _listitems.clear();
            return;
        }

        const size_t n = _system->size();

        // Walked in index order, not in the order `_for_each_selected` hands
        // them over: a caller's `_included_indices` need not be sorted, and the
        // offsets only make sense if they increase.
        _selected.assign(n, false);
        _for_each_selected([&](size_t i) { _selected[i] = true; });

        const float radius = _search_radius();
        _listoffsets.assign(n + 1, 0);
        _listitems.clear();

        for (size_t i = 0; i < n; i++)
        {
            if (_selected[i])
            {
                const auto & element = _system->at(i);
                _for_each_neighbor_cell(concepts::locatable::get_position(element), [&](size_t cell) {
                    for (size_t slot = _cellstart[cell]; slot < _cellstart[cell + 1]; slot++)
                    {
                        const size_t candidate = _cellitems[slot];
                        if (candidate == i)
                            continue;
                        if (measure::distance(element, _system->at(candidate)) < radius)
                            _listitems.push_back(static_cast<uint32_t>(candidate));
                    }
                });
            }
            // Written for every particle, selected or not, so a particle the
            // grid does not serve owns an empty slice rather than a broken one.
            _listoffsets[i + 1] = _listitems.size();
        }

        _listisvalid = true;
    }

    // The cell of the particle at `index`.
    size_t _compute_cell_of(size_t index) const
    {
        return _compute_cell(concepts::locatable::get_position(_system->at(index)));
    }

    // Returns true when the grid must be rebuilt: either no skin was
    // requested (always rebuild), or a tracked particle has moved more than
    // `_skin` since the last rebuild.
    //
    // Correctness note: a query point's own cell is always recomputed from
    // its current position (see `for_each_neighbor`), so only the staleness
    // of *candidate* cell placement threatens correctness. For a true
    // neighbor pair (current distance < `_cutoff`), the distance between a
    // candidate's build-time position and the query's current position is at
    // most `candidate_drift + _cutoff`. Cells built with radius
    // `_cutoff + _skin` are guaranteed to still find that candidate as long
    // as `candidate_drift <= _skin`, which is exactly what this check
    // enforces (per particle, not per pair, so it is safe for every query).
    bool _exceeds_skin() const
    {
        if (_skin <= 0.0f)
            return true;

        const auto has_drifted = [&](size_t i) {
            if (_excluded_index && i == *_excluded_index)
                return false;

            const auto & current = concepts::locatable::get_position(_system->at(i));
            const auto & reference = _referencePositions[i];
            const double dx = current[0] - reference[0];
            const double dy = current[1] - reference[1];
            const double dz = current[2] - reference[2];
            // HALF the skin. The cell grid alone could afford the whole of it,
            // because a query recomputes its own cell from its current position
            // and only the candidate was stale. A pair list has no fresh end:
            // both particles were placed when it was built, so two of them
            // drifting by skin/2 towards each other close the whole skin.
            return std::sqrt(dx * dx + dy * dy + dz * dz) > 0.5 * static_cast<double>(_skin);
        };

        if (_included_indices.empty())
        {
            for (size_t i = 0; i < _system->size(); i++)
                if (has_drifted(i))
                    return true;
        }
        else
        {
            for (size_t i : _included_indices)
                if (i < _system->size() && has_drifted(i))
                    return true;
        }
        return false;
    }
};

} // namespace nsearch
} // namespace biospring

#endif // __NSEARCH_HPP__