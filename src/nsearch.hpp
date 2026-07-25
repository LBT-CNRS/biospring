// Neighbor search.
//
// This is a simple implementation of neighbor search. It is not
// optimized for performance, but it is easy to use and understand.
//
// The neighbor search is performed in two steps. First, the
// `NeighborSearch` class is initialized with the particles.
// Then, the `find_neighbors` method is called to find the
// neighbors of a given particle. The neighbors are returned as a
// vector of indices.
//
// Example:
//
//     // Create a neighbor search object.
//     NeighborSearch nsearch;
//
//     // Initialize the neighbor search object with the particles.
//     nsearch.init(particles);
//
//     // Find the neighbors of particle 0.
//     auto neighbors = nsearch.find_neighbors(0);
//
//     // Print the neighbors.
//     for (auto i : neighbors) {
//         std::cout << i << std::endl;
//     }
//

#ifndef __NSEARCH_HPP__
#define __NSEARCH_HPP__

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

#include "box.hpp"
#include "concepts.hpp"
#include "grid/InfiniteGrid.hpp"
#include "measure.hpp"

namespace biospring
{
namespace nsearch
{

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

template <concepts::LocatableContainer ContainerType> class NeighborSearch2 : public NeighborSearchBase<ContainerType>
{
  protected:
    using NeighborSearchBase<ContainerType>::_cutoff;
    using NeighborSearchBase<ContainerType>::_system;

    grid::InfiniteGridOfContainers<size_t> _grid;

  public:
    NeighborSearch2(const ContainerType & container, float cutoff)
        : NeighborSearchBase<ContainerType>(container, cutoff)
    {
        _populate();
    }

    template <concepts::Locatable T> std::vector<size_t> get_neighbors(const T & element) const
    {
        std::vector<size_t> neighbors;
        auto cells = _grid.cells_within_radius(concepts::locatable::get_position(element), _cutoff);

        for (auto cell : cells)
        {
            // This is an infinite grid, so we need to check if the cell exists.
            if (!_grid.has_cell(cell))
                continue;

            for (auto index : _grid.at(cell))
            {
                const T & candidate = _system->at(index);

                if (&candidate != &element && measure::distance(element, candidate) < _cutoff)
                    neighbors.push_back(index);
            }
        }
        return neighbors;
    }

  protected:
    void _populate()
    {
        _grid.clear();
        _grid.initialize({_cutoff, _cutoff, _cutoff});

        for (size_t i = 0; i < _system->size(); i++)
        {
            const auto & position = concepts::locatable::get_position(_system->at(i));
            _grid.add(position, i);
        }
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

    // The cell list, which maps each cell to the particles it contains.
    std::unordered_map<size_t, std::vector<size_t>> _cells;

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

  public:
    // Initializes the neighbor search object with the particles.

    NeighborSearch(const ContainerType & container, float cutoff, float skin = 0.0f)
        : NeighborSearchBase<ContainerType>(container, cutoff), _skin(skin)
    {
        _build_grid();
    }

    NeighborSearch(const ContainerType & container, float cutoff, std::vector<size_t> included_indices,
                    float skin = 0.0f)
        : NeighborSearchBase<ContainerType>(container, cutoff), _included_indices(std::move(included_indices)),
          _skin(skin)
    {
        _build_grid();
    }

    // Applies a callback to each neighbor of the given element.
    template <concepts::Locatable T, typename Callback> void for_each_neighbor(const T & element, Callback && callback) const
    {
        // Loops over the particles in the cell of the given particle and in the neighboring cells.
        _for_each_neighbor_cell(concepts::locatable::get_position(element), [&](size_t neighbor_cell_id) {
            // If cell does not exist (aka is empty), skip it.
            const auto cell = _cells.find(neighbor_cell_id);
            if (cell == _cells.end())
                return;

            for (size_t particle_index : cell->second)
            {
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

    // Returns the neighbors of the element located at `index` in `_system`.
    std::vector<size_t> get_neighbors(size_t i) const { return get_neighbors(_system->at(i)); }

    // Rebuilds the cell list based on the system coordinates. When a positive
    // skin was requested, the rebuild is skipped until a tracked particle has
    // drifted far enough since the last rebuild that a neighbor could
    // otherwise be missed (see `_exceeds_skin`).
    void update()
    {
        if (!_exceeds_skin())
            return;
        _build_grid();
    }

    // Excludes one particle index from the grid. The grid is rebuilt
    // unconditionally, bypassing the skin check, so future neighbor queries
    // immediately reflect the exclusion.
    void exclude_index(size_t index)
    {
        _excluded_index = index;
        _build_grid();
    }

  protected:
    // Returns the search radius used to size grid cells: the physical cutoff
    // used to filter pairs, plus the skin margin. Equal to `_cutoff` when no
    // skin was requested.
    float _search_radius() const { return _cutoff + _skin; }

    // Returns the total number of cells.
    size_t _number_of_cells() const { return _ncells_x * _ncells_y * _ncells_z; }

    // Returns the cell id of the given position.
    size_t _compute_cell(const std::array<double, 3> & position) const
    {
        const float radius = _search_radius();

        // Calculate grid cell coordinates for the given position, considering negative coordinates
        size_t cell_y = static_cast<size_t>((position[1] - _box.min_y()) / radius);
        size_t cell_z = static_cast<size_t>((position[2] - _box.min_z()) / radius);
        size_t cell_x = static_cast<size_t>((position[0] - _box.min_x()) / radius);

        // Ensure that the cell coordinates are within bounds
        cell_x = std::min(cell_x, _ncells_x - 1);
        cell_y = std::min(cell_y, _ncells_y - 1);
        cell_z = std::min(cell_z, _ncells_z - 1);

        // Calculate a unique cell ID for the position
        return cell_x + cell_y * _ncells_x + cell_z * _ncells_x * _ncells_y;
    }

    template <typename Callback> void _for_each_neighbor_cell(size_t cell_id, Callback && callback) const
    {
        for (int dx = -1; dx <= 1; dx++)
        {
            for (int dy = -1; dy <= 1; dy++)
            {
                for (int dz = -1; dz <= 1; dz++)
                {
                    int neighbor_cell_x = static_cast<int>(cell_id % _ncells_x) + dx;
                    int neighbor_cell_y = static_cast<int>((cell_id / _ncells_x) % _ncells_y) + dy;
                    int neighbor_cell_z = static_cast<int>(cell_id / (_ncells_x * _ncells_y)) + dz;

                    if (neighbor_cell_x < 0 || neighbor_cell_x >= static_cast<int>(_ncells_x))
                        continue;
                    if (neighbor_cell_y < 0 || neighbor_cell_y >= static_cast<int>(_ncells_y))
                        continue;
                    if (neighbor_cell_z < 0 || neighbor_cell_z >= static_cast<int>(_ncells_z))
                        continue;

                    callback(static_cast<size_t>(neighbor_cell_x) + static_cast<size_t>(neighbor_cell_y) * _ncells_x +
                             static_cast<size_t>(neighbor_cell_z) * _ncells_x * _ncells_y);
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

    // Builds the cell list.
    void _build_grid()
    {
        // Clears the cell list.
        _cells.clear();

        // Computes the simulation box.
        _box = measure::box(*_system);

        // Computes the number of cells in each direction.
        const float radius = _search_radius();
        _ncells_x = size_t(std::ceil(_box.length()[0] / radius) + 1);
        _ncells_y = size_t(std::ceil(_box.length()[1] / radius) + 1);
        _ncells_z = size_t(std::ceil(_box.length()[2] / radius) + 1);

        if (_skin > 0.0f)
            _referencePositions.assign(_system->size(), std::array<double, 3>{});

        const auto add_particle_to_cell = [&](size_t i) {
            if (_excluded_index && i == *_excluded_index)
                return;

            const auto & position = concepts::locatable::get_position(_system->at(i));

            // Computes the cell number of the particle.
            size_t cell_id = _compute_cell(position);

            // Adds the particle to the appropriate cell.
            _cells[cell_id].push_back(i);

            if (_skin > 0.0f)
                _referencePositions[i] = position;
        };

        // Finds the cell of each selected particle.
        if (_included_indices.empty())
        {
            for (size_t i = 0; i < _system->size(); i++)
                add_particle_to_cell(i);
        }
        else
        {
            for (size_t i : _included_indices)
            {
                if (i >= _system->size())
                    continue;
                add_particle_to_cell(i);
            }
        }
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
            return std::sqrt(dx * dx + dy * dy + dz * dz) > static_cast<double>(_skin);
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

// Dynamic neighbor search.
//
// This class is aimed to be used when computing the neighbor list of elements
// multiple times.
//
// Implementation details:
//   Compared to `NeighborSearch`, this class stores the indices of the neighbor cells
//   of each cell.
template <concepts::LocatableContainer ContainerType> class NeighborSearchDynamic : public NeighborSearch<ContainerType>
{
  protected:
    using NeighborSearch<ContainerType>::_cutoff;
    using NeighborSearch<ContainerType>::_system;
    using NeighborSearch<ContainerType>::_cells;

    // The neighbor cells of each cell.
    std::unordered_map<size_t, std::vector<size_t>> _neighbor_cells;

  public:
    // using NeighborSearch<ContainerType>::NeighborSearch;

    NeighborSearchDynamic(const ContainerType & container, float cutoff)
        : NeighborSearch<ContainerType>(container, cutoff)
    {
        _build_grid();
    }

    NeighborSearchDynamic(const ContainerType & container, float cutoff, std::vector<size_t> included_indices)
        : NeighborSearch<ContainerType>(container, cutoff, std::move(included_indices))
    {
        _build_grid();
    }

    // Returns the neighbors of the given element.
    template <concepts::Locatable T> std::vector<size_t> get_neighbors(const T & element) const
    {
        return NeighborSearch<ContainerType>::get_neighbors(element);
    }

    // Rebuilds the cell list and the cached neighbor-cell topology.
    void update() { _build_grid(); }

    // Returns the neighbors of the element located at `index` in `_system`.
    std::vector<size_t> get_neighbors(size_t index) const
    {
        size_t cell_id = _compute_cell(concepts::locatable::get_position(_system->at(index)));

        std::vector<size_t> neighbors;
        for (auto neighbor_cell_id : _neighbor_cells.at(cell_id))
        {
            const auto cell = _cells.find(neighbor_cell_id);
            if (cell == _cells.end())
                continue;

            for (size_t particle_index : cell->second)
            {
                const auto & candidate = _system->at(particle_index);

                if (&candidate != &_system->at(index) && measure::distance(_system->at(index), candidate) < _cutoff)
                    neighbors.push_back(particle_index);
            }
        }
        return neighbors;
    }

  protected:
    using NeighborSearch<ContainerType>::_build_grid;
    using NeighborSearch<ContainerType>::_compute_cell;
    using NeighborSearch<ContainerType>::_compute_neighbor_cells;
    using NeighborSearch<ContainerType>::_number_of_cells;

    // Constructs the cell list and the neighbor cells of each cell.
    void _build_grid()
    {
        NeighborSearch<ContainerType>::_build_grid();

        // Computes the neighbor cells of each cell.
        _neighbor_cells.clear();
        for (size_t cell_id = 0; cell_id < _number_of_cells(); cell_id++)
            _neighbor_cells[cell_id] = _compute_neighbor_cells(cell_id);
    }
};

} // namespace nsearch
} // namespace biospring

#endif // __NSEARCH_HPP__