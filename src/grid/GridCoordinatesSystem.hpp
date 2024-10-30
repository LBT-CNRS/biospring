#ifndef __GRID_COORDINATES_SYSTEM_HPP__
#define __GRID_COORDINATES_SYSTEM_HPP__

#include <iostream>

#include <array>
#include <initializer_list>
#include <vector>

#include "measure.hpp"

namespace biospring
{
namespace grid
{

template <typename T> struct _coordinates;
using discrete_coordinates = _coordinates<int>;
using real_coordinates = _coordinates<double>;

template <typename T> struct _coordinates
{
    T x, y, z;

    // Constructs from 3 numbers.
    _coordinates(T x, T y, T z) : x(x), y(y), z(z) {}

    // Constructs from any container with 3 elements.
    template <typename Container> _coordinates(const Container & container)
    {
        if (container.size() != 3)
            throw std::invalid_argument("Container must have 3 elements.");
        x = container[0];
        y = container[1];
        z = container[2];
    }

    // Constructs from initializer list.
    _coordinates(std::initializer_list<T> list)
    {
        if (list.size() != 3)
            throw std::invalid_argument("Initializer list must have 3 elements.");
        x = *(list.begin());
        y = *(list.begin() + 1);
        z = *(list.begin() + 2);
    }

    // Provides the `size()` method necessary for container initialization.
    size_t size() const { return 3; }

    // Provides subscript operator.
    T operator[](size_t i) const
    {
        if (i >= 3)
            throw std::out_of_range("Index out of range.");
        return i == 0 ? x : i == 1 ? y : z;
    }

    // Equality operator.
    bool operator==(const _coordinates & other) const { return x == other.x && y == other.y && z == other.z; }
};

// Calculates 3D grid coordinates.
// No boundary checks are performed in this class.
class InfiniteGridCoordinatesSystem
{
  protected:
    std::array<double, 3> _cell_size;
    std::array<double, 3> _origin;

  public:
    InfiniteGridCoordinatesSystem() { _origin = {0.0, 0.0, 0.0}; }
    InfiniteGridCoordinatesSystem(const std::array<double, 3> & cell_size)
        : _cell_size(cell_size), _origin({0.0, 0.0, 0.0})
    {
    }

    const std::array<double, 3> & cell_size() const { return _cell_size; }
    const std::array<double, 3> & origin() const { return _origin; }

    // Initializes the grid coordinates system from a cell size.
    void initialize(const std::array<double, 3> & cell_size) { _cell_size = cell_size; }

    // Returns the cell coordinates for the given position.
    discrete_coordinates cell_coordinates(const real_coordinates & position) const
    {
        // Calculates grid cell coordinates for the given position, considering negative coordinates
        int cell_x = static_cast<int>((position.x - origin()[0]) / cell_size()[0]);
        int cell_y = static_cast<int>((position.y - origin()[1]) / cell_size()[1]);
        int cell_z = static_cast<int>((position.z - origin()[2]) / cell_size()[2]);

        // Returns the cell coordinates for the position
        return {cell_x, cell_y, cell_z};
    }

    // Returns the cells that are within a given offset of a given cell.
    // The current cell is returned as well, for API convenience reasons.
    std::vector<discrete_coordinates> cells_within_offset(const discrete_coordinates & cell,
                                                          const discrete_coordinates & offset) const
    {
        std::vector<discrete_coordinates> cells;

        for (int i = -offset.x; i <= offset.x; i++)
            for (int j = -offset.y; j <= offset.y; j++)
                for (int k = -offset.z; k <= offset.z; k++)
                    cells.push_back(discrete_coordinates(cell.x + i, cell.y + j, cell.z + k));
        return cells;
    }

    // Returns the coordinates of the cells in a given radius around a given cell.
    // Returns all cells where elements within the given radius could be found.
    // `cell` is returned as well, for API convenience reasons.
    std::vector<discrete_coordinates> cells_within_radius(const discrete_coordinates & cell, double radius) const
    {
        // Number of cells in each direction.
        int n_x = static_cast<int>(std::ceil(radius / _cell_size[0]));
        int n_y = static_cast<int>(std::ceil(radius / _cell_size[1]));
        int n_z = static_cast<int>(std::ceil(radius / _cell_size[2]));

        // List of cells within the given radius.
        std::vector<discrete_coordinates> cells;

        // Loops over the cells within the given radius.
        for (int i = -n_x; i <= n_x; i++)
            for (int j = -n_y; j <= n_y; j++)
                for (int k = -n_z; k <= n_z; k++)
                {
                    // Computes the cell coordinates.
                    discrete_coordinates neighbor_cell{cell.x + i, cell.y + j, cell.z + k};

                    // Adds the cell to the list.
                    cells.push_back(neighbor_cell);
                }
        return cells;
    }

    // Returns the coordinates of the cells in a given radius around a given position.
    // Returns all cells where elements within the given radius could be found.
    // The cell of `position` is returned as well, for API convenience reasons.
    std::vector<discrete_coordinates> cells_within_radius(const real_coordinates & position, double radius) const
    {
        return cells_within_radius(cell_coordinates(position), radius);
    }
};

// Calculates 3D grid coordinates.
class GridCoordinatesSystem : public InfiniteGridCoordinatesSystem
{
  protected:
    std::array<size_t, 3> _shape;
    Box _boundaries;
    bool _initialized = false;

  public:
    // Initializes the grid coordinates system from a box and a cell size.
    GridCoordinatesSystem(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
        : InfiniteGridCoordinatesSystem(cell_size), _boundaries(box), _initialized(true)
    {
        // Computes the grid dimensions.
        _shape[0] = static_cast<size_t>(std::ceil(_boundaries.length()[0] / _cell_size[0]));
        _shape[1] = static_cast<size_t>(std::ceil(_boundaries.length()[1] / _cell_size[1]));
        _shape[2] = static_cast<size_t>(std::ceil(_boundaries.length()[2] / _cell_size[2]));
        _origin = _boundaries.origin();
    }

    // Initializes the grid coordinates system from an origin, a shape and a cell size.
    GridCoordinatesSystem(const std::array<double, 3> & origin, const std::array<size_t, 3> & shape,
                          const std::array<double, 3> & cell_size)
        : InfiniteGridCoordinatesSystem(cell_size), _shape(shape), _initialized(true)
    {
        // Initializes the grid boundaries using origin and shape.
        _boundaries = Box(origin, {shape[0] * cell_size[0], shape[1] * cell_size[1], shape[2] * cell_size[2]});
        _origin = _boundaries.origin();
    }

    GridCoordinatesSystem() {}

    // Initialize function: allows to use the empty constructor.
    void initialize(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
    {
        InfiniteGridCoordinatesSystem::initialize(cell_size);
        _boundaries = Box(box);

        // Computes the grid dimensions.
        _shape[0] = static_cast<size_t>(std::ceil(_boundaries.length()[0] / _cell_size[0]));
        _shape[1] = static_cast<size_t>(std::ceil(_boundaries.length()[1] / _cell_size[1]));
        _shape[2] = static_cast<size_t>(std::ceil(_boundaries.length()[2] / _cell_size[2]));

        _origin = _boundaries.origin();
        _initialized = true;
    }

    // Getters.
    const std::array<size_t, 3> & shape() const { return _shape; }
    const Box & boundaries() const { return _boundaries; }

    // Returns the maximum number of cells in the system.
    size_t max_size() const
    {
        if (!_initialized)
            return 0;
        return _shape[0] * _shape[1] * _shape[2];
    }

    // Returns wheters a position is outside the grid boundaries.
    bool is_out_of_grid(const real_coordinates & position) const
    {
        double min_x = _boundaries.min_x();
        double min_y = _boundaries.min_y();
        double min_z = _boundaries.min_z();

        double max_x = _boundaries.max_x() - 1e-6;
        double max_y = _boundaries.max_y() - 1e-6;
        double max_z = _boundaries.max_z() - 1e-6;

        return position.x < min_x || position.x > max_x || position.y < min_y || position.y > max_y ||
               position.z < min_z || position.z > max_z;
    }

    // Returns wheters a cell in outside the grid boundaries.
    bool is_out_of_grid(const discrete_coordinates & cell) const
    {
        return cell.x < 0 || cell.x >= _shape[0] || cell.y < 0 || cell.y >= _shape[1] || cell.z < 0 ||
               cell.z >= _shape[2];
    }

    // =============================================================================
    // Safe cell coordinates retrieval methods.
    //
    // Boundary checks are performed in these methods.
    // =============================================================================

    // Returns the cell coordinates for the given position.
    discrete_coordinates cell_coordinates(const real_coordinates & position) const
    {
        if (is_out_of_grid(position))
            throw std::out_of_range("Position is out of grid boundaries.");
        return InfiniteGridCoordinatesSystem::cell_coordinates(position);
    }

    // Returns the cell coordinates for the given cell.
    discrete_coordinates at(const discrete_coordinates & cell)
    {
        if (is_out_of_grid(cell))
            throw std::out_of_range("Cell is out of grid boundaries.");
        return cell;
    }

    // Returns the cell coordinates for the given position.
    discrete_coordinates at(const real_coordinates & position)
    {
        if (is_out_of_grid(position))
            throw std::out_of_range("Position is out of grid boundaries.");
        return cell_coordinates(position);
    }

    // Returns the cells that are within a given offset of a given cell.
    // The current cell is returned as well, for API convenience reasons.
    // Does not return cells that are out of grid boundaries.
    std::vector<discrete_coordinates> cells_within_offset(const discrete_coordinates & cell,
                                                          const discrete_coordinates & offset) const
    {
        auto cells = InfiniteGridCoordinatesSystem::cells_within_offset(cell, offset);
        std::erase_if(cells, [&](const discrete_coordinates & cell) { return is_out_of_grid(cell); });
        return cells;
    }

    // Returns the coordinates of the cells in a given radius around a given position.
    // The current cell is returned as well, for API convenience reasons.
    std::vector<discrete_coordinates> cells_within_radius(const discrete_coordinates & cell, double radius) const
    {
        auto cells = InfiniteGridCoordinatesSystem::cells_within_radius(cell, radius);
        std::erase_if(cells, [&](const discrete_coordinates & cell) { return is_out_of_grid(cell); });
        return cells;
    }

    // Iteration over the grid cell discretized coordinates.
    class iterator
    {
      protected:
        const GridCoordinatesSystem & _grid;
        discrete_coordinates _cell;

      public:
        iterator(const GridCoordinatesSystem & grid, const discrete_coordinates & cell) : _grid(grid), _cell(cell) {}

        iterator & operator++()
        {
            // Increment to the next cell
            _cell.z++;
            if (_cell.z >= _grid._shape[2])
            {
                _cell.z = 0;
                _cell.y++;
                if (_cell.y >= _grid._shape[1])
                {
                    _cell.y = 0;
                    _cell.x++;
                }
            }
            return *this;
        }

        bool operator!=(const iterator & other) const { return _cell != other._cell; }

        discrete_coordinates operator*() const { return _cell; }
    };

    iterator begin() { return iterator(*this, {0, 0, 0}); }
    iterator end() { return iterator(*this, {static_cast<int>(_shape[0]), 0, 0}); }
};

} // namespace grid
} // namespace biospring

// Defines hash function for discrete_coordinates.
namespace std
{

// Define hash_combine function
template <class T> void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <> struct hash<biospring::grid::discrete_coordinates>
{
    size_t operator()(const biospring::grid::discrete_coordinates & t) const
    {
        size_t hashVal = 0;
        hash_combine(hashVal, t.x);
        hash_combine(hashVal, t.y);
        hash_combine(hashVal, t.z);
        return hashVal;
    }
};
} // namespace std

#endif // __GRID_COORDINATES_SYSTEM_HPP__