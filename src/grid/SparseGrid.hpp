#ifndef __SPARSE_GRID_HPP__
#define __SPARSE_GRID_HPP__

#include <array>
#include <unordered_map>
#include <vector>

#include "GridCoordinatesSystem.hpp"

namespace biospring
{
namespace grid
{

template <typename T> class SparseGrid
{
  protected:
    // Grid data.
    std::unordered_map<discrete_coordinates, T> _data;

    // Grid coordinates system.
    GridCoordinatesSystem _coordinates_system;

  public:
    // Initializes the grid from a box and a cell size.
    SparseGrid(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
        : _coordinates_system(box, cell_size)
    {
    }

    // Initializes the grid from an origin, a shape and a cell size.
    SparseGrid(const std::array<double, 3> & origin, const std::array<size_t, 3> & shape,
               const std::array<double, 3> & cell_size)
        : _coordinates_system(origin, shape, cell_size)
    {
    }

    SparseGrid() : _coordinates_system() {}

    // `reshape` function: allows to use the empty constructor.
    void reshape(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
    {
        _coordinates_system.initialize(box, cell_size);
    }

    // == Shortcuts to coordinates system methods =================================

    const std::array<size_t, 3> & shape() const { return _coordinates_system.shape(); }
    const std::array<double, 3> & cell_size() const { return _coordinates_system.cell_size(); }
    const Box & boundaries() const { return _coordinates_system.boundaries(); }
    std::array<double, 3> origin() const { return _coordinates_system.boundaries().origin(); }

    bool is_out_of_grid(const discrete_coordinates & cell) const { return _coordinates_system.is_out_of_grid(cell); }
    bool is_out_of_grid(const real_coordinates & position) const
    {
        return _coordinates_system.is_out_of_grid(position);
    }

    // ============================================================================

    // Returns the number of cells in the grid.
    size_t size() const { return _data.size(); }

    // Returns the maximum number of cells in the grid.
    size_t max_size() const { return _coordinates_system.max_size(); }

    // Empties the grid.
    void clear() { _data.clear(); }

    // Returns the cell coordinates of the given position.
    discrete_coordinates cell_coordinates(const real_coordinates & position) const
    {
        return _coordinates_system.cell_coordinates(position);
    }

    // Returns true if a cell with given coordinates exists.
    bool has_cell(const discrete_coordinates & cell) const { return _data.contains(cell); }

    // Returns the element at the given cell coordinates.
    // Throws an exception if the cell is out of grid boundaries.
    T & at(const discrete_coordinates & cell)
    {
        // Checks if the cell is out of grid boundaries
        if (_coordinates_system.is_out_of_grid(cell))
            throw std::out_of_range("Cell is out of grid boundaries.");

        // Returns the element at the given cell coordinates
        return _data[cell];
    }

    // Returns the element at the given cell coordinates.
    T & at(const real_coordinates & position) { return at(_coordinates_system.cell_coordinates(position)); }

    // Adds an element to the grid.
    void add(const discrete_coordinates & cell, const T & value)
    {
        if (_coordinates_system.is_out_of_grid(cell))
            throw std::out_of_range("Cell is out of grid boundaries.");
        _data[cell] = value;
    }

    void add(const real_coordinates & position, const T & value) { add(_coordinates_system.at(position), value); }

    // Returns the cell indices of the cells within a given offset from a given cell.
    std::vector<discrete_coordinates> cells_within_offset(const discrete_coordinates & cell,
                                                          const discrete_coordinates & offset) const
    {
        return _coordinates_system.cells_within_offset(cell, offset);
    }
};

template <typename T> class SparseGridOfContainers : public SparseGrid<std::vector<T>>
{
  protected:
    using SparseGrid<std::vector<T>>::_data;
    using SparseGrid<std::vector<T>>::_coordinates_system;

  public:
    using SparseGrid<std::vector<T>>::SparseGrid;
    using SparseGrid<std::vector<T>>::reshape;
    using SparseGrid<std::vector<T>>::at;
    using SparseGrid<std::vector<T>>::is_out_of_grid;

    SparseGridOfContainers(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
        : SparseGrid<std::vector<T>>(box, cell_size)
    {
    }

    SparseGridOfContainers() : SparseGrid<std::vector<T>>() {}

    // Adds an element to the grid.
    void add(const discrete_coordinates & cell, const T & value) { _data[cell].push_back(value); }
    void add(const real_coordinates & position, const T & value) { add(_coordinates_system.at(position), value); }

    // Returns the number of elements in the grid.
    size_t number_of_elements() const
    {
        size_t size = 0;
        for (const auto & v : _data)
            for (const auto & w : v)
                for (const auto & x : w)
                    size += x.size();
        return size;
    }
};

} // namespace grid
} // namespace biospring

#endif // __SPARSE_GRID_HPP__