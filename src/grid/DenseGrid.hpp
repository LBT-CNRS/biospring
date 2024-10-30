#ifndef __DENSE_GRID_HPP__
#define __DENSE_GRID_HPP__

#include <array>
#include <unordered_map>
#include <vector>

#include "GridCoordinatesSystem.hpp"

namespace biospring
{
namespace grid
{

// A DenseGrid is a 3D grid with regular cells.
//
// The grid is stored as a 3D array of elements of type T.
// Therefore, it is optimized for fast access to the elements but can be memory consuming.
template <typename T> class DenseGrid
{
  protected:
    // Grid data.
    std::vector<std::vector<std::vector<T>>> _data;

    // Grid coordinates system.
    GridCoordinatesSystem _coordinates_system;

  public:
    // Initializes the grid from a box and a cell size.
    DenseGrid(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
        : _coordinates_system(box, cell_size)
    {
        _reshape();
    }

    // Initializes the grid from an origin, a shape and a cell size.
    DenseGrid(const std::array<double, 3> & origin, const std::array<size_t, 3> & shape,
              const std::array<double, 3> & cell_size)
        : _coordinates_system(origin, shape, cell_size)
    {
        _reshape();
    }

    DenseGrid() : _coordinates_system() {}

    // `reshape` function: allows to use the empty constructor.
    void reshape(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
    {
        _coordinates_system.initialize(box, cell_size);
        _reshape();
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

    // Returns true if a cell with given coordinates exists.
    bool has_cell(const discrete_coordinates & cell) const { return !is_out_of_grid(cell); }

    // ============================================================================

    // Returns the number of cells in the grid.
    size_t size() const { return _coordinates_system.max_size(); }

    // Empties the grid.
    void clear()
    {
        for (auto & v : _data)
            for (auto & w : v)
                w.clear();
        _reshape();
    }

    // Returns the cell coordinates of the given position.
    discrete_coordinates cell_coordinates(const real_coordinates & position) const
    {
        return _coordinates_system.cell_coordinates(position);
    }

    // Returns the element at the given cell coordinates.
    // Throws an exception if the cell is out of grid boundaries.
    T & at(const discrete_coordinates & cell)
    {
        if (_coordinates_system.is_out_of_grid(cell))
            throw std::out_of_range("Cell is out of grid boundaries.");
        return _data[cell[0]][cell[1]][cell[2]];
    }

    const T & at(const discrete_coordinates & cell) const
    {
        if (_coordinates_system.is_out_of_grid(cell))
            throw std::out_of_range("Cell is out of grid boundaries.");
        return _data[cell[0]][cell[1]][cell[2]];
    }

    // Returns the element at the given position.
    T & at(const real_coordinates & position) { return at(_coordinates_system.cell_coordinates(position)); }
    const T & at(const real_coordinates & position) const { return at(_coordinates_system.cell_coordinates(position)); }
    T & get(const float x, const float y, const float z) { return at(real_coordinates(x, y, z)); }
    const T & get(const float x, const float y, const float z) const { return at(real_coordinates(x, y, z)); }

    // Returns the cell indices of the cells within a given offset from a given cell.
    std::vector<discrete_coordinates> cells_within_offset(const discrete_coordinates & cell,
                                                          const discrete_coordinates & offset) const
    {
        return _coordinates_system.cells_within_offset(cell, offset);
    }

  protected:
    // Resizes the grid data.
    void _reshape()
    {
        _data.resize(_coordinates_system.shape()[0]);
        for (auto & v : _data)
        {
            v.resize(_coordinates_system.shape()[1]);
            for (auto & w : v)
                w.resize(_coordinates_system.shape()[2]);
        }
    }
};

template <typename T> class DenseGridOfContainers : public DenseGrid<std::vector<T>>
{
  protected:
    using DenseGrid<std::vector<T>>::_data;
    using DenseGrid<std::vector<T>>::_coordinates_system;

  public:
    using DenseGrid<std::vector<T>>::DenseGrid;
    using DenseGrid<std::vector<T>>::reshape;
    using DenseGrid<std::vector<T>>::at;
    using DenseGrid<std::vector<T>>::is_out_of_grid;

    DenseGridOfContainers(const std::array<double, 6> & box, const std::array<double, 3> & cell_size)
        : DenseGrid<std::vector<T>>(box, cell_size)
    {
    }

    DenseGridOfContainers() : DenseGrid<std::vector<T>>() {}

    // Adds an element to the grid.
    void add(const discrete_coordinates & cell, const T & value) { at(cell).push_back(value); }
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

#endif // __DENSE_GRID_HPP__