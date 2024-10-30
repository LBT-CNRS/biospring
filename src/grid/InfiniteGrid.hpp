#ifndef __INFINITE_GRID_HPP__
#define __INFINITE_GRID_HPP__

#include <array>
#include <unordered_map>
#include <vector>

#include "GridCoordinatesSystem.hpp"

namespace biospring
{
namespace grid
{

template <typename T> class InfiniteGrid
{
  protected:
    std::unordered_map<discrete_coordinates, T> _data;
    InfiniteGridCoordinatesSystem _coordinates;

  public:
    InfiniteGrid() {}
    InfiniteGrid(const std::array<double, 3> & cell_size) : _coordinates(cell_size) {}

    void initialize(const std::array<double, 3> & cell_size)
    {
        if (!_data.empty())
            throw std::runtime_error("Trying to initialize grid which is not empty.");
        _coordinates.initialize(cell_size);
    }

    const std::array<double, 3> & cell_size() const { return _coordinates.cell_size(); }

    // Returns the number of cells in the grid.
    size_t size() const { return _data.size(); }

    // Remove data from the grid.
    void clear() { _data.clear(); }

    // Returns the cell coordinates of the given position.
    discrete_coordinates cell_coordinates(const real_coordinates & position) const
    {
        return _coordinates.cell_coordinates(position);
    }

    // Returns all the cells that can contain an element within a given radius.
    std::vector<discrete_coordinates> cells_within_radius(const discrete_coordinates & cell, double radius) const
    {
        return _coordinates.cells_within_radius(cell, radius);
    }

    // Returns true if a cell with given coordinates exists.
    bool has_cell(const discrete_coordinates & cell) const { return _data.contains(cell); }

    // Returns the element at the given cell coordinates.
    T & at(const discrete_coordinates & cell) { return _data[cell]; }
    T & operator()(const discrete_coordinates & cell) { return _data[cell]; }

    // Returns the element at the given position.
    T & at(const real_coordinates & position) { return at(_coordinates.cell_coordinates(position)); }
    T & operator()(const real_coordinates & cell) { return _data[cell]; }

    // Adds an element to the grid.
    void add(const discrete_coordinates & cell, const T & value) { _data[cell] = value; }
    void add(const real_coordinates & position, const T & value)
    {
        add(_coordinates.cell_coordinates(position), value);
    }
};

// Specialization for when T is std::vector
template <typename ElementType> class InfiniteGridOfContainers : public InfiniteGrid<std::vector<ElementType>>
{
  protected:
    using InfiniteGrid<std::vector<ElementType>>::_data;
    using InfiniteGrid<std::vector<ElementType>>::_coordinates;

  public:
    void add(const discrete_coordinates & cell, const ElementType & element) { _data[cell].push_back(element); }
    void add(const real_coordinates & position, const ElementType & element)
    {
        add(_coordinates.cell_coordinates(position), element);
    }
    void add(const std::array<double, 3> & position, const ElementType & element)
    {
        add(_coordinates.cell_coordinates(position), element);
    }

    // Returns all the cells that can contain an element within a given radius.
    std::vector<discrete_coordinates> cells_within_radius(const discrete_coordinates & cell, double radius) const
    {
        return _coordinates.cells_within_radius(cell, radius);
    }

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

#endif // __INFINITE_GRID_HPP__