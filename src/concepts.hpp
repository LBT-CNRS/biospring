#ifndef __CONCEPTS_HPP__
#define __CONCEPTS_HPP__

#include <array>
#include <concepts>
#include <functional>
#include <vector>

#include "Particle.h"
#include "Vector3f.h"
#include "type_traits.hpp"

//
// Defines concepts used in the library to constrain template parameters.
//
// Also defines `getX()`, `getY()` and `getZ()` functions that can be used to get the
// coordinates of a `Locatable` object / pointer / reference.
//

namespace biospring
{

namespace concepts
{

// ======================================================================================
//
// `Locatable` concept.
//
// Defines a type that can be located in 3D space.
// Requires that the type has a `getX()`, `getY()` and `getZ()` methods that returns a
// floating point number.
//
// ======================================================================================

// =====================================================================================
// Locatable concept
// =====================================================================================

// clang-format off

// Object version.
template <typename T>
concept LocatableObject = requires(T a)
{
    { std::invoke(&T::type::getX, a) } -> std::convertible_to<double>;
    { std::invoke(&T::type::getY, a) } -> std::convertible_to<double>;
    { std::invoke(&T::type::getZ, a) } -> std::convertible_to<double>;
};

// Pointer version.
template <typename T>
concept LocatablePointer = requires(T a)
{
    requires std::is_pointer_v<T>;
    { a->getX() } -> std::convertible_to<double>;
    { a->getY() } -> std::convertible_to<double>;
    { a->getZ() } -> std::convertible_to<double>;
};

// Reference version.
template <typename T>
concept LocatableReference = requires(T & a)
{
    { a.getX() } -> std::convertible_to<double>;
    { a.getY() } -> std::convertible_to<double>;
    { a.getZ() } -> std::convertible_to<double>;
};

// clang-format on

// Locatable version.
template <typename T>
concept Locatable = LocatableObject<T> || LocatablePointer<T> || LocatableReference<T>;

// Makes sure that the `Particle` and `Vector3f` classes are `Locatable` (pointers too).
static_assert(Locatable<biospring::spn::Particle>);
static_assert(Locatable<biospring::spn::Particle *>);
static_assert(Locatable<Vector3f>);
static_assert(Locatable<Vector3f *>);
static_assert(Locatable<biospring::spn::Particle &>);
static_assert(Locatable<std::reference_wrapper<biospring::spn::Particle>>);

// ======================================================================================
//
// Use these functions to get the coordinates of a `Locatable` object / pointer / reference.
//
// ======================================================================================

namespace locatable
{

template <Locatable T> double getX(const T & element)
{
    if constexpr (is_reference_wrapper_v<T>)
        return element.get().getX();

    else if constexpr (std::is_pointer_v<T>)
        return element->getX();

    else
        return element.getX();
}

template <Locatable T> double getY(const T & element)
{
    if constexpr (is_reference_wrapper_v<T>)
        return element.get().getY();

    else if constexpr (std::is_pointer_v<T>)
        return element->getY();

    else
        return element.getY();
}

template <Locatable T> double getZ(const T & element)
{
    if constexpr (is_reference_wrapper_v<T>)
        return element.get().getZ();

    else if constexpr (std::is_pointer_v<T>)
        return element->getZ();

    else
        return element.getZ();
}

template <Locatable T> std::array<double, 3> get_position(const T & element)
{
    return {getX(element), getY(element), getZ(element)};
}

} // namespace locatable

// ======================================================================================
//
// `LocatableContainer` concept.
//
// Defines a type that contains `Locatable` objects.
//
// ======================================================================================
template <class ContainerType>
concept LocatableContainer =
    requires(ContainerType a, const ContainerType b) {
        // clang-format off

    // Container requirements.
    typename ContainerType::value_type; // Ensure that the container has a value type.
    requires Locatable<typename ContainerType::value_type>; // Ensure that the value type is Locatable.

    requires std::destructible<typename ContainerType::value_type>;
    requires std::same_as<typename ContainerType::reference, typename ContainerType::value_type &>;
    requires std::same_as<typename ContainerType::const_reference, const typename ContainerType::value_type &>;
    requires std::forward_iterator<typename ContainerType::iterator>;
    requires std::forward_iterator<typename ContainerType::const_iterator>;

    { a.at(std::declval<typename ContainerType::size_type>()) } -> std::same_as<typename ContainerType::reference>;
    { b.at(std::declval<typename ContainerType::size_type>()) } -> std::same_as<typename ContainerType::const_reference>;
    { a.begin() } -> std::same_as<typename ContainerType::iterator>; // Ensure the container has a begin() method returning an iterator
    { a.end() } -> std::same_as<typename ContainerType::iterator>; // Ensure the container has an end() method returning an iterator
    { b.begin() } -> std::same_as<typename ContainerType::const_iterator>; // Ensure the container has a begin() method returning a const_iterator
    { b.end() } -> std::same_as<typename ContainerType::const_iterator>; // Ensure the container has an end() method returning a const_iterator
    { a.size() } -> std::same_as<typename ContainerType::size_type>; // Ensure the container has a size() method returning a size_t type
    { a.empty() } -> std::same_as<bool>;  // Ensure the container has an empty() method returning a bool type

        // clang-format on
    };

static_assert(LocatableContainer<std::vector<biospring::spn::Particle>>);
static_assert(LocatableContainer<std::vector<biospring::spn::Particle *>>);
static_assert(LocatableContainer<std::array<Vector3f, 0>>);
static_assert(LocatableContainer<std::array<Vector3f *, 0>>);
static_assert(LocatableContainer<std::vector<std::reference_wrapper<biospring::spn::Particle>>>);

} // namespace concepts
} // namespace biospring

#endif // __CONCEPTS_HPP__