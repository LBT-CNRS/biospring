#ifndef __TYPE_TRAITS_HPP__
#define __TYPE_TRAITS_HPP__

//
// Defines a type trait to check if a type is a std::reference_wrapper.
//

#include <type_traits>

namespace biospring
{
struct _A {};

template<typename T>
struct is_reference_wrapper : std::false_type {};

template<typename T>
struct is_reference_wrapper<std::reference_wrapper<T>> : std::true_type {};

template< class T >
inline constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;

static_assert(
       ! is_reference_wrapper_v<_A>
    && ! is_reference_wrapper_v<_A *>
    && ! is_reference_wrapper_v<_A &>
    &&   is_reference_wrapper_v<std::reference_wrapper<_A>>
    &&   is_reference_wrapper_v<std::reference_wrapper<_A const>>
);


} // namespace biospring


#endif // __TYPE_TRAITS_HPP__

