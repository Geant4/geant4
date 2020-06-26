//
// MIT License
// Copyright (c) 2020 Jonathan R. Madsen
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
// "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#pragma once

#ifndef FALSE
#    define FALSE 0
#endif

#ifndef TRUE
#    define TRUE 1
#endif

#include <algorithm>  // Retrieve definitions of min/max

// Include base types
#include "PTL/Types.hh"

// Global utility functions
#include "PTL/Utility.hh"

#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>

#if !defined(PTL_FOLD_EXPRESSION)
#    define PTL_FOLD_EXPRESSION(...)                                                     \
        ::PTL::details::consume_parameters(                                              \
            ::std::initializer_list<int>{ (__VA_ARGS__, 0)... })
#endif

namespace PTL
{
template <typename T>
using decay_t = typename std::decay<T>::type;

template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

// for pre-C++14 tuple expansion to arguments
namespace details
{
//--------------------------------------------------------------------------------------//

template <typename... Args>
void
consume_parameters(Args&&...)
{}

//--------------------------------------------------------------------------------------//

namespace impl
{
//--------------------------------------------------------------------------------------//
// Stores a tuple of indices.  Used by tuple and pair, and by bind() to
// extract the elements in a tuple.
template <size_t... _Indexes>
struct _Index_tuple
{};

// Concatenates two _Index_tuples.
template <typename _Itup1, typename _Itup2>
struct _Itup_cat;

template <size_t... _Ind1, size_t... _Ind2>
struct _Itup_cat<_Index_tuple<_Ind1...>, _Index_tuple<_Ind2...>>
{
    using __type = _Index_tuple<_Ind1..., (_Ind2 + sizeof...(_Ind1))...>;
};

// Builds an _Index_tuple<0, 1, 2, ..., _Num-1>.
template <size_t _Num>
struct _Build_index_tuple
: _Itup_cat<typename _Build_index_tuple<_Num / 2>::__type,
            typename _Build_index_tuple<_Num - _Num / 2>::__type>
{};

template <>
struct _Build_index_tuple<1>
{
    typedef _Index_tuple<0> __type;
};

template <>
struct _Build_index_tuple<0>
{
    typedef _Index_tuple<> __type;
};

/// Class template integer_sequence
template <typename _Tp, _Tp... _Idx>
struct integer_sequence
{
    typedef _Tp             value_type;
    static constexpr size_t size() noexcept { return sizeof...(_Idx); }
};

template <typename _Tp, _Tp _Num,
          typename _ISeq = typename _Build_index_tuple<_Num>::__type>
struct _Make_integer_sequence;

template <typename _Tp, _Tp _Num, size_t... _Idx>
struct _Make_integer_sequence<_Tp, _Num, _Index_tuple<_Idx...>>
{
    static_assert(_Num >= 0, "Cannot make integer sequence of negative length");

    typedef integer_sequence<_Tp, static_cast<_Tp>(_Idx)...> __type;
};

/// Alias template make_integer_sequence
template <typename _Tp, _Tp _Num>
using make_integer_sequence = typename _Make_integer_sequence<_Tp, _Num>::__type;

/// Alias template index_sequence
template <size_t... _Idx>
using index_sequence = integer_sequence<size_t, _Idx...>;

/// Alias template make_index_sequence
template <size_t _Num>
using make_index_sequence = make_integer_sequence<size_t, _Num>;

/// Alias template index_sequence_for
template <typename... _Types>
using index_sequence_for = make_index_sequence<sizeof...(_Types)>;

template <size_t I, typename Tup>
using index_type_t = decay_t<decltype(std::get<I>(std::declval<Tup>()))>;

template <typename _Fn, typename _Tuple, size_t... _Idx>
static inline void
apply(_Fn&& __f, _Tuple&& __t, impl::index_sequence<_Idx...>)
{
    std::forward<_Fn>(__f)(std::get<_Idx>(std::forward<_Tuple>(__t))...);
    // __f(std::get<_Idx>(std::forward<_Tuple>(__t))...);
}

//--------------------------------------------------------------------------------------//

}  // namespace impl

//--------------------------------------------------------------------------------------//

/// Alias template index_sequence
template <size_t... _Idx>
using index_sequence = impl::integer_sequence<size_t, _Idx...>;

/// Alias template make_index_sequence
template <size_t _Num>
using make_index_sequence = impl::make_integer_sequence<size_t, _Num>;

/// Alias template index_sequence_for
template <typename... _Types>
using index_sequence_for = impl::make_index_sequence<sizeof...(_Types)>;

template <typename _Fn, typename _Tuple>
static inline void
apply(_Fn&& __f, _Tuple&& __t)
{
    constexpr auto _N = std::tuple_size<_Tuple>::value;
    impl::apply(std::forward<_Fn>(__f), std::forward<_Tuple>(__t),
                impl::make_index_sequence<_N>{});
}

//--------------------------------------------------------------------------------------//

}  // namespace details

namespace thread_pool
{
namespace state
{
static const short STARTED = 0;
static const short PARTIAL = 1;
static const short STOPPED = 2;
static const short NONINIT = 3;

}  // namespace state

}  // namespace thread_pool

//--------------------------------------------------------------------------------------//

}  // namespace PTL
