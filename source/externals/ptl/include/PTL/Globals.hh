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
        ::PTL::mpl::consume_parameters(                                                  \
            ::std::initializer_list<int>{ (__VA_ARGS__, 0)... })
#endif

namespace PTL
{
template <typename T>
using decay_t = typename std::decay<T>::type;

template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

// for pre-C++14 tuple expansion to arguments
namespace mpl
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
template <size_t... Indexes>
struct Index_tuple
{};

// Concatenates two Index_tuples.
template <typename Itup1, typename Itup2>
struct Itup_cat;

template <size_t... Ind1, size_t... Ind2>
struct Itup_cat<Index_tuple<Ind1...>, Index_tuple<Ind2...>>
{
    using __type = Index_tuple<Ind1..., (Ind2 + sizeof...(Ind1))...>;
};

// Builds an Index_tuple<0, 1, 2, ..., NumT-1>.
template <size_t NumT>
struct Build_index_tuple
: Itup_cat<typename Build_index_tuple<NumT / 2>::__type,
           typename Build_index_tuple<NumT - NumT / 2>::__type>
{};

template <>
struct Build_index_tuple<1>
{
    typedef Index_tuple<0> __type;
};

template <>
struct Build_index_tuple<0>
{
    typedef Index_tuple<> __type;
};

/// Class template integer_sequence
template <typename Tp, Tp... Idx>
struct integer_sequence
{
    typedef Tp              value_type;
    static constexpr size_t size() noexcept { return sizeof...(Idx); }
};

template <typename Tp, Tp NumT, typename ISeq = typename Build_index_tuple<NumT>::__type>
struct Make_integer_sequence;

template <typename Tp, Tp NumT, size_t... Idx>
struct Make_integer_sequence<Tp, NumT, Index_tuple<Idx...>>
{
    static_assert(NumT >= 0, "Cannot make integer sequence of negative length");

    typedef integer_sequence<Tp, static_cast<Tp>(Idx)...> __type;
};

/// Alias template make_integer_sequence
template <typename Tp, Tp NumT>
using make_integer_sequence = typename Make_integer_sequence<Tp, NumT>::__type;

/// Alias template index_sequence
template <size_t... Idx>
using index_sequence = integer_sequence<size_t, Idx...>;

/// Alias template make_index_sequence
template <size_t NumT>
using make_index_sequence = make_integer_sequence<size_t, NumT>;

/// Alias template index_sequence_for
template <typename... Types>
using index_sequence_for = make_index_sequence<sizeof...(Types)>;

template <size_t Idx, typename Tup>
using index_type_t = decay_t<decltype(std::get<Idx>(std::declval<Tup>()))>;

template <typename FnT, typename TupleT, size_t... Idx>
static inline auto
apply(FnT&& __f, TupleT&& __t, impl::index_sequence<Idx...>)
    -> decltype(std::forward<FnT>(__f)(std::get<Idx>(__t)...))
{
    return std::forward<FnT>(__f)(std::get<Idx>(__t)...);
}

//--------------------------------------------------------------------------------------//

}  // namespace impl

//--------------------------------------------------------------------------------------//

/// Alias template index_sequence
template <size_t... Idx>
using index_sequence = impl::integer_sequence<size_t, Idx...>;

/// Alias template make_index_sequence
template <size_t NumT>
using make_index_sequence = impl::make_integer_sequence<size_t, NumT>;

/// Alias template index_sequence_for
template <typename... Types>
using index_sequence_for = impl::make_index_sequence<sizeof...(Types)>;

template <typename FnT, typename TupleT>
static inline void
apply(FnT&& __f, TupleT&& __t)
{
    constexpr auto N = std::tuple_size<TupleT>::value;
    impl::apply(std::forward<FnT>(__f), std::forward<TupleT>(__t),
                impl::make_index_sequence<N>{});
}

//--------------------------------------------------------------------------------------//

}  // namespace mpl

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
