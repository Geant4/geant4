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

#include <algorithm>  // Retrieve definitions of min/max

// Include base types
#include "PTL/Types.hh"

// Global utility functions
#include "PTL/Utility.hh"

#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>

#if !defined(PTL_NO_SANITIZE_THREAD)
// expect that sanitizer is from compiler which supports __has_attribute
#    if defined(__has_attribute)
#        if __has_attribute(no_sanitize)
#            define PTL_NO_SANITIZE_THREAD __attribute__((no_sanitize("thread")))
#        else
#            define PTL_NO_SANITIZE_THREAD
#        endif
#    elif defined(__clang__) || defined(__GNUC__)
#        define PTL_NO_SANITIZE_THREAD __attribute__((no_sanitize("thread")))
#    else
// otherwise, make blank
#        define PTL_NO_SANITIZE_THREAD
#    endif
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
    using __type = Index_tuple<0>;
};

template <>
struct Build_index_tuple<0>
{
    using __type = Index_tuple<>;
};

/// Class template integer_sequence
template <typename Tp, Tp... Idx>
struct integer_sequence
{
    using value_type = Tp;
    static constexpr size_t size() noexcept { return sizeof...(Idx); }
};

template <typename Tp, Tp NumT, typename ISeq = typename Build_index_tuple<NumT>::__type>
struct Make_integer_sequence;

template <typename Tp, Tp NumT, size_t... Idx>
struct Make_integer_sequence<Tp, NumT, Index_tuple<Idx...>>
{
    static_assert(NumT >= 0, "Cannot make integer sequence of negative length");

    using __type = integer_sequence<Tp, static_cast<Tp>(Idx)...>;
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
apply(FnT&& _func, TupleT _args, impl::index_sequence<Idx...>)
    -> decltype(std::forward<FnT>(_func)(std::get<Idx>(std::move(_args))...))
{
    // GCC 5.3 warns about unused variable _args when the index sequence is empty
#if defined(__GNUC__) && (__GNUC__ < 6)
    if(sizeof...(Idx) == 0)
    {
        ConsumeParameters(_args);
    }
#endif
    return std::forward<FnT>(_func)(std::get<Idx>(std::move(_args))...);
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
apply(FnT&& _func, TupleT&& _args)
{
    using tuple_type = typename std::decay<TupleT>::type;
    constexpr auto N = std::tuple_size<tuple_type>::value;
    impl::apply(std::forward<FnT>(_func), std::forward<TupleT>(_args),
                impl::make_index_sequence<N>{});
}

//--------------------------------------------------------------------------------------//

}  // namespace mpl

}  // namespace PTL
