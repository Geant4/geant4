//  MIT License
//  Copyright (c) 2020 Jonathan R. Madsen
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//  SOFTWARE.

#pragma once

#include <cstring>
#include <functional>
#include <tuple>
#include <type_traits>

namespace PTL
{
template <typename Tp, typename Up>
struct SmallerThanT
{
    static constexpr bool value = sizeof(Tp) < sizeof(Up);
};

//======================================================================================//

// CTValue == compile-time value
//
// useful to work with sequences of compile-time values, such as the bounds of a
// multidimensional array or indices into another typelist.

template <typename Tp, Tp Value>
struct CTValue
{
    static constexpr Tp value = Value;
};

//======================================================================================//

template <std::size_t Height, typename Tp,
          bool = std::is_class<Tp>::value && !std::is_final<Tp>::value>
class TupleElt;

//--------------------------------------------------------------------------------------//
//  specialization that does not satisfy is_class<> and not is_final<>
//
template <std::size_t Height, typename Tp>
class TupleElt<Height, Tp, false>
{
    Tp value;

public:
    TupleElt() = default;

    template <typename Up>
    TupleElt(Up&& other)
    : value(std::forward<Up>(other))
    {}

    Tp&       get() { return value; }
    Tp const& get() const { return value; }
};

//--------------------------------------------------------------------------------------//
//  specialization that does satisfy is_class<> and not is_final<>
//
template <std::size_t Height, typename Tp>
class TupleElt<Height, Tp, true> : private Tp
{
public:
    TupleElt() = default;

    template <typename Up>
    explicit TupleElt(Up&& other)
    : Tp(std::forward<Up>(other))
    {}

    Tp&       get() { return *this; }
    const Tp& get() const { return *this; }
};

//======================================================================================//

template <unsigned H, typename T>
T&
get_height(TupleElt<H, T>& te)
{
    return te.get();
}

//======================================================================================//

template <typename... Types>
class Tuple;

//--------------------------------------------------------------------------------------//

// recursive case:
template <typename Head, typename... Tail>
class Tuple<Head, Tail...>
: private TupleElt<sizeof...(Tail), Head>
, private Tuple<Tail...>
{
    template <unsigned I, typename... Elements>
    friend auto get(Tuple<Elements...>& t)
        -> decltype(get_height<sizeof...(Elements) - I - 1>(t));

public:
    Head&           head() { return static_cast<HeadElt*>(this)->get(); }
    const Head&     head() const { return static_cast<HeadElt const*>(this)->get(); }
    Tuple<Tail...>& tail() { return *this; }
    const Tuple<Tail...>& tail() const { return *this; }

private:
    using HeadElt = TupleElt<sizeof...(Tail), Head>;
};

//--------------------------------------------------------------------------------------//
// basis case:
template <>
class Tuple<>
{
    // no storage required
};

//--------------------------------------------------------------------------------------//

template <unsigned I, typename... Elements>
auto
get(Tuple<Elements...>& t) -> decltype(get_height<sizeof...(Elements) - I - 1>(t))
{
    return get_height<sizeof...(Elements) - I - 1>(t);
}

template <typename... Tp>
using TypeList = Tuple<Tp...>;

//======================================================================================//

template <typename List>
class IsEmpty
{
public:
    static constexpr bool value = false;
};

template <>
class IsEmpty<Tuple<>>
{
public:
    static constexpr bool value = true;
};

template <typename List, typename NewElement>
class PushBackT;

template <typename... Elements, typename NewElement>
class PushBackT<Tuple<Elements...>, NewElement>
{
public:
    using Type = Tuple<Elements..., NewElement>;
};

template <typename List, typename NewElement>
using PushBack = typename PushBackT<List, NewElement>::Type;

template <typename List>
class PopFrontT;

template <typename Head, typename... Tail>
class PopFrontT<std::tuple<Head, Tail...>>
{
public:
    using Type = std::tuple<Tail...>;
};

template <typename List>
using PopFront = typename PopFrontT<List>::Type;

template <typename List, typename Element>
class PushFrontT;

template <typename... Types, typename Element>
class PushFrontT<std::tuple<Types...>, Element>
{
public:
    using Type = std::tuple<Element, Types...>;
};

template <typename List, typename NewElement>
using PushFront = typename PushFrontT<List, NewElement>::Type;

template <typename... Types, typename V>
PushFront<std::tuple<Types...>, V>
pushFront(std::tuple<Types...> const& tuple, V const& value)
{
    return PushFront<std::tuple<Types...>, V>(value, tuple);
}

template <typename T, T... Values>
struct Valuelist
{};

template <typename... Types>
struct Front;

template <typename FrontT, typename... Types>
struct Front<FrontT, Types...>
{
    using Type = FrontT;
};

template <typename... Types>
struct Back;

template <typename BackT, typename... Types>
struct Back<Types..., BackT>
{                        // ERROR: pack expansion not at the end of
    using Type = BackT;  //       template argument list
};

//======================================================================================//

template <typename List, template <typename T> class MetaFun,
          bool Empty = IsEmpty<List>::value>
class TransformT;

// recursive case:
template <typename List, template <typename T> class MetaFun>
class TransformT<List, MetaFun, false>
: public PushFrontT<typename TransformT<PopFront<List>, MetaFun>::Type,
                    typename MetaFun<Front<List>>::Type>
{};

// basis case:
template <typename List, template <typename T> class MetaFun>
class TransformT<List, MetaFun, true>
{
public:
    using Type = List;
};

template <typename List, template <typename T> class MetaFun>
using Transform = typename TransformT<List, MetaFun>::Type;

template <typename... Elements, template <typename T> class MetaFun>
class TransformT<Tuple<Elements...>, MetaFun, false>
{
public:
    using Type = Tuple<typename MetaFun<Elements>::Type...>;
};

template <size_t N>
struct ForwardTupleAsArgs
{
    template <typename Func, typename Head, typename... Tail>
    static inline auto forward(Func&& func, Head&& head, Tail&&... tail)
        -> decltype(ForwardTupleAsArgs<N - 1>::forward(
            std::forward<Func>(func), std::forward<Head>(head),
            std::get<N - 1>(std::forward<Head>(head)), std::forward<Tail>(tail)...))
    {
        return ForwardTupleAsArgs<N - 1>::forward(
            std::forward<Func>(func), std::forward<Head>(head),
            std::get<N - 1>(std::forward<Head>(head)), std::forward<Tail>(tail)...);
    }
};

//======================================================================================//

template <>
struct ForwardTupleAsArgs<0>
{
    template <typename Func, typename Head, typename... Tail>
    static inline auto forward(Func&& func, Head&&, Tail&&... tail)
        -> decltype(std::forward<Func>(func)(std::forward<Tail>(tail)...))
    {
        return std::forward<Func>(func)(std::forward<Tail>(tail)...);
    }
};

//======================================================================================//

template <size_t N>
struct ForEachTupleArg
{
    template <typename Func, typename Head>
    static inline auto apply(Func&& func, Head&& head)
    {
        std::forward<Func>(func)(std::forward<Head>(head));
    }

    template <typename Func, typename Head, typename... Tail>
    static inline void apply(Func&& func, Head&& head, Tail&&... tail)
    {
        std::forward<Func>(func)(std::forward<Head>(head));
        ForEachTupleArg<N - 1>::apply(std::forward<Func>(func),
                                      std::forward<Tail>(tail)...);
    }
};

//======================================================================================//

template <typename Func, typename Tuple>
void
for_each_tuple_arg(Func&& func, Tuple&& _tuple)
{
    ForEachTupleArg<std::tuple_size<Tuple>::value>::apply(
        std::forward<Func>(func), std::forward<std::tuple>(_tuple));
}

//======================================================================================//

template <typename Func, typename Tuple, std::size_t Head>
inline auto
InvokeSequence_impl(const Func& func, const Tuple& data)
{
    func(std::get<Head>(data));
}

//======================================================================================//

template <typename Func, typename Tuple, std::size_t Head, std::size_t... Tail>
inline auto
InvokeSequence_impl(const Func& func, const Tuple& data)
{
    func(std::get<Head>(data));
    InvokeSequence_impl<Func, Tuple, Tail...>(func, data);
}

//======================================================================================//

template <typename Func, typename Tuple, std::size_t N = std::tuple_size<Tuple>::value,
          typename Indices = std::make_index_sequence<N>>
inline auto
InvokeSequence(const Func& func, const Tuple& data)
{
    return InvokeSequence_impl(func, data);
}

//======================================================================================//

// Convert array into a tuple
template <typename Container, std::size_t... N>
inline auto
ContainerToTuple_impl(const Container& tasks, std::index_sequence<N...>)
{
    return std::make_tuple(tasks[N]...);
}

//======================================================================================//

template <std::size_t N, typename Container,
          typename Indices = std::make_index_sequence<N>>
inline auto
ContainerToTuple(const Container& tasks)
{
    return ContainerToTuple_impl(tasks, Indices{});
}

//======================================================================================//

template <typename... Args>
struct tuple_subset
{
    static std::tuple<> get(const std::tuple<>&) { return std::tuple<>{}; }

    template <typename... SubArgs>
    static std::tuple<SubArgs...> get(const std::tuple<Args...>& t)
    {
        return std::tuple<SubArgs...>{ std::get<SubArgs>(t)... };
    }
};

template <typename Head>
inline void
tuple_transform(const std::function<void(const Head&)>& pred,
                const std::tuple<Head>&                 data)
{
    pred(std::get<0>(data));
}

template <typename Head, typename... Tail>
inline void
tuple_transform(const std::function<void(const Head&)>& pred,
                const std::tuple<Head, Tail...>&        data)
{
    pred(std::get<0>(data));
    auto subset = tuple_subset<Head, Tail...>::template get<Tail...>(data);
    tuple_transform<Tail...>(pred, std::forward<decltype(subset)>(subset));
}

template <typename Head, typename... Tail>
struct transform_tuple
{
    using Function = std::function<void(Head)>;
    template <typename TupleType>
    static void apply(const Function& func, const TupleType& t)
    {
        func(std::get<0>(t));
        PopFront<TupleType> nt = std::tuple<Tail...>{ std::get<Tail>(t)... };
        transform_tuple<Tail...>::apply(func, nt);
    }
};

template <typename Head>
struct transform_tuple<Head>
{
    using Function = std::function<void(Head)>;
    template <typename TupleType>
    static void apply(const Function& func, const TupleType& t)
    {
        func(std::get<0>(t));
    }
};

//======================================================================================//

template <typename Func, typename... Elements, unsigned... Indices>
auto
applyImpl(Func func, std::tuple<Elements...> const& t, Valuelist<unsigned, Indices...>)
    -> decltype(func(std::get<Indices>(t)...))
{
    return func(std::get<Indices>(t)...);
}

template <typename Func, typename... Elements, unsigned N = sizeof...(Elements)>
auto
apply(Func func, std::tuple<Elements...> const& t)
    -> decltype(applyImpl(func, t, std::make_index_sequence<N>()))
{
    return applyImpl(func, t, std::make_index_sequence<N>());
}

}  // namespace PTL
