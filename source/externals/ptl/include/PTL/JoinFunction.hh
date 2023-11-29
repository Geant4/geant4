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
//
// ---------------------------------------------------------------
// Tasking class header file
//
// Class Description:
//
// This is the join function used by task groups
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/Types.hh"

#include <functional>
#include <utility>

namespace PTL
{
template <typename JoinT, typename JoinArg>
struct JoinFunction
{
public:
    using Type = std::function<JoinT(JoinT&, JoinArg&&)>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    template <typename... Args>
    JoinT& operator()(Args&&... args)
    {
        return std::move(m_func(std::forward<Args>(args)...));
    }

private:
    Type m_func = [](JoinT& lhs, JoinArg&&) { return lhs; };
};

//--------------------------------------------------------------------------------------//

template <typename JoinArg>
struct JoinFunction<void, JoinArg>
{
public:
    using Type = std::function<void(JoinArg)>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    template <typename... Args>
    void operator()(Args&&... args)
    {
        m_func(std::forward<Args>(args)...);
    }

private:
    Type m_func = [](JoinArg) {};
};

//--------------------------------------------------------------------------------------//

template <>
struct JoinFunction<void, void>
{
public:
    using Type = std::function<void()>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    void operator()() { m_func(); }

private:
    Type m_func = []() {};
};

}  // namespace PTL
