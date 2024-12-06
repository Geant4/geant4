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

#include <functional>
#include <utility>

namespace PTL
{
struct ScopeDestructor
{
    template <typename FuncT>
    ScopeDestructor(FuncT&& _func)
    : m_functor(std::forward<FuncT>(_func))
    {}

    // delete copy operations
    ScopeDestructor(const ScopeDestructor&) = delete;
    ScopeDestructor& operator=(const ScopeDestructor&) = delete;

    // allow move operations
    ScopeDestructor(ScopeDestructor&& rhs) noexcept
    : m_functor(std::move(rhs.m_functor))
    {
        rhs.m_functor = []() {};
    }

    ScopeDestructor& operator=(ScopeDestructor&& rhs) noexcept
    {
        if(this != &rhs)
        {
            m_functor     = std::move(rhs.m_functor);
            rhs.m_functor = []() {};
        }
        return *this;
    }

    ~ScopeDestructor() { m_functor(); }

private:
    std::function<void()> m_functor = []() {};
};

}  // namespace PTL
