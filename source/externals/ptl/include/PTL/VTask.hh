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
// Tasking class header file
//
// Class Description:
//
// This file creates an abstract base class for the thread-pool tasking
// system
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/AutoLock.hh"
#include "PTL/Globals.hh"
#include "PTL/Threading.hh"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <string>
#include <thread>
#include <tuple>
#include <utility>

namespace PTL
{
class ThreadPool;

//======================================================================================//

/// \brief VTask is the abstract class stored in thread_pool
class VTask
{
public:
    typedef std::thread::id       tid_type;
    typedef size_t                size_type;
    typedef std::function<void()> void_func_t;

public:
    VTask(bool _is_native, intmax_t _depth);

    VTask()          = default;
    virtual ~VTask() = default;

    VTask(const VTask&) = delete;
    VTask& operator=(const VTask&) = delete;

    VTask(VTask&&) = default;
    VTask& operator=(VTask&&) = default;

public:
    // execution operator
    virtual void operator()() = 0;

public:
    bool     is_native_task() const { return m_is_native; }
    intmax_t depth() const { return m_depth; }

protected:
    bool        m_is_native = false;
    intmax_t    m_depth     = 0;
    void_func_t m_func      = []() {};
};

//======================================================================================//

}  // namespace PTL
