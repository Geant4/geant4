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
#include "PTL/TaskAllocator.hh"
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
class VTaskGroup;
class ThreadPool;

//======================================================================================//

/// \brief VTask is the abstract class stored in thread_pool
class VTask
{
public:
    typedef std::thread::id       tid_type;
    typedef size_t                size_type;
    typedef VTask                 this_type;
    typedef std::atomic_uintmax_t count_t;
    typedef VTask*                iterator;
    typedef const VTask*          const_iterator;
    typedef std::function<void()> void_func_t;

public:
    VTask();
    explicit VTask(VTaskGroup* task_group);
    explicit VTask(ThreadPool* pool);
    virtual ~VTask();

public:
    // execution operator
    virtual void operator()() = 0;

public:
    // used by thread_pool
    void                operator--();
    virtual bool        is_native_task() const;
    virtual ThreadPool* pool() const;
    VTaskGroup*         group() const { return m_group; }

public:
    // used by task tree
    iterator begin() { return this; }
    iterator end() { return this + 1; }

    const_iterator begin() const { return this; }
    const_iterator end() const { return this + 1; }

    const_iterator cbegin() const { return this; }
    const_iterator cend() const { return this + 1; }

    intmax_t&       depth() { return m_depth; }
    const intmax_t& depth() const { return m_depth; }

protected:
    static tid_type this_tid() { return std::this_thread::get_id(); }

protected:
    intmax_t    m_depth;
    VTaskGroup* m_group;
    ThreadPool* m_pool;
    void_func_t m_func = []() {};
};

//======================================================================================//

}  // namespace PTL
