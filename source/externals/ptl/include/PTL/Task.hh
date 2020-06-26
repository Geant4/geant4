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
// This file defines the task types for TaskManager and ThreadPool
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "Globals.hh"
#include "TaskAllocator.hh"
#include "VTask.hh"

#include <cstdint>
#include <functional>
#include <stdexcept>

namespace PTL
{
class VTaskGroup;
class ThreadPool;

//======================================================================================//

/// \brief The task class is supplied to thread_pool.
template <typename _Ret, typename... _Args>
class PackagedTask : public VTask
{
public:
    typedef PackagedTask<_Ret, _Args...>       this_type;
    typedef std::promise<_Ret>                 promise_type;
    typedef std::future<_Ret>                  future_type;
    typedef std::packaged_task<_Ret(_Args...)> packaged_task_type;
    typedef _Ret                               result_type;
    typedef std::tuple<_Args...>               tuple_type;

public:
    // pass a free function pointer
    template <typename _Func>
    PackagedTask(_Func&& func, _Args&&... args)
    : VTask()
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    template <typename _Func>
    PackagedTask(VTaskGroup* tg, _Func&& func, _Args&&... args)
    : VTask(tg)
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    template <typename _Func>
    PackagedTask(ThreadPool* _pool, _Func&& func, _Args&&... args)
    : VTask(_pool)
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    virtual ~PackagedTask() {}

public:
    // execution operator
    virtual void operator()() override
    {
        details::apply(std::forward<packaged_task_type>(m_ptask),
                       std::forward<tuple_type>(m_args));
    }
    future_type  get_future() { return m_ptask.get_future(); }
    virtual bool is_native_task() const override { return true; }

private:
    packaged_task_type m_ptask;
    tuple_type         m_args;
};

//======================================================================================//

/// \brief The task class is supplied to thread_pool.
template <typename _Ret, typename... _Args>
class Task : public VTask
{
public:
    typedef Task<_Ret, _Args...>               this_type;
    typedef std::promise<_Ret>                 promise_type;
    typedef std::future<_Ret>                  future_type;
    typedef std::packaged_task<_Ret(_Args...)> packaged_task_type;
    typedef _Ret                               result_type;
    typedef std::tuple<_Args...>               tuple_type;

public:
    template <typename _Func>
    Task(_Func&& func, _Args&&... args)
    : VTask()
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    template <typename _Func>
    Task(VTaskGroup* tg, _Func&& func, _Args&&... args)
    : VTask(tg)
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    template <typename _Func>
    Task(ThreadPool* tp, _Func&& func, _Args&&... args)
    : VTask(tp)
    , m_ptask(std::forward<_Func>(func))
    , m_args(std::forward<_Args>(args)...)
    {}

    virtual ~Task() {}

public:
    // execution operator
    virtual void operator()() final
    {
        details::apply(std::forward<packaged_task_type>(m_ptask),
                       std::forward<tuple_type>(m_args));
        // decrements the task-group counter on active tasks
        // when the counter is < 2, if the thread owning the task group is
        // sleeping at the TaskGroup::wait(), it signals the thread to wake
        // up and check if all tasks are finished, proceeding if this
        // check returns as true
        this_type::operator--();
    }

    virtual bool is_native_task() const override { return true; }
    future_type  get_future() { return m_ptask.get_future(); }

private:
    packaged_task_type m_ptask;
    tuple_type         m_args;
};

//======================================================================================//

/// \brief The task class is supplied to thread_pool.
template <typename _Ret>
class Task<_Ret, void> : public VTask
{
public:
    typedef Task<_Ret>                 this_type;
    typedef std::promise<_Ret>         promise_type;
    typedef std::future<_Ret>          future_type;
    typedef std::packaged_task<_Ret()> packaged_task_type;
    typedef _Ret                       result_type;

public:
    template <typename _Func>
    Task(_Func&& func)
    : VTask()
    , m_ptask(std::forward<_Func>(func))
    {}

    template <typename _Func>
    Task(VTaskGroup* tg, _Func&& func)
    : VTask(tg)
    , m_ptask(std::forward<_Func>(func))
    {}

    template <typename _Func>
    Task(ThreadPool* tp, _Func&& func)
    : VTask(tp)
    , m_ptask(std::forward<_Func>(func))
    {}

    virtual ~Task() {}

public:
    // execution operator
    virtual void operator()() final
    {
        m_ptask();
        // decrements the task-group counter on active tasks
        // when the counter is < 2, if the thread owning the task group is
        // sleeping at the TaskGroup::wait(), it signals the thread to wake
        // up and check if all tasks are finished, proceeding if this
        // check returns as true
        this_type::operator--();
    }

    virtual bool is_native_task() const override { return true; }
    future_type  get_future() { return m_ptask.get_future(); }

private:
    packaged_task_type m_ptask;
};

//======================================================================================//

/// \brief The task class is supplied to thread_pool.
template <>
class Task<void, void> : public VTask
{
public:
    typedef void                       _Ret;
    typedef Task<void, void>           this_type;
    typedef std::promise<_Ret>         promise_type;
    typedef std::future<_Ret>          future_type;
    typedef std::packaged_task<_Ret()> packaged_task_type;
    typedef _Ret                       result_type;

public:
    template <typename _Func>
    explicit Task(_Func&& func)
    : VTask()
    , m_ptask(std::forward<_Func>(func))
    {}

    template <typename _Func>
    Task(VTaskGroup* tg, _Func&& func)
    : VTask(tg)
    , m_ptask(std::forward<_Func>(func))
    {}

    template <typename _Func>
    Task(ThreadPool* tp, _Func&& func)
    : VTask(tp)
    , m_ptask(std::forward<_Func>(func))
    {}

    virtual ~Task() {}

public:
    // execution operator
    virtual void operator()() final
    {
        m_ptask();
        // decrements the task-group counter on active tasks
        // when the counter is < 2, if the thread owning the task group is
        // sleeping at the TaskGroup::wait(), it signals the thread to wake
        // up and check if all tasks are finished, proceeding if this
        // check returns as true
        this_type::operator--();
    }

    virtual bool is_native_task() const override { return true; }
    future_type  get_future() { return m_ptask.get_future(); }

private:
    packaged_task_type m_ptask;
};

//======================================================================================//

}  // namespace PTL
