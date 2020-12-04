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
// ---------------------------------------------------------------
// Tasking class header file
//
// Class Description:
//
// This file creates a class for handling the wrapping of functions
// into task objects and submitting to thread pool
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/Globals.hh"
#include "PTL/TBBTaskGroup.hh"
#include "PTL/Task.hh"
#include "PTL/TaskGroup.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Threading.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>

namespace PTL
{
//======================================================================================//

class TaskManager
{
public:
    typedef TaskManager           this_type;
    typedef ThreadPool::size_type size_type;

public:
    // Constructor and Destructors
    explicit TaskManager(ThreadPool*);
    virtual ~TaskManager();

    TaskManager(const this_type&) = delete;
    TaskManager(this_type&&)      = default;
    this_type& operator=(const this_type&) = delete;
    this_type& operator=(this_type&&) = default;

public:
    /// get the singleton pointer
    static TaskManager* GetInstance();
    static TaskManager* GetInstanceIfExists();
    static unsigned     ncores() { return std::thread::hardware_concurrency(); }

public:
    //------------------------------------------------------------------------//
    // return the thread pool
    inline ThreadPool* thread_pool() const { return m_pool; }

    //------------------------------------------------------------------------//
    // return the number of threads in the thread pool
    inline size_type size() const { return m_pool->size(); }

    //------------------------------------------------------------------------//
    // kill all the threads
    inline void finalize() { m_pool->destroy_threadpool(); }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    // direct insertion of a task
    //------------------------------------------------------------------------//
    template <typename... Args>
    void exec(Task<Args...>* _task)
    {
        m_pool->add_task(_task);
    }

    //------------------------------------------------------------------------//
    // direct insertion of a packaged_task
    //------------------------------------------------------------------------//
    template <typename RetT, typename FuncT, typename... Args>
    std::future<RetT> async(FuncT&& func, Args&&... args)
    {
        typedef PackagedTask<RetT, Args...>  task_type;
        typedef task_type*                   task_pointer;

        task_pointer _ptask =
            new task_type(std::forward<FuncT>(func), std::forward<Args>(args)...);
        std::future<RetT> _f = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename FuncT>
    std::future<RetT> async(FuncT&& func)
    {
        typedef PackagedTask<RetT> task_type;
        typedef task_type*         task_pointer;

        task_pointer      _ptask = new task_type(std::forward<FuncT>(func));
        std::future<RetT> _f     = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//
    template <typename FuncT, typename... Args>
    auto async(FuncT&& func, Args... args)
        -> std::future<decay_t<decltype(func(args...))>>
    {
        using RetT = decay_t<decltype(func(args...))>;
        typedef PackagedTask<RetT, Args...> task_type;

        auto _ptask =
            new task_type(std::forward<FuncT>(func), std::forward<Args>(args)...);
        auto _f = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    // public wrap functions
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT, typename... Args>
    Task<RetT, ArgT, Args...>* wrap(TaskGroup<RetT, ArgT>& tg, FuncT&& func,
                                    Args&&... args)
    {
        return tg.wrap(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT>
    Task<RetT, ArgT>* wrap(TaskGroup<RetT, ArgT>& tg, FuncT&& func)
    {
        return tg.wrap(std::forward<FuncT>(func));
    }

public:
    //------------------------------------------------------------------------//
    // public exec functions
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT, typename... Args>
    void exec(TaskGroup<RetT, ArgT>& tg, FuncT&& func, Args&&... args)
    {
        tg.exec(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT>
    void exec(TaskGroup<RetT, ArgT>& tg, FuncT&& func)
    {
        tg.exec(std::forward<FuncT>(func));
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT, typename... Args>
    void rexec(TaskGroup<RetT, ArgT>& tg, FuncT&& func, Args&&... args)
    {
        tg.exec(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT>
    void rexec(TaskGroup<RetT, ArgT>& tg, FuncT&& func)
    {
        tg.exec(std::forward<FuncT>(func));
    }
    //------------------------------------------------------------------------//
    // public exec functions (void specializations)
    //------------------------------------------------------------------------//
    template <typename FuncT, typename... Args>
    void rexec(TaskGroup<void, void>& tg, FuncT&& func, Args&&... args)
    {
        tg.exec(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename FuncT>
    void rexec(TaskGroup<void, void>& tg, FuncT&& func)
    {
        tg.exec(std::forward<FuncT>(func));
    }
    //------------------------------------------------------------------------//

#if defined(PTL_USE_TBB)
    //------------------------------------------------------------------------//
    // public wrap functions using TBB tasks
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT, typename... Args>
    Task<RetT, ArgT, Args...>* wrap(TBBTaskGroup<RetT, ArgT>& tg, FuncT&& func,
                                    Args&&... args)
    {
        return tg.wrap(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT>
    Task<RetT, ArgT>* wrap(TBBTaskGroup<RetT, ArgT>& tg, FuncT&& func)
    {
        return tg.wrap(std::forward<FuncT>(func));
    }

    //------------------------------------------------------------------------//
    // public exec functions using TBB tasks
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT, typename... Args>
    void exec(TBBTaskGroup<RetT, ArgT>& tg, FuncT&& func, Args&&... args)
    {
        tg.exec(std::forward<FuncT>(func), std::forward<Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename RetT, typename ArgT, typename FuncT>
    void exec(TBBTaskGroup<RetT, ArgT>& tg, FuncT&& func)
    {
        tg.exec(std::forward<FuncT>(func));
    }
    //------------------------------------------------------------------------//
#endif

protected:
    // Protected variables
    ThreadPool* m_pool = nullptr;

private:
    static TaskManager*& fgInstance();
};

}  // namespace PTL
//======================================================================================//

#include "TaskRunManager.hh"

//--------------------------------------------------------------------------------------//

inline PTL::TaskManager*&
PTL::TaskManager::fgInstance()
{
    static thread_local TaskManager* _instance = nullptr;
    return _instance;
}

//--------------------------------------------------------------------------------------//

inline PTL::TaskManager*
PTL::TaskManager::GetInstance()
{
    if(!fgInstance())
    {
        auto nthreads = std::thread::hardware_concurrency();
        std::cout << "Allocating mad::TaskManager with " << nthreads << " thread(s)..."
                  << std::endl;
        new TaskManager(TaskRunManager::GetMasterRunManager()->GetThreadPool());
    }
    return fgInstance();
}

//--------------------------------------------------------------------------------------//

inline PTL::TaskManager*
PTL::TaskManager::GetInstanceIfExists()
{
    return fgInstance();
}

//--------------------------------------------------------------------------------------//

inline PTL::TaskManager::TaskManager(ThreadPool* _pool)
: m_pool(_pool)
{
    if(!fgInstance())
        fgInstance() = this;
}

//--------------------------------------------------------------------------------------//

inline PTL::TaskManager::~TaskManager()
{
    finalize();
    if(fgInstance() == this)
        fgInstance() = nullptr;
}

//======================================================================================//
