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

    template <bool _Bp, typename _Tp = void>
    using enable_if_t = typename std::enable_if<_Bp, _Tp>::type;

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
    template <typename... _Args>
    void exec(Task<_Args...>* _task)
    {
        m_pool->add_task(_task);
    }

    //------------------------------------------------------------------------//
    // direct insertion of a packaged_task
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Func, typename... _Args>
    std::future<_Ret> async(_Func&& func, _Args&&... args)
    {
        typedef PackagedTask<_Ret, _Args...> task_type;
        typedef task_type*                   task_pointer;

        task_pointer _ptask =
            new task_type(std::forward<_Func>(func), std::forward<_Args>(args)...);
        std::future<_Ret> _f = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Func>
    std::future<_Ret> async(_Func&& func)
    {
        typedef PackagedTask<_Ret> task_type;
        typedef task_type*         task_pointer;

        task_pointer      _ptask = new task_type(std::forward<_Func>(func));
        std::future<_Ret> _f     = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//
    template <typename _Func, typename... _Args,
              typename _Ret = typename std::result_of<_Func(_Args&&...)>::type>
    auto async(_Func&& func, _Args&&... args) -> std::future<_Ret>
    {
        typedef PackagedTask<_Ret, _Args...> task_type;

        auto _ptask =
            new task_type(std::forward<_Func>(func), std::forward<_Args>(args)...);
        auto _f = _ptask->get_future();
        m_pool->add_task(_ptask);
        return _f;
    }
    //------------------------------------------------------------------------//

public:
    //------------------------------------------------------------------------//
    // public wrap functions
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func, typename... _Args>
    Task<_Ret, _Arg, _Args...>* wrap(TaskGroup<_Ret, _Arg>& tg, _Func&& func,
                                     _Args&&... args)
    {
        return tg.wrap(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func>
    Task<_Ret, _Arg>* wrap(TaskGroup<_Ret, _Arg>& tg, _Func&& func)
    {
        return tg.wrap(std::forward<_Func>(func));
    }

public:
    //------------------------------------------------------------------------//
    // public exec functions
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func, typename... _Args>
    void exec(TaskGroup<_Ret, _Arg>& tg, _Func&& func, _Args&&... args)
    {
        tg.exec(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func>
    void exec(TaskGroup<_Ret, _Arg>& tg, _Func&& func)
    {
        tg.exec(std::forward<_Func>(func));
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func, typename... _Args>
    void rexec(TaskGroup<_Ret, _Arg>& tg, _Func&& func, _Args&&... args)
    {
        tg.exec(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func>
    void rexec(TaskGroup<_Ret, _Arg>& tg, _Func&& func)
    {
        tg.exec(std::forward<_Func>(func));
    }
    //------------------------------------------------------------------------//
    // public exec functions (void specializations)
    //------------------------------------------------------------------------//
    template <typename _Func, typename... _Args>
    void rexec(TaskGroup<void, void>& tg, _Func&& func, _Args&&... args)
    {
        tg.exec(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Func>
    void rexec(TaskGroup<void, void>& tg, _Func&& func)
    {
        tg.exec(std::forward<_Func>(func));
    }
    //------------------------------------------------------------------------//

#if defined(PTL_USE_TBB)
    //------------------------------------------------------------------------//
    // public wrap functions using TBB tasks
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func, typename... _Args>
    Task<_Ret, _Arg, _Args...>* wrap(TBBTaskGroup<_Ret, _Arg>& tg, _Func&& func,
                                     _Args&&... args)
    {
        return tg.wrap(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func>
    Task<_Ret, _Arg>* wrap(TBBTaskGroup<_Ret, _Arg>& tg, _Func&& func)
    {
        return tg.wrap(std::forward<_Func>(func));
    }

    //------------------------------------------------------------------------//
    // public exec functions using TBB tasks
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func, typename... _Args>
    void exec(TBBTaskGroup<_Ret, _Arg>& tg, _Func&& func, _Args&&... args)
    {
        tg.exec(std::forward<_Func>(func), std::forward<_Args>(args)...);
    }
    //------------------------------------------------------------------------//
    template <typename _Ret, typename _Arg, typename _Func>
    void exec(TBBTaskGroup<_Ret, _Arg>& tg, _Func&& func)
    {
        tg.exec(std::forward<_Func>(func));
    }
    //------------------------------------------------------------------------//
#endif

protected:
    // Protected variables
    ThreadPool* m_pool;

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
