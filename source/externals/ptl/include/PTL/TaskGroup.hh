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
// This file creates the a class for handling a group of tasks that
// can be independently joined
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/AutoLock.hh"
#ifndef G4GMAKE
#include "PTL/Config.hh"
#endif
#include "PTL/Globals.hh"
#include "PTL/JoinFunction.hh"
#include "PTL/Task.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Threading.hh"
#include "PTL/Utility.hh"
#include "PTL/VTask.hh"
#include "PTL/VUserTaskQueue.hh"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(PTL_USE_TBB)
#    include <tbb/task_group.h>
#endif

namespace PTL
{
namespace internal
{
std::atomic_uintmax_t&
task_group_counter();

ThreadPool*
get_default_threadpool();

intmax_t
get_task_depth();
}  // namespace internal

template <typename Tp, typename Arg = Tp, intmax_t MaxDepth = 0>
class TaskGroup
{
public:
    //------------------------------------------------------------------------//
    template <typename Up>
    using container_type = std::vector<Up>;

    using tid_type               = std::thread::id;
    using size_type              = uintmax_t;
    using lock_t                 = Mutex;
    using atomic_int             = std::atomic_intmax_t;
    using atomic_uint            = std::atomic_uintmax_t;
    using condition_t            = Condition;
    using ArgTp                  = decay_t<Arg>;
    using result_type            = Tp;
    using task_pointer           = std::shared_ptr<TaskFuture<ArgTp>>;
    using task_list_t            = container_type<task_pointer>;
    using this_type              = TaskGroup<Tp, Arg, MaxDepth>;
    using promise_type           = std::promise<ArgTp>;
    using future_type            = std::future<ArgTp>;
    using packaged_task_type     = std::packaged_task<ArgTp()>;
    using future_list_t          = container_type<future_type>;
    using join_type              = typename JoinFunction<Tp, Arg>::Type;
    using iterator               = typename future_list_t::iterator;
    using reverse_iterator       = typename future_list_t::reverse_iterator;
    using const_iterator         = typename future_list_t::const_iterator;
    using const_reverse_iterator = typename future_list_t::const_reverse_iterator;
    //------------------------------------------------------------------------//
    template <typename... Args>
    using task_type = Task<ArgTp, decay_t<Args>...>;
    //------------------------------------------------------------------------//

public:
    // Constructor
    template <typename Func>
    TaskGroup(Func&& _join, ThreadPool* _tp = internal::get_default_threadpool());

    template <typename Up = Tp>
    TaskGroup(ThreadPool* _tp = internal::get_default_threadpool(),
              enable_if_t<std::is_void<Up>::value, int> = 0);

    // Destructor
    ~TaskGroup();

    // delete copy-construct
    TaskGroup(const this_type&) = delete;
    // define move-construct
    // NOLINTNEXTLINE(performance-noexcept-move-constructor)
    TaskGroup(this_type&& rhs) = default;
    // delete copy-assign
    TaskGroup& operator=(const this_type& rhs) = delete;
    // define move-assign
    // NOLINTNEXTLINE(performance-noexcept-move-constructor)
    TaskGroup& operator=(this_type&& rhs) = default;

public:
    template <typename Up>
    std::shared_ptr<Up> operator+=(std::shared_ptr<Up>&& _task);

    // wait to finish
    void wait();

    // increment (prefix)
    intmax_t operator++() { return ++(m_tot_task_count); }
    intmax_t operator++(int) { return (m_tot_task_count)++; }
    intmax_t operator--() { return --(m_tot_task_count); }
    intmax_t operator--(int) { return (m_tot_task_count)--; }

    // size
    intmax_t size() const { return m_tot_task_count.load(); }

    // get the locks/conditions
    lock_t&      task_lock() { return m_task_lock; }
    condition_t& task_cond() { return m_task_cond; }

    // identifier
    uintmax_t id() const { return m_id; }

    // thread pool
    void         set_pool(ThreadPool* tp) { m_pool = tp; }
    ThreadPool*& pool() { return m_pool; }
    ThreadPool*  pool() const { return m_pool; }

    bool is_native_task_group() const { return (m_tbb_task_group) == nullptr; }
    bool is_main() const { return this_tid() == m_main_tid; }

    // check if any tasks are still pending
    intmax_t pending() { return m_tot_task_count.load(); }

    static void set_verbose(int level) { f_verbose = level; }

    ScopeDestructor get_scope_destructor();

    void notify();
    void notify_all();

    void reserve(size_t _n)
    {
        m_task_list.reserve(_n);
        m_future_list.reserve(_n);
    }

public:
    template <typename Func, typename... Args>
    std::shared_ptr<task_type<Args...>> wrap(Func func, Args... args)
    {
        return operator+=(std::make_shared<task_type<Args...>>(
            is_native_task_group(), m_depth, std::move(func), std::move(args)...));
    }

    template <typename Func, typename... Args, typename Up = Tp>
    enable_if_t<std::is_void<Up>::value, void> exec(Func func, Args... args);

    template <typename Func, typename... Args, typename Up = Tp>
    enable_if_t<!std::is_void<Up>::value, void> exec(Func func, Args... args);

    template <typename Func, typename... Args>
    void run(Func func, Args... args)
    {
        exec(std::move(func), std::move(args)...);
    }

protected:
    template <typename Up, typename Func, typename... Args>
    enable_if_t<std::is_void<Up>::value, void> local_exec(Func func, Args... args);

    template <typename Up, typename Func, typename... Args>
    enable_if_t<!std::is_void<Up>::value, void> local_exec(Func func, Args... args);

    // shorter typedefs
    using itr_t   = iterator;
    using citr_t  = const_iterator;
    using ritr_t  = reverse_iterator;
    using critr_t = const_reverse_iterator;

public:
    //------------------------------------------------------------------------//
    // Get tasks with non-void return types
    //
    future_list_t&       get_tasks() { return m_future_list; }
    const future_list_t& get_tasks() const { return m_future_list; }

    //------------------------------------------------------------------------//
    // iterate over tasks with return type
    //
    itr_t   begin() { return m_future_list.begin(); }
    itr_t   end() { return m_future_list.end(); }
    citr_t  begin() const { return m_future_list.begin(); }
    citr_t  end() const { return m_future_list.end(); }
    citr_t  cbegin() const { return m_future_list.begin(); }
    citr_t  cend() const { return m_future_list.end(); }
    ritr_t  rbegin() { return m_future_list.rbegin(); }
    ritr_t  rend() { return m_future_list.rend(); }
    critr_t rbegin() const { return m_future_list.rbegin(); }
    critr_t rend() const { return m_future_list.rend(); }

    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, enable_if_t<!std::is_void<Up>::value, int> = 0>
    inline Up join(Up accum = {});
    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, typename Rp = Arg,
              enable_if_t<std::is_void<Up>::value && std::is_void<Rp>::value, int> = 0>
    inline void join();
    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, typename Rp = Arg,
              enable_if_t<std::is_void<Up>::value && !std::is_void<Rp>::value, int> = 0>
    inline void join();
    //------------------------------------------------------------------------//
    // clear the task result history
    void clear();

protected:
    //------------------------------------------------------------------------//
    // get the thread id
    static tid_type this_tid() { return std::this_thread::get_id(); }

    //------------------------------------------------------------------------//
    // get the task count
    atomic_int&       task_count() { return m_tot_task_count; }
    const atomic_int& task_count() const { return m_tot_task_count; }

protected:
    static int f_verbose;
    // Private variables
    uintmax_t         m_id       = internal::task_group_counter()++;
    intmax_t          m_depth    = internal::get_task_depth();
    tid_type          m_main_tid = std::this_thread::get_id();
    atomic_int        m_tot_task_count{ 0 };
    lock_t            m_task_lock = {};
    condition_t       m_task_cond = {};
    join_type         m_join{};
    ThreadPool*       m_pool           = internal::get_default_threadpool();
    tbb_task_group_t* m_tbb_task_group = nullptr;
    task_list_t       m_task_list      = {};
    future_list_t     m_future_list    = {};

private:
    void internal_update();
};

}  // namespace PTL
namespace PTL
{
template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Func>
TaskGroup<Tp, Arg, MaxDepth>::TaskGroup(Func&& _join, ThreadPool* _tp)
: m_join{ std::forward<Func>(_join) }
, m_pool{ _tp }
{
    internal_update();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up>
TaskGroup<Tp, Arg, MaxDepth>::TaskGroup(ThreadPool* _tp,
                                        enable_if_t<std::is_void<Up>::value, int>)
: m_join{ []() {} }
, m_pool{ _tp }
{
    internal_update();
}

// Destructor
template <typename Tp, typename Arg, intmax_t MaxDepth>
TaskGroup<Tp, Arg, MaxDepth>::~TaskGroup()
{
    {
        // task will decrement counter and then acquire the lock to notify
        // condition variable so acquiring lock here will prevent the
        // task group from being destroyed before this is completed
        AutoLock _lk{ m_task_lock, std::defer_lock };
        if(!_lk.owns_lock())
            _lk.lock();
    }

    if(m_tbb_task_group)
    {
        auto* _arena = m_pool->get_task_arena();
        _arena->execute([this]() { this->m_tbb_task_group->wait(); });
    }
    delete m_tbb_task_group;
    this->clear();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up>
std::shared_ptr<Up>
TaskGroup<Tp, Arg, MaxDepth>::operator+=(std::shared_ptr<Up>&& _task)
{
    // thread-safe increment of tasks in task group
    operator++();
    // copy the shared pointer to abstract instance
    m_task_list.push_back(_task);
    // return the derived instance
    return std::move(_task);
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
void
TaskGroup<Tp, Arg, MaxDepth>::wait()
{
    auto _dtor = ScopeDestructor{ [&]() {
        if(m_tbb_task_group)
        {
            auto* _arena = m_pool->get_task_arena();
            _arena->execute([this]() { this->m_tbb_task_group->wait(); });
        }
    } };

    ThreadData* data = ThreadData::GetInstance();
    if(!data)
        return;

    // if no pool was initially present at creation
    if(!m_pool)
    {
        // check for master MT run-manager
        m_pool = internal::get_default_threadpool();

        // if no thread pool created
        if(!m_pool)
        {
            if(f_verbose > 0)
            {
                fprintf(stderr, "%s @ %i :: Warning! nullptr to thread-pool (%p)\n",
                        __FUNCTION__, __LINE__, static_cast<void*>(m_pool));
                std::cerr << __FUNCTION__ << "@" << __LINE__ << " :: Warning! "
                          << "nullptr to thread pool!" << std::endl;
            }
            return;
        }
    }

    ThreadPool*     tpool = (m_pool) ? m_pool : data->thread_pool;
    VUserTaskQueue* taskq = (tpool) ? tpool->get_queue() : data->current_queue;

    bool _is_main     = data->is_main;
    bool _within_task = data->within_task;

    auto is_active_state = [&]() {
        return (tpool->state()->load(std::memory_order_relaxed) !=
                thread_pool::state::STOPPED);
    };

    auto execute_this_threads_tasks = [&]() {
        if(!taskq)
            return;

        // only want to process if within a task
        if((!_is_main || tpool->size() < 2) && _within_task)
        {
            int bin = static_cast<int>(taskq->GetThreadBin());
            // const auto nitr = (tpool) ? tpool->size() :
            // Thread::hardware_concurrency();
            while(this->pending() > 0)
            {
                if(!taskq->empty())
                {
                    auto _task = taskq->GetTask(bin);
                    if(_task)
                        (*_task)();
                }
            }
        }
    };

    // checks for validity
    if(!is_native_task_group())
    {
        // for external threads
        if(!_is_main || tpool->size() < 2)
            return;
    }
    else if(f_verbose > 0)
    {
        if(!tpool || !taskq)
        {
            // something is wrong, didn't create thread-pool?
            fprintf(stderr,
                    "%s @ %i :: Warning! nullptr to thread data (%p) or task-queue "
                    "(%p)\n",
                    __FUNCTION__, __LINE__, static_cast<void*>(tpool),
                    static_cast<void*>(taskq));
        }
        // return if thread pool isn't built
        else if(is_native_task_group() && !tpool->is_alive())
        {
            fprintf(stderr, "%s @ %i :: Warning! thread-pool is not alive!\n",
                    __FUNCTION__, __LINE__);
        }
        else if(!is_active_state())
        {
            fprintf(stderr, "%s @ %i :: Warning! thread-pool is not active!\n",
                    __FUNCTION__, __LINE__);
        }
    }

    intmax_t wake_size = 2;
    AutoLock _lock(m_task_lock, std::defer_lock);

    while(is_active_state())
    {
        execute_this_threads_tasks();

        // while loop protects against spurious wake-ups
        while(_is_main && pending() > 0 && is_active_state())
        {
            // auto _wake = [&]() { return (wake_size > pending() ||
            // !is_active_state());
            // };

            // lock before sleeping on condition
            if(!_lock.owns_lock())
                _lock.lock();

            // Wait until signaled that a task has been competed
            // Unlock mutex while wait, then lock it back when signaled
            // when true, this wakes the thread
            if(pending() >= wake_size)
            {
                m_task_cond.wait(_lock);
            }
            else
            {
                m_task_cond.wait_for(_lock, std::chrono::microseconds(100));
            }
            // unlock
            if(_lock.owns_lock())
                _lock.unlock();
        }

        // if pending is not greater than zero, we are joined
        if(pending() <= 0)
            break;
    }

    if(_lock.owns_lock())
        _lock.unlock();

    intmax_t ntask = this->task_count().load();
    if(ntask > 0)
    {
        std::stringstream ss;
        ss << "\nWarning! Join operation issue! " << ntask << " tasks still "
           << "are running!" << std::endl;
        std::cerr << ss.str();
        this->wait();
    }
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
ScopeDestructor
TaskGroup<Tp, Arg, MaxDepth>::get_scope_destructor()
{
    auto& _counter   = m_tot_task_count;
    auto& _task_cond = task_cond();
    auto& _task_lock = task_lock();
    return ScopeDestructor{ [&_task_cond, &_task_lock, &_counter]() {
        auto _count = --(_counter);
        if(_count < 1)
        {
            AutoLock _lk{ _task_lock };
            _task_cond.notify_all();
        }
    } };
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
void
TaskGroup<Tp, Arg, MaxDepth>::notify()
{
    AutoLock _lk{ m_task_lock };
    m_task_cond.notify_one();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
void
TaskGroup<Tp, Arg, MaxDepth>::notify_all()
{
    AutoLock _lk{ m_task_lock };
    m_task_cond.notify_all();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Func, typename... Args, typename Up>
enable_if_t<std::is_void<Up>::value, void>
TaskGroup<Tp, Arg, MaxDepth>::exec(Func func, Args... args)
{
    if(MaxDepth > 0 && !m_tbb_task_group && ThreadData::GetInstance() &&
       ThreadData::GetInstance()->task_depth > MaxDepth)
    {
        local_exec<Tp>(std::move(func), std::move(args)...);
    }
    else
    {
        auto& _counter   = m_tot_task_count;
        auto& _task_cond = task_cond();
        auto& _task_lock = task_lock();
        auto  _task      = wrap([&_task_cond, &_task_lock, &_counter, func, args...]() {
            auto* _tdata = ThreadData::GetInstance();
            if(_tdata)
                ++(_tdata->task_depth);
            func(args...);
            auto _count = --(_counter);
            if(_tdata)
                --(_tdata->task_depth);
            if(_count < 1)
            {
                AutoLock _lk{ _task_lock };
                _task_cond.notify_all();
            }
        });

        if(m_tbb_task_group)
        {
            auto* _arena          = m_pool->get_task_arena();
            auto* _tbb_task_group = m_tbb_task_group;
            auto* _ptask          = _task.get();
            _arena->execute([_tbb_task_group, _ptask]() {
                _tbb_task_group->run([_ptask]() { (*_ptask)(); });
            });
        }
        else
        {
            m_pool->add_task(std::move(_task));
        }
    }
}
template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Func, typename... Args, typename Up>
enable_if_t<!std::is_void<Up>::value, void>
TaskGroup<Tp, Arg, MaxDepth>::exec(Func func, Args... args)
{
    if(MaxDepth > 0 && !m_tbb_task_group && ThreadData::GetInstance() &&
       ThreadData::GetInstance()->task_depth > MaxDepth)
    {
        local_exec<Tp>(std::move(func), std::move(args)...);
    }
    else
    {
        auto& _counter   = m_tot_task_count;
        auto& _task_cond = task_cond();
        auto& _task_lock = task_lock();
        auto  _task      = wrap([&_task_cond, &_task_lock, &_counter, func, args...]() {
            auto* _tdata = ThreadData::GetInstance();
            if(_tdata)
                ++(_tdata->task_depth);
            auto&& _ret   = func(args...);
            auto   _count = --(_counter);
            if(_tdata)
                --(_tdata->task_depth);
            if(_count < 1)
            {
                AutoLock _lk{ _task_lock };
                _task_cond.notify_all();
            }
            return std::forward<decltype(_ret)>(_ret);
        });

        if(m_tbb_task_group)
        {
            auto* _arena          = m_pool->get_task_arena();
            auto* _tbb_task_group = m_tbb_task_group;
            auto* _ptask          = _task.get();
            _arena->execute([_tbb_task_group, _ptask]() {
                _tbb_task_group->run([_ptask]() { (*_ptask)(); });
            });
        }
        else
        {
            m_pool->add_task(std::move(_task));
        }
    }
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up, typename Func, typename... Args>
enable_if_t<std::is_void<Up>::value, void>
TaskGroup<Tp, Arg, MaxDepth>::local_exec(Func func, Args... args)
{
    auto* _tdata = ThreadData::GetInstance();
    if(_tdata)
        ++(_tdata->task_depth);
    promise_type _p{};
    m_future_list.emplace_back(_p.get_future());
    func(args...);
    _p.set_value();
    if(_tdata)
        --(_tdata->task_depth);
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up, typename Func, typename... Args>
enable_if_t<!std::is_void<Up>::value, void>
TaskGroup<Tp, Arg, MaxDepth>::local_exec(Func func, Args... args)
{
    auto* _tdata = ThreadData::GetInstance();
    if(_tdata)
        ++(_tdata->task_depth);
    promise_type _p{};
    m_future_list.emplace_back(_p.get_future());
    _p.set_value(func(args...));
    if(_tdata)
        --(_tdata->task_depth);
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up, enable_if_t<!std::is_void<Up>::value, int>>
inline Up
TaskGroup<Tp, Arg, MaxDepth>::join(Up accum)
{
    this->wait();
    for(auto& itr : m_task_list)
    {
        using RetT = decay_t<decltype(itr->get())>;
        accum      = std::move(m_join(std::ref(accum), std::forward<RetT>(itr->get())));
    }
    for(auto& itr : m_future_list)
    {
        using RetT = decay_t<decltype(itr.get())>;
        accum      = std::move(m_join(std::ref(accum), std::forward<RetT>(itr.get())));
    }
    this->clear();
    return accum;
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up, typename Rp,
          enable_if_t<std::is_void<Up>::value && std::is_void<Rp>::value, int>>
inline void
TaskGroup<Tp, Arg, MaxDepth>::join()
{
    this->wait();
    for(auto& itr : m_task_list)
        itr->get();
    for(auto& itr : m_future_list)
        itr.get();
    m_join();
    this->clear();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
template <typename Up, typename Rp,
          enable_if_t<std::is_void<Up>::value && !std::is_void<Rp>::value, int>>
inline void
TaskGroup<Tp, Arg, MaxDepth>::join()
{
    this->wait();
    for(auto& itr : m_task_list)
    {
        using RetT = decay_t<decltype(itr->get())>;
        m_join(std::forward<RetT>(itr->get()));
    }
    for(auto& itr : m_future_list)
    {
        using RetT = decay_t<decltype(itr.get())>;
        m_join(std::forward<RetT>(itr.get()));
    }
    this->clear();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
void
TaskGroup<Tp, Arg, MaxDepth>::clear()
{
    m_future_list.clear();
    m_task_list.clear();
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
void
TaskGroup<Tp, Arg, MaxDepth>::internal_update()
{
    if(!m_pool)
        m_pool = internal::get_default_threadpool();

    if(!m_pool)
    {
        std::stringstream ss{};
        ss << "[TaskGroup]> " << __FUNCTION__ << "@" << __LINE__
           << " :: nullptr to thread pool";
        throw std::runtime_error(ss.str());
    }

    if(m_pool->is_tbb_threadpool())
    {
        m_tbb_task_group = new tbb_task_group_t{};
    }
}

template <typename Tp, typename Arg, intmax_t MaxDepth>
int TaskGroup<Tp, Arg, MaxDepth>::f_verbose = GetEnv<int>("PTL_VERBOSE", 0);

}  // namespace PTL
