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

#include "PTL/JoinFunction.hh"
#include "PTL/Task.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Utility.hh"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <future>
#include <iostream>
#include <memory>
#include <vector>

#if defined(PTL_USE_TBB)
#    include <tbb/tbb.h>
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
    TaskGroup(this_type&& rhs) = default;
    // delete copy-assign
    TaskGroup& operator=(const this_type& rhs) = delete;
    // define move-assign
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

    bool is_native_task_group() const { return (m_tbb_task_group) ? false : true; }
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

#include "PTL/TaskGroup.icc"
