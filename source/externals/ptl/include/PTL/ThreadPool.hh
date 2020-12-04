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
// This file creates a class for an efficient thread-pool that
// accepts work in the form of tasks.
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/AutoLock.hh"
#include "PTL/ThreadData.hh"
#include "PTL/Threading.hh"
#include "PTL/VTask.hh"
#include "PTL/VTaskGroup.hh"
#include "PTL/VUserTaskQueue.hh"

#ifdef PTL_USE_TBB
#    include <tbb/global_control.h>
#    include <tbb/tbb.h>
#endif

// C
#include <cstdint>
#include <cstdlib>
#include <cstring>
// C++
#include <atomic>
#include <deque>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <vector>

namespace PTL
{
class ThreadPool
{
public:
    template <typename KeyT, typename MappedT, typename HashT = KeyT>
    using uomap = std::unordered_map<KeyT, MappedT, std::hash<HashT>>;

    // pod-types
    using size_type        = size_t;
    using task_count_type  = std::shared_ptr<std::atomic_uintmax_t>;
    using atomic_int_type  = std::shared_ptr<std::atomic_uintmax_t>;
    using pool_state_type  = std::shared_ptr<std::atomic_short>;
    using atomic_bool_type = std::shared_ptr<std::atomic_bool>;
    // objects
    using task_type    = VTask;
    using lock_t       = std::shared_ptr<Mutex>;
    using condition_t  = std::shared_ptr<Condition>;
    using task_pointer = task_type*;
    using task_queue_t = VUserTaskQueue;
    // containers
    typedef std::deque<ThreadId>          thread_list_t;
    typedef std::vector<bool>             bool_list_t;
    typedef std::map<ThreadId, uintmax_t> thread_id_map_t;
    typedef std::map<uintmax_t, ThreadId> thread_index_map_t;
    using thread_vec_t = std::vector<Thread>;
    // functions
    typedef std::function<void()>             initialize_func_t;
    typedef std::function<intmax_t(intmax_t)> affinity_func_t;

public:
    // Constructor and Destructors
    ThreadPool(
        const size_type& pool_size, VUserTaskQueue* task_queue = nullptr,
        bool _use_affinity     = GetEnv<bool>("PTL_CPU_AFFINITY", false),
        const affinity_func_t& = [](intmax_t) {
            static std::atomic<intmax_t> assigned;
            intmax_t                     _assign = assigned++;
            return _assign % Thread::hardware_concurrency();
        });
    // Virtual destructors are required by abstract classes
    // so add it by default, just in case
    virtual ~ThreadPool();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&)      = default;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&) = default;

public:
    // Public functions
    size_type initialize_threadpool(size_type);  // start the threads
    size_type destroy_threadpool();              // destroy the threads
    size_type stop_thread();

    template <typename FuncT>
    void execute_on_all_threads(FuncT&& _func);

public:
    // Public functions related to TBB
    static bool using_tbb();
    // enable using TBB if available
    static void set_use_tbb(bool val);

public:
    // add tasks for threads to process
    size_type add_task(task_pointer task, int bin = -1);
    // size_type add_thread_task(ThreadId id, task_pointer&& task);
    // add a generic container with iterator
    template <typename ListT>
    size_type add_tasks(ListT&);

    Thread* get_thread(size_type _n) const;
    Thread* get_thread(std::thread::id id) const;

    task_queue_t* get_queue() const { return m_task_queue; }

    // only relevant when compiled with PTL_USE_TBB
    static tbb_global_control_t*& tbb_global_control();

    void set_initialization(initialize_func_t f) { m_init_func = f; }
    void reset_initialization()
    {
        auto f      = []() {};
        m_init_func = f;
    }

public:
    // get the pool state
    const pool_state_type& state() const { return m_pool_state; }
    // see how many main task threads there are
    size_type size() const { return m_pool_size; }
    // set the thread pool size
    void resize(size_type _n);
    // affinity assigns threads to cores, assignment at constructor
    bool using_affinity() const { return m_use_affinity; }
    bool is_alive() { return m_alive_flag->load(); }
    void notify();
    void notify_all();
    void notify(size_type);
    bool is_initialized() const;
    int  get_active_threads_count() const
    {
        return (m_thread_awake) ? m_thread_awake->load() : 0;
    }

    void set_affinity(affinity_func_t f) { m_affinity_func = f; }
    void set_affinity(intmax_t i, Thread&);

    void set_verbose(int n) { m_verbose = n; }
    int  get_verbose() const { return m_verbose; }
    bool is_master() const { return ThisThread::get_id() == m_master_tid; }

public:
    // read FORCE_NUM_THREADS environment variable
    static const thread_id_map_t& get_thread_ids();
    static uintmax_t              get_this_thread_id();

protected:
    void execute_thread(VUserTaskQueue*);  // function thread sits in
    int  insert(const task_pointer&, int = -1);
    int  run_on_this(task_pointer);

protected:
    // called in THREAD INIT
    static void start_thread(ThreadPool*, intmax_t = -1);

    void record_entry()
    {
        if(m_thread_active)
            ++(*m_thread_active);
    }

    void record_exit()
    {
        if(m_thread_active)
            --(*m_thread_active);
    }

private:
    // Private variables
    // random
    bool             m_use_affinity;
    bool             m_tbb_tp;
    int              m_verbose   = 0;
    size_type        m_pool_size = 0;
    ThreadId         m_master_tid;
    atomic_bool_type m_alive_flag   = std::make_shared<std::atomic_bool>(false);
    pool_state_type  m_pool_state   = std::make_shared<std::atomic_short>(0);
    atomic_int_type  m_thread_awake = std::make_shared<std::atomic_uintmax_t>();
    atomic_int_type  m_thread_active = std::make_shared<std::atomic_uintmax_t>();

    // locks
    lock_t m_task_lock = std::make_shared<Mutex>();
    // conditions
    condition_t m_task_cond = std::make_shared<Condition>();

    // containers
    bool_list_t   m_is_joined;     // join list
    bool_list_t   m_is_stopped;    // lets thread know to stop
    thread_list_t m_main_threads;  // storage for active threads
    thread_list_t m_stop_threads;  // storage for stopped threads
    thread_vec_t  m_threads;

    // task queue
    task_queue_t*     m_task_queue;
    tbb_task_group_t* m_tbb_task_group;

    // functions
    initialize_func_t m_init_func;
    affinity_func_t   m_affinity_func;

private:
    // Private static variables
    PTL_DLL static thread_id_map_t f_thread_ids;
    PTL_DLL static bool            f_use_tbb;
};

//--------------------------------------------------------------------------------------//
inline void
ThreadPool::notify()
{
    // wake up one thread that is waiting for a task to be available
    if(m_thread_awake && m_thread_awake->load() < m_pool_size)
    {
        AutoLock l(*m_task_lock);
        m_task_cond->notify_one();
    }
}
//--------------------------------------------------------------------------------------//
inline void
ThreadPool::notify_all()
{
    // wake all threads
    AutoLock l(*m_task_lock);
    m_task_cond->notify_all();
}
//--------------------------------------------------------------------------------------//
inline void
ThreadPool::notify(size_type ntasks)
{
    if(ntasks == 0)
        return;

    // wake up as many threads that tasks just added
    if(m_thread_awake && m_thread_awake->load() < m_pool_size)
    {
        AutoLock l(*m_task_lock);
        if(ntasks < this->size())
        {
            for(size_type i = 0; i < ntasks; ++i)
                m_task_cond->notify_one();
        }
        else
            m_task_cond->notify_all();
    }
}
//--------------------------------------------------------------------------------------//
// local function for getting the tbb task scheduler
inline tbb_global_control_t*&
ThreadPool::tbb_global_control()
{
    static thread_local tbb_global_control_t* _instance = nullptr;
    return _instance;
}
//--------------------------------------------------------------------------------------//
inline void
ThreadPool::resize(size_type _n)
{
    if(_n == m_pool_size)
        return;
    initialize_threadpool(_n);
    m_task_queue->resize(static_cast<intmax_t>(_n));
}
//--------------------------------------------------------------------------------------//
inline int
ThreadPool::run_on_this(task_pointer task)
{
    auto _func = [=]() {
        (*task)();
        if(!task->group())
            delete task;
    };

    if(m_tbb_tp && m_tbb_task_group)
    {
        m_tbb_task_group->run(_func);
    }
    else
    {
        _func();
    }
    // return the number of tasks added to task-list
    return 0;
}
//--------------------------------------------------------------------------------------//
inline int
ThreadPool::insert(const task_pointer& task, int bin)
{
    static thread_local ThreadData* _data = ThreadData::GetInstance();

    // pass the task to the queue
    auto ibin = m_task_queue->InsertTask(task, _data, bin);
    notify();
    return ibin;
}
//--------------------------------------------------------------------------------------//
inline ThreadPool::size_type
ThreadPool::add_task(task_pointer task, int bin)
{
    // if not native (i.e. TBB) then return
    if(!task->is_native_task())
        return 0;

    // if we haven't built thread-pool, just execute
    if(!m_alive_flag->load())
        return static_cast<size_type>(run_on_this(task));

    return static_cast<size_type>(insert(task, bin));
}
//--------------------------------------------------------------------------------------//
template <typename ListT>
inline ThreadPool::size_type
ThreadPool::add_tasks(ListT& c)
{
    if(!m_alive_flag)  // if we haven't built thread-pool, just execute
    {
        for(auto& itr : c)
            run(itr);
        c.clear();
        return 0;
    }

    // TODO: put a limit on how many tasks can be added at most
    auto c_size = c.size();
    for(auto& itr : c)
    {
        if(!itr->is_native_task())
            --c_size;
        else
        {
            //++(m_task_queue);
            m_task_queue->InsertTask(itr);
        }
    }
    c.clear();

    // notify sleeping threads
    notify(c_size);

    return c_size;
}
//--------------------------------------------------------------------------------------//
template <typename FuncT>
inline void
ThreadPool::execute_on_all_threads(FuncT&& _func)
{
    if(m_tbb_tp && m_tbb_task_group)
    {
#if defined(PTL_USE_TBB)
        // TBB lazily activates threads to process tasks and the master thread
        // participates in processing the tasks so getting a specific
        // function to execute only on the worker threads requires some trickery
        //
        auto                      master_tid = ThisThread::get_id();
        std::set<std::thread::id> _first;
        Mutex                     _mutex;
        // init function which executes function and returns 1 only once
        auto _init = [&]() {
            static thread_local int _once = 0;
            _mutex.lock();
            if(_first.find(std::this_thread::get_id()) == _first.end())
            {
                // we need to reset this thread-local static for multiple invocations
                // of the same template instantiation
                _once = 0;
                _first.insert(std::this_thread::get_id());
            }
            _mutex.unlock();
            if(_once++ == 0)
            {
                _func();
                return 1;
            }
            return 0;
        };
        // consumes approximately N milliseconds of cpu time
        auto _consume = [](long n) {
            using stl_mutex_t   = std::mutex;
            using unique_lock_t = std::unique_lock<stl_mutex_t>;
            // a mutex held by one lock
            stl_mutex_t mutex;
            // acquire lock
            unique_lock_t hold_lk(mutex);
            // associate but defer
            unique_lock_t try_lk(mutex, std::defer_lock);
            // get current time
            auto now = std::chrono::steady_clock::now();
            // try until time point
            while(std::chrono::steady_clock::now() < (now + std::chrono::milliseconds(n)))
                try_lk.try_lock();
        };
        // this will collect the number of threads which have
        // executed the _init function above
        std::atomic<size_t> _total_init{ 0 };
        // this is the task passed to the task-group
        auto _init_task = [&]() {
            int _ret = 0;
            // don't let the master thread execute the function
            if(ThisThread::get_id() != master_tid)
            {
                // execute the function
                _ret = _init();
                // add the result
                _total_init += _ret;
            }
            // if the function did not return anything, put it to sleep
            // so TBB will wake other threads to execute the remaining tasks
            if(_ret == 0)
                _consume(100);
        };

        // TBB won't oversubscribe so we need to limit by ncores - 1
        size_t nitr  = 0;
        size_t _maxp = tbb_global_control()->active_value(
            tbb::global_control::max_allowed_parallelism);
        size_t _sz         = size();
        size_t _ncore      = Threading::GetNumberOfCores() - 1;
        size_t _num        = std::min(_maxp, std::min(_sz, _ncore));
        auto   _fname      = __FUNCTION__;
        auto   _write_info = [&]() {
            std::cerr << "[" << _fname << "]> Total initalized: " << _total_init
                      << ", expected: " << _num << ", max-parallel: " << _maxp
                      << ", size: " << _sz << ", ncore: " << _ncore << std::endl;
        };
        while(_total_init < _num)
        {
            auto _n = _num;
            while(--_n > 0)
                m_tbb_task_group->run(_init_task);
            m_tbb_task_group->wait();
            // don't loop infinitely but use a strict condition
            if(nitr++ > 2 * (_num + 1) && (_total_init - 1) == _num)
            {
                _write_info();
                break;
            }
            // at this point we need to exit
            if(nitr > 4 * (_ncore + 1))
            {
                _write_info();
                break;
            }
        }
        if(get_verbose() > 3)
            _write_info();
#endif
    }
    else if(get_queue())
    {
        get_queue()->ExecuteOnAllThreads(this, std::forward<FuncT>(_func));
    }
}
//======================================================================================//

}  // namespace PTL
