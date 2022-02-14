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
#include "PTL/VUserTaskQueue.hh"

#if defined(PTL_USE_TBB)
#    if !defined(TBB_SUPPRESS_DEPRECATED_MESSAGES)
#        define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#    endif
#    if !defined(TBB_PREVIEW_GLOBAL_CONTROL)
#        define TBB_PREVIEW_GLOBAL_CONTROL 1
#    endif
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
    using task_pointer = std::shared_ptr<task_type>;
    using task_queue_t = VUserTaskQueue;
    // containers
    typedef std::deque<ThreadId>          thread_list_t;
    typedef std::vector<bool>             bool_list_t;
    typedef std::map<ThreadId, uintmax_t> thread_id_map_t;
    typedef std::map<uintmax_t, ThreadId> thread_index_map_t;
    using thread_vec_t  = std::vector<Thread>;
    using thread_data_t = std::vector<std::shared_ptr<ThreadData>>;
    // functions
    typedef std::function<void()>             initialize_func_t;
    typedef std::function<intmax_t(intmax_t)> affinity_func_t;

public:
    // Constructor and Destructors
    ThreadPool(const size_type& pool_size, VUserTaskQueue* task_queue = nullptr,
               bool _use_affinity = GetEnv<bool>("PTL_CPU_AFFINITY", false),
               affinity_func_t    = [](intmax_t) {
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

    template <typename FuncT>
    void execute_on_specific_threads(const std::set<std::thread::id>& _tid,
                                     FuncT&&                          _func);

    task_queue_t*  get_queue() const { return m_task_queue; }
    task_queue_t*& get_valid_queue(task_queue_t*&) const;

    bool is_tbb_threadpool() const { return m_tbb_tp; }

public:
    // Public functions related to TBB
    static bool using_tbb();
    // enable using TBB if available
    static void set_use_tbb(bool val);

public:
    // add tasks for threads to process
    size_type add_task(task_pointer&& task, int bin = -1);
    // size_type add_thread_task(ThreadId id, task_pointer&& task);
    // add a generic container with iterator
    template <typename ListT>
    size_type add_tasks(ListT&);

    Thread* get_thread(size_type _n) const;
    Thread* get_thread(std::thread::id id) const;

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
    bool is_main() const { return ThisThread::get_id() == m_main_tid; }

    tbb_task_arena_t* get_task_arena();

public:
    // read FORCE_NUM_THREADS environment variable
    static const thread_id_map_t& get_thread_ids();
    static uintmax_t              get_this_thread_id();
    static uintmax_t              add_thread_id()
    {
        AutoLock lock(TypeMutex<ThreadPool>(), std::defer_lock);
        if(!lock.owns_lock())
            lock.lock();
        auto _tid = ThisThread::get_id();
        if(f_thread_ids.find(_tid) == f_thread_ids.end())
        {
            auto _idx          = f_thread_ids.size();
            f_thread_ids[_tid] = _idx;
            Threading::SetThreadId(_idx);
        }
        return f_thread_ids.at(_tid);
    }

protected:
    void execute_thread(VUserTaskQueue*);  // function thread sits in
    int  insert(task_pointer&&, int = -1);
    int  run_on_this(task_pointer&&);

protected:
    // called in THREAD INIT
    static void start_thread(ThreadPool*, thread_data_t*, intmax_t = -1);

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
    bool             m_use_affinity      = false;
    bool             m_tbb_tp            = false;
    bool             m_delete_task_queue = false;
    int              m_verbose           = GetEnv<int>("PTL_VERBOSE", 0);
    size_type        m_pool_size         = 0;
    ThreadId         m_main_tid          = ThisThread::get_id();
    atomic_bool_type m_alive_flag        = std::make_shared<std::atomic_bool>(false);
    pool_state_type  m_pool_state        = std::make_shared<std::atomic_short>(0);
    atomic_int_type  m_thread_awake      = std::make_shared<std::atomic_uintmax_t>();
    atomic_int_type  m_thread_active     = std::make_shared<std::atomic_uintmax_t>();

    // locks
    lock_t m_task_lock = std::make_shared<Mutex>();
    // conditions
    condition_t m_task_cond = std::make_shared<Condition>();

    // containers
    bool_list_t   m_is_joined{};     // join list
    bool_list_t   m_is_stopped{};    // lets thread know to stop
    thread_list_t m_main_threads{};  // storage for active threads
    thread_list_t m_stop_threads{};  // storage for stopped threads
    thread_vec_t  m_threads{};
    thread_data_t m_thread_data{};

    // task queue
    task_queue_t*     m_task_queue     = nullptr;
    tbb_task_arena_t* m_tbb_task_arena = nullptr;
    tbb_task_group_t* m_tbb_task_group = nullptr;

    // functions
    initialize_func_t m_init_func = []() {};
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
        {
            m_task_cond->notify_all();
        }
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
// task arena
inline tbb_task_arena_t*
ThreadPool::get_task_arena()
{
#if defined(PTL_USE_TBB)
    // create a task arena
    if(!m_tbb_task_arena)
    {
        auto _sz = (tbb_global_control())
                       ? tbb_global_control()->active_value(
                             tbb::global_control::max_allowed_parallelism)
                       : size();
        m_tbb_task_arena = new tbb_task_arena_t(::tbb::task_arena::attach{});
        m_tbb_task_arena->initialize(_sz, 1);
    }
#else
    if(!m_tbb_task_arena)
        m_tbb_task_arena = new tbb_task_arena_t{};
#endif
    return m_tbb_task_arena;
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
ThreadPool::run_on_this(task_pointer&& _task)
{
    auto&& _func = [_task]() { (*_task)(); };

    if(m_tbb_tp && m_tbb_task_group)
    {
        auto _arena = get_task_arena();
        _arena->execute([this, _func]() { this->m_tbb_task_group->run(_func); });
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
ThreadPool::insert(task_pointer&& task, int bin)
{
    static thread_local ThreadData* _data = ThreadData::GetInstance();

    // pass the task to the queue
    auto ibin = get_valid_queue(m_task_queue)->InsertTask(std::move(task), _data, bin);
    notify();
    return ibin;
}
//--------------------------------------------------------------------------------------//
inline ThreadPool::size_type
ThreadPool::add_task(task_pointer&& task, int bin)
{
    // if not native (i.e. TBB) or we haven't built thread-pool, just execute
    if(m_tbb_tp || !task->is_native_task() || !m_alive_flag->load())
        return static_cast<size_type>(run_on_this(std::move(task)));

    return static_cast<size_type>(insert(std::move(task), bin));
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
            get_valid_queue(m_task_queue)->InsertTask(itr);
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
        // TBB lazily activates threads to process tasks and the main thread
        // participates in processing the tasks so getting a specific
        // function to execute only on the worker threads requires some trickery
        //
        std::set<std::thread::id> _first{};
        Mutex                     _mutex{};
        // init function which executes function and returns 1 only once
        auto _init = [&]() {
            int _once = 0;
            _mutex.lock();
            if(_first.find(std::this_thread::get_id()) == _first.end())
            {
                // we need to reset this thread-local static for multiple invocations
                // of the same template instantiation
                _once = 1;
                _first.insert(std::this_thread::get_id());
            }
            _mutex.unlock();
            if(_once != 0)
            {
                _func();
                return 1;
            }
            return 0;
        };
        // this will collect the number of threads which have
        // executed the _init function above
        std::atomic<size_t> _total_init{ 0 };
        // max parallelism by TBB
        size_t _maxp = tbb_global_control()->active_value(
            tbb::global_control::max_allowed_parallelism);
        // create a task arean
        auto _arena = get_task_arena();
        // size of the thread-pool
        size_t _sz = size();
        // number of cores
        size_t _ncore = Threading::GetNumberOfCores();
        // maximum depth for recursion
        size_t _dmax = std::max<size_t>(_ncore, 4);
        // how many threads we need to initialize
        size_t _num = std::min(_maxp, std::min(_sz, _ncore));
        // this is the task passed to the task-group
        std::function<void()> _init_task;
        _init_task = [&]() {
            add_thread_id();
            static thread_local size_type _depth = 0;
            int                           _ret   = 0;
            // don't let the main thread execute the function
            if(!is_main())
            {
                // execute the function
                _ret = _init();
                // add the result
                _total_init += _ret;
            }
            // if the function did not return anything, recursively execute
            // two more tasks
            ++_depth;
            if(_ret == 0 && _depth < _dmax && _total_init.load() < _num)
            {
                tbb::task_group tg{};
                tg.run([&]() { _init_task(); });
                tg.run([&]() { _init_task(); });
                ThisThread::sleep_for(std::chrono::milliseconds{ 1 });
                tg.wait();
            }
            --_depth;
        };

        // TBB won't oversubscribe so we need to limit by ncores - 1
        size_t nitr        = 0;
        auto   _fname      = __FUNCTION__;
        auto   _write_info = [&]() {
            std::cout << "[" << _fname << "]> Total initalized: " << _total_init
                      << ", expected: " << _num << ", max-parallel: " << _maxp
                      << ", size: " << _sz << ", ncore: " << _ncore << std::endl;
        };
        while(_total_init < _num)
        {
            auto _n = 2 * _num;
            while(--_n > 0)
            {
                _arena->execute(
                    [&]() { m_tbb_task_group->run([&]() { _init_task(); }); });
            }
            _arena->execute([&]() { m_tbb_task_group->wait(); });
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

//--------------------------------------------------------------------------------------//

template <typename FuncT>
inline void
ThreadPool::execute_on_specific_threads(const std::set<std::thread::id>& _tids,
                                        FuncT&&                          _func)
{
    if(m_tbb_tp && m_tbb_task_group)
    {
#if defined(PTL_USE_TBB)
        // TBB lazily activates threads to process tasks and the main thread
        // participates in processing the tasks so getting a specific
        // function to execute only on the worker threads requires some trickery
        //
        std::set<std::thread::id> _first{};
        Mutex                     _mutex{};
        // init function which executes function and returns 1 only once
        auto _exec = [&]() {
            int _once = 0;
            _mutex.lock();
            if(_first.find(std::this_thread::get_id()) == _first.end())
            {
                // we need to reset this thread-local static for multiple invocations
                // of the same template instantiation
                _once = 1;
                _first.insert(std::this_thread::get_id());
            }
            _mutex.unlock();
            if(_once != 0)
            {
                _func();
                return 1;
            }
            return 0;
        };
        // this will collect the number of threads which have
        // executed the _exec function above
        std::atomic<size_t> _total_exec{ 0 };
        // number of cores
        size_t _ncore = Threading::GetNumberOfCores();
        // maximum depth for recursion
        size_t _dmax = std::max<size_t>(_ncore, 4);
        // how many threads we need to initialize
        size_t _num = _tids.size();
        // create a task arena
        auto _arena = get_task_arena();
        // this is the task passed to the task-group
        std::function<void()> _exec_task;
        _exec_task = [&]() {
            add_thread_id();
            static thread_local size_type _depth    = 0;
            int                           _ret      = 0;
            auto                          _this_tid = std::this_thread::get_id();
            // don't let the main thread execute the function
            if(_tids.count(_this_tid) > 0)
            {
                // execute the function
                _ret = _exec();
                // add the result
                _total_exec += _ret;
            }
            // if the function did not return anything, recursively execute
            // two more tasks
            ++_depth;
            if(_ret == 0 && _depth < _dmax && _total_exec.load() < _num)
            {
                tbb::task_group tg{};
                tg.run([&]() { _exec_task(); });
                tg.run([&]() { _exec_task(); });
                ThisThread::sleep_for(std::chrono::milliseconds{ 1 });
                tg.wait();
            }
            --_depth;
        };

        // TBB won't oversubscribe so we need to limit by ncores - 1
        size_t nitr        = 0;
        auto   _fname      = __FUNCTION__;
        auto   _write_info = [&]() {
            std::cout << "[" << _fname << "]> Total executed: " << _total_exec
                      << ", expected: " << _num << ", size: " << size() << std::endl;
        };
        while(_total_exec < _num)
        {
            auto _n = 2 * _num;
            while(--_n > 0)
            {
                _arena->execute(
                    [&]() { m_tbb_task_group->run([&]() { _exec_task(); }); });
            }
            _arena->execute([&]() { m_tbb_task_group->wait(); });
            // don't loop infinitely but use a strict condition
            if(nitr++ > 2 * (_num + 1) && (_total_exec - 1) == _num)
            {
                _write_info();
                break;
            }
            // at this point we need to exit
            if(nitr > 8 * (_num + 1))
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
        get_queue()->ExecuteOnSpecificThreads(_tids, this, std::forward<FuncT>(_func));
    }
}

//======================================================================================//

}  // namespace PTL
