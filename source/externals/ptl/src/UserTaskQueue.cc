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
//  Tasking class implementation
//  Class Description:
//  ---------------------------------------------------------------
//  Author: Jonathan Madsen
//  ---------------------------------------------------------------

#include "PTL/UserTaskQueue.hh"

#include "PTL/AutoLock.hh"
#include "PTL/TaskGroup.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Utility.hh"

#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <utility>

using namespace PTL;

//======================================================================================//

UserTaskQueue::UserTaskQueue(intmax_t nworkers, UserTaskQueue* parent)
: VUserTaskQueue(nworkers)
, m_is_clone((parent) != nullptr)
, m_thread_bin((parent) ? (ThreadPool::get_this_thread_id() % (nworkers + 1)) : 0)
, m_insert_bin((parent) ? (ThreadPool::get_this_thread_id() % (nworkers + 1)) : 0)
, m_hold((parent) ? parent->m_hold : new std::atomic_bool(false))
, m_ntasks((parent) ? parent->m_ntasks : new std::atomic_uintmax_t(0))
, m_mutex((parent) ? parent->m_mutex : new Mutex{})
, m_subqueues((parent) ? parent->m_subqueues : new TaskSubQueueContainer())
{
    // create nthreads + 1 subqueues so there is always a subqueue available
    if(!parent)
    {
        for(intmax_t i = 0; i < nworkers + 1; ++i)
            m_subqueues->emplace_back(new TaskSubQueue(m_ntasks));
    }

#if defined(DEBUG)
    if(GetEnv<int>("PTL_VERBOSE", 0) > 3)
    {
        RecursiveAutoLock l(TypeMutex<decltype(std::cout), RecursiveMutex>());
        std::stringstream ss;
        ss << ThreadPool::get_this_thread_id() << "> " << ThisThread::get_id() << " ["
           << __FUNCTION__ << ":" << __LINE__ << "] "
           << "this = " << this << ", "
           << "clone = " << std::boolalpha << m_is_clone << ", "
           << "thread = " << m_thread_bin << ", "
           << "insert = " << m_insert_bin << ", "
           << "hold = " << m_hold->load() << " @ " << m_hold << ", "
           << "tasks = " << m_ntasks->load() << " @ " << m_ntasks << ", "
           << "subqueue = " << m_subqueues << ", "
           << "size = " << true_size() << ", "
           << "empty = " << true_empty();
        std::cout << ss.str() << std::endl;
    }
#endif
}

//======================================================================================//

UserTaskQueue::~UserTaskQueue()
{
    if(!m_is_clone)
    {
        for(auto& itr : *m_subqueues)
        {
            assert(itr->empty());
            delete itr;
        }
        m_subqueues->clear();
        delete m_hold;
        delete m_ntasks;
        delete m_mutex;
        delete m_subqueues;
    }
}

//======================================================================================//

void
UserTaskQueue::resize(intmax_t n)
{
    if(!m_mutex)
        throw std::runtime_error("nullptr to mutex");
    AutoLock lk(m_mutex);
    if(m_workers < n)
    {
        while(m_workers < n)
        {
            m_subqueues->emplace_back(new TaskSubQueue(m_ntasks));
            ++m_workers;
        }
    }
    else if(m_workers > n)
    {
        while(m_workers > n)
        {
            delete m_subqueues->back();
            m_subqueues->pop_back();
            --m_workers;
        }
    }
}

//======================================================================================//

VUserTaskQueue*
UserTaskQueue::clone()
{
    return new UserTaskQueue(workers(), this);
}
//======================================================================================//

intmax_t
UserTaskQueue::GetThreadBin() const
{
    // get a thread id number
    static thread_local intmax_t tl_bin =
        (m_thread_bin + ThreadPool::get_this_thread_id()) % (m_workers + 1);
    return tl_bin;
}

//======================================================================================//

intmax_t
UserTaskQueue::GetInsertBin() const
{
    return (++m_insert_bin % (m_workers + 1));
}

//======================================================================================//

UserTaskQueue::task_pointer
UserTaskQueue::GetThreadBinTask()
{
    intmax_t      tbin      = GetThreadBin();
    TaskSubQueue* task_subq = (*m_subqueues)[tbin % (m_workers + 1)];
    task_pointer  _task     = nullptr;

    //------------------------------------------------------------------------//
    auto get_task = [&]() {
        if(task_subq->AcquireClaim())
        {
            // run task
            _task = task_subq->PopTask(true);
            // release the claim on the bin
            task_subq->ReleaseClaim();
        }
        if(_task)
            --(*m_ntasks);
        // return success if valid pointer
        return (_task != nullptr);
    };
    //------------------------------------------------------------------------//

    // while not empty
    while(!task_subq->empty())
    {
        if(get_task())
            break;
    }
    return _task;
}

//======================================================================================//

UserTaskQueue::task_pointer
UserTaskQueue::GetTask(intmax_t subq, intmax_t nitr)
{
    // exit if empty
    if(this->true_empty())
        return nullptr;

    // ensure the thread has a bin assignment
    intmax_t tbin = GetThreadBin();
    intmax_t n    = (subq < 0) ? tbin : subq;
    if(nitr < 1)
        nitr = (m_workers + 1);  // * m_ntasks->load(std::memory_order_relaxed);

    if(m_hold->load(std::memory_order_relaxed))
    {
        return GetThreadBinTask();
    }

    task_pointer _task = nullptr;
    //------------------------------------------------------------------------//
    auto get_task = [&](intmax_t _n) {
        TaskSubQueue* task_subq = (*m_subqueues)[_n % (m_workers + 1)];
        // try to acquire a claim for the bin
        // if acquired, no other threads will access bin until claim is released
        if(!task_subq->empty() && task_subq->AcquireClaim())
        {
            // pop task out of bin
            _task = task_subq->PopTask(n == tbin);
            // release the claim on the bin
            task_subq->ReleaseClaim();
        }
        if(_task)
            --(*m_ntasks);
        // return success if valid pointer
        return (_task != nullptr);
    };
    //------------------------------------------------------------------------//

    // there are num_workers+1 bins so there is always a bin that is open
    // execute num_workers+2 iterations so the thread checks its bin twice
    // while(!empty())
    {
        for(intmax_t i = 0; i < nitr; ++i, ++n)
        {
            if(get_task(n % (m_workers + 1)))
                return _task;
        }
    }

    // only reached if looped over all bins (and looked in own bin twice)
    // and found no work so return an empty task and the thread will be put to
    // sleep if there is still no work by the time it reaches its
    // condition variable
    return _task;
}

//======================================================================================//

intmax_t
UserTaskQueue::InsertTask(task_pointer&& task, ThreadData* data, intmax_t subq)
{
    // increment number of tasks
    ++(*m_ntasks);

    bool     spin = m_hold->load(std::memory_order_relaxed);
    intmax_t tbin = GetThreadBin();

    if(data && data->within_task)
    {
        subq = tbin;
        // spin = true;
    }

    // subq is -1 unless specified so unless specified
    // GetInsertBin() call increments a counter and returns
    // counter % (num_workers + 1) so that tasks are distributed evenly
    // among the bins
    intmax_t n = (subq < 0) ? GetInsertBin() : subq;

    //------------------------------------------------------------------------//
    auto insert_task = [&](intmax_t _n) {
        TaskSubQueue* task_subq = (*m_subqueues)[_n];
        // TaskSubQueue* next_subq = (*m_subqueues)[(_n + 1) % (m_workers + 1)];
        // if not threads bin and size difference, insert into smaller
        // if(n != tbin && next_subq->size() < task_subq->size())
        //    task_subq = next_subq;
        // try to acquire a claim for the bin
        // if acquired, no other threads will access bin until claim is released
        if(task_subq->AcquireClaim())
        {
            // push the task into the bin
            task_subq->PushTask(std::move(task));
            // release the claim on the bin
            task_subq->ReleaseClaim();
            // return success
            return true;
        }
        return false;
    };
    //------------------------------------------------------------------------//

    // if not in "hold/spin mode", where thread only inserts tasks into
    // specified bin, then move onto next bin
    //
    if(spin)
    {
        n = n % (m_workers + 1);
        while(!insert_task(n))
            ;
        return n;
    }

    // there are num_workers+1 bins so there is always a bin that is open
    // execute num_workers+2 iterations so the thread checks its bin twice
    while(true)
    {
        auto _n = (n++) % (m_workers + 1);
        if(insert_task(_n))
            return _n;
    }
    return GetThreadBin();
}

//======================================================================================//

void
UserTaskQueue::ExecuteOnAllThreads(ThreadPool* tp, function_type func)
{
    using task_group_type      = TaskGroup<int, int>;
    using thread_execute_map_t = std::map<int64_t, bool>;

    if(!tp->is_alive())
    {
        func();
        return;
    }

    task_group_type tg{ [](int& ref, int i) { return (ref += i); }, tp };

    // wait for all threads to finish any work
    // NOTE: will cause deadlock if called from a task
    while(tp->get_active_threads_count() > 0)
        ThisThread::sleep_for(std::chrono::milliseconds(10));

    thread_execute_map_t                thread_execute_map{};
    std::vector<std::shared_ptr<VTask>> _tasks{};
    _tasks.reserve(m_workers + 1);

    AcquireHold();
    for(int i = 0; i < (m_workers + 1); ++i)
    {
        if(i == GetThreadBin())
            continue;

        //--------------------------------------------------------------------//
        auto thread_specific_func = [&]() {
            ScopeDestructor _dtor = tg.get_scope_destructor();
            static Mutex    _mtx;
            _mtx.lock();
            bool& _executed = thread_execute_map[GetThreadBin()];
            _mtx.unlock();
            if(!_executed)
            {
                func();
                _executed = true;
                return 1;
            }
            return 0;
        };
        //--------------------------------------------------------------------//

        InsertTask(tg.wrap(thread_specific_func), ThreadData::GetInstance(), i);
    }

    tp->notify_all();
    int nexecuted = tg.join();
    if(nexecuted != m_workers)
    {
        std::stringstream msg;
        msg << "Failure executing routine on all threads! Only " << nexecuted
            << " threads executed function out of " << m_workers << " workers";
        std::cerr << msg.str() << std::endl;
    }
    ReleaseHold();
}

//======================================================================================//

void
UserTaskQueue::ExecuteOnSpecificThreads(ThreadIdSet tid_set, ThreadPool* tp,
                                        function_type func)
{
    using task_group_type      = TaskGroup<int, int>;
    using thread_execute_map_t = std::map<int64_t, bool>;

    task_group_type tg{ [](int& ref, int i) { return (ref += i); }, tp };

    // wait for all threads to finish any work
    // NOTE: will cause deadlock if called from a task
    while(tp->get_active_threads_count() > 0)
        ThisThread::sleep_for(std::chrono::milliseconds(10));

    if(!tp->is_alive())
    {
        func();
        return;
    }

    thread_execute_map_t thread_execute_map{};

    //========================================================================//
    // wrap the function so that it will only be executed if the thread
    // has an ID in the set
    auto thread_specific_func = [&]() {
        ScopeDestructor _dtor = tg.get_scope_destructor();
        static Mutex    _mtx;
        _mtx.lock();
        bool& _executed = thread_execute_map[GetThreadBin()];
        _mtx.unlock();
        if(!_executed && tid_set.count(ThisThread::get_id()) > 0)
        {
            func();
            _executed = true;
            return 1;
        }
        return 0;
    };
    //========================================================================//

    if(tid_set.count(ThisThread::get_id()) > 0)
        func();

    AcquireHold();
    for(int i = 0; i < (m_workers + 1); ++i)
    {
        if(i == GetThreadBin())
            continue;

        InsertTask(tg.wrap(thread_specific_func), ThreadData::GetInstance(), i);
    }
    tp->notify_all();
    decltype(tid_set.size()) nexecuted = tg.join();
    if(nexecuted != tid_set.size())
    {
        std::stringstream msg;
        msg << "Failure executing routine on specific threads! Only " << nexecuted
            << " threads executed function out of " << tid_set.size() << " workers";
        std::cerr << msg.str() << std::endl;
    }
    ReleaseHold();
}

//======================================================================================//

void
UserTaskQueue::AcquireHold()
{
    bool _hold;
    while(!(_hold = m_hold->load(std::memory_order_relaxed)))
    {
        m_hold->compare_exchange_strong(_hold, true, std::memory_order_release,
                                        std::memory_order_relaxed);
    }
}

//======================================================================================//

void
UserTaskQueue::ReleaseHold()
{
    bool _hold;
    while((_hold = m_hold->load(std::memory_order_relaxed)))
    {
        m_hold->compare_exchange_strong(_hold, false, std::memory_order_release,
                                        std::memory_order_relaxed);
    }
}

//======================================================================================//
