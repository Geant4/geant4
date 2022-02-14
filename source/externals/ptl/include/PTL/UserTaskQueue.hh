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
//  ---------------------------------------------------------------
//  Tasking class header
//  Class Description:
//  ---------------------------------------------------------------
//  Author: Jonathan Madsen
//  ---------------------------------------------------------------

#pragma once

#include "PTL/Globals.hh"
#include "PTL/Threading.hh"
#include "PTL/Types.hh"
#include "PTL/VUserTaskQueue.hh"

#include <atomic>
#include <deque>
#include <list>
#include <memory>
#include <queue>
#include <random>
#include <set>
#include <stack>

namespace PTL
{
class VTask;
class TaskSubQueue;  // definition in UserTaskQueue.icc

class UserTaskQueue : public VUserTaskQueue
{
public:
    typedef std::shared_ptr<VTask>             task_pointer;
    typedef std::vector<TaskSubQueue*>         TaskSubQueueContainer;
    typedef std::default_random_engine         random_engine_t;
    typedef std::uniform_int_distribution<int> int_dist_t;

public:
    // Constructor and Destructors
    UserTaskQueue(intmax_t nworkers = -1, UserTaskQueue* = nullptr);
    // Virtual destructors are required by abstract classes
    // so add it by default, just in case
    virtual ~UserTaskQueue() override;

public:
    // Virtual  function for getting a task from the queue
    virtual task_pointer GetTask(intmax_t subq = -1, intmax_t nitr = -1) override;
    // Virtual function for inserting a task into the queue
    virtual intmax_t InsertTask(task_pointer&&, ThreadData* = nullptr,
                                intmax_t subq = -1) override PTL_NO_SANITIZE_THREAD;

    // if executing only tasks in threads bin
    task_pointer GetThreadBinTask();

    // Overload this function to hold threads
    virtual void Wait() override {}
    virtual void resize(intmax_t) override;

    virtual bool      empty() const override;
    virtual size_type size() const override;

    virtual size_type bin_size(size_type bin) const override;
    virtual bool      bin_empty(size_type bin) const override;

    inline bool      true_empty() const override;
    inline size_type true_size() const override;

    virtual void ExecuteOnAllThreads(ThreadPool* tp, function_type f) override;

    virtual void ExecuteOnSpecificThreads(ThreadIdSet tid_set, ThreadPool* tp,
                                          function_type f) override;

    virtual VUserTaskQueue* clone() override;

    virtual intmax_t GetThreadBin() const override;

protected:
    intmax_t GetInsertBin() const;

private:
    void AcquireHold();
    void ReleaseHold();

private:
    bool                       m_is_clone;
    intmax_t                   m_thread_bin;
    mutable intmax_t           m_insert_bin;
    std::atomic_bool*          m_hold;
    std::atomic_uintmax_t*     m_ntasks;
    Mutex*                     m_mutex;
    TaskSubQueueContainer*     m_subqueues;
    std::vector<int>           m_rand_list;
    std::vector<int>::iterator m_rand_itr;
};

}  // namespace PTL

//======================================================================================//

#include "PTL/UserTaskQueue.icc"

//======================================================================================//

inline bool
PTL::UserTaskQueue::empty() const
{
    return (m_ntasks->load(std::memory_order_relaxed) == 0);
}

//======================================================================================//

inline PTL::UserTaskQueue::size_type
PTL::UserTaskQueue::size() const
{
    return m_ntasks->load(std::memory_order_relaxed);
}

//======================================================================================//

inline PTL::UserTaskQueue::size_type
PTL::UserTaskQueue::bin_size(size_type bin) const
{
    return (*m_subqueues)[bin]->size();
}

//======================================================================================//

inline bool
PTL::UserTaskQueue::bin_empty(size_type bin) const
{
    return (*m_subqueues)[bin]->empty();
}

//======================================================================================//

inline bool
PTL::UserTaskQueue::true_empty() const
{
    for(const auto& itr : *m_subqueues)
        if(!itr->empty())
            return false;
    return true;
}

//======================================================================================//

inline PTL::UserTaskQueue::size_type
PTL::UserTaskQueue::true_size() const
{
    size_type _n = 0;
    for(const auto& itr : *m_subqueues)
        _n += itr->size();
    return _n;
}

//======================================================================================//
