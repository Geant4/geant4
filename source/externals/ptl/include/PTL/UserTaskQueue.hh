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
#include "PTL/TaskSubQueue.hh"
#include "PTL/Threading.hh"
#include "PTL/VTask.hh"
#include "PTL/VUserTaskQueue.hh"

#include <atomic>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

namespace PTL
{
class UserTaskQueue : public VUserTaskQueue
{
public:
    using task_pointer          = std::shared_ptr<VTask>;
    using TaskSubQueueContainer = std::vector<TaskSubQueue*>;
    using random_engine_t       = std::default_random_engine;
    using int_dist_t            = std::uniform_int_distribution<int>;

public:
    // Constructor and Destructors
    UserTaskQueue(intmax_t nworkers = -1, UserTaskQueue* = nullptr);
    // Virtual destructors are required by abstract classes
    // so add it by default, just in case
    ~UserTaskQueue() override;

public:
    // Virtual  function for getting a task from the queue
    task_pointer GetTask(intmax_t subq = -1, intmax_t nitr = -1) override;
    // Virtual function for inserting a task into the queue
    intmax_t InsertTask(task_pointer&&, ThreadData* = nullptr,
                        intmax_t subq = -1) override PTL_NO_SANITIZE_THREAD;

    // if executing only tasks in threads bin
    task_pointer GetThreadBinTask();

    // Overload this function to hold threads
    void Wait() override {}
    void resize(intmax_t) override;

    bool      empty() const override;
    size_type size() const override;

    size_type bin_size(size_type bin) const override;
    bool      bin_empty(size_type bin) const override;

    bool      true_empty() const override;
    size_type true_size() const override;

    void ExecuteOnAllThreads(ThreadPool* tp, function_type f) override;

    void ExecuteOnSpecificThreads(ThreadIdSet tid_set, ThreadPool* tp,
                                  function_type f) override;

    VUserTaskQueue* clone() override;

    intmax_t GetThreadBin() const override;

protected:
    intmax_t GetInsertBin() const;

private:
    void AcquireHold();
    void ReleaseHold();

private:
    bool                       m_is_clone;
    intmax_t                   m_thread_bin;
    mutable intmax_t           m_insert_bin;
    std::atomic_bool*          m_hold      = nullptr;
    std::atomic_uintmax_t*     m_ntasks    = nullptr;
    Mutex*                     m_mutex     = nullptr;
    TaskSubQueueContainer*     m_subqueues = nullptr;
    std::vector<int>           m_rand_list = {};
    std::vector<int>::iterator m_rand_itr  = {};
};

//======================================================================================//

inline bool
UserTaskQueue::empty() const
{
    return (m_ntasks->load(std::memory_order_relaxed) == 0);
}

//======================================================================================//

inline UserTaskQueue::size_type
UserTaskQueue::size() const
{
    return m_ntasks->load(std::memory_order_relaxed);
}

//======================================================================================//

inline UserTaskQueue::size_type
UserTaskQueue::bin_size(size_type bin) const
{
    return (*m_subqueues)[bin]->size();
}

//======================================================================================//

inline bool
UserTaskQueue::bin_empty(size_type bin) const
{
    return (*m_subqueues)[bin]->empty();
}

//======================================================================================//

inline bool
UserTaskQueue::true_empty() const
{
    for(const auto& itr : *m_subqueues)
        if(!itr->empty())
            return false;
    return true;
}

//======================================================================================//

inline UserTaskQueue::size_type
UserTaskQueue::true_size() const
{
    size_type _n = 0;
    for(const auto& itr : *m_subqueues)
        _n += itr->size();
    return _n;
}

//======================================================================================//
}  // namespace PTL
