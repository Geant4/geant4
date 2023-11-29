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
//      Abstract base class for creating a task queue used by
//      ThreadPool
//  ---------------------------------------------------------------
//  Author: Jonathan Madsen
//  ---------------------------------------------------------------

#pragma once

#include "PTL/Globals.hh"
#include "PTL/Threading.hh"

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <set>

namespace PTL
{
class VTask;
class ThreadPool;
class ThreadData;

class VUserTaskQueue
{
public:
    using task_pointer  = std::shared_ptr<VTask>;
    using AtomicInt     = std::atomic<intmax_t>;
    using size_type     = uintmax_t;
    using function_type = std::function<void()>;
    using ThreadIdSet   = std::set<ThreadId>;

public:
    // Constructor - accepting the number of workers
    explicit VUserTaskQueue(intmax_t nworkers = -1);
    // Virtual destructors are required by abstract classes
    // so add it by default, just in case
    virtual ~VUserTaskQueue() = default;

public:
    // Virtual function for getting a task from the queue
    // parameters:
    //      1. int - get from specific sub-queue
    //      2. int - number of iterations
    // returns:
    //      VTask* - a task or nullptr
    virtual task_pointer GetTask(intmax_t subq = -1, intmax_t nitr = -1) = 0;

    // Virtual function for inserting a task into the queue
    // parameters:
    //      1. VTask* - task to insert
    //      2. int - sub-queue to inserting into
    // return:
    //      int - subqueue inserted into
    virtual intmax_t InsertTask(task_pointer&&, ThreadData* = nullptr,
                                intmax_t subq = -1) PTL_NO_SANITIZE_THREAD = 0;

    // Overload this function to hold threads
    virtual void     Wait()               = 0;
    virtual intmax_t GetThreadBin() const = 0;

    virtual void resize(intmax_t) = 0;

    // these are used for stanard checking
    virtual size_type size() const  = 0;
    virtual bool      empty() const = 0;

    virtual size_type bin_size(size_type bin) const  = 0;
    virtual bool      bin_empty(size_type bin) const = 0;

    // these are for slower checking, default to returning normal size()/empty
    virtual size_type true_size() const { return size(); }
    virtual bool      true_empty() const { return empty(); }

    // a method of executing a specific function on all threads
    virtual void ExecuteOnAllThreads(ThreadPool* tp, function_type f) = 0;

    virtual void ExecuteOnSpecificThreads(ThreadIdSet tid_set, ThreadPool* tp,
                                          function_type f) = 0;

    intmax_t workers() const { return m_workers; }

    virtual VUserTaskQueue* clone() = 0;

protected:
    intmax_t m_workers = 0;
};

}  // namespace PTL
