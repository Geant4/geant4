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

#include <cstddef>
#include <cstdint>
#include <deque>

#if defined(PTL_USE_TBB)
#    include <tbb/global_control.h>
#    include <tbb/task_group.h>
#    include <tbb/task_scheduler_init.h>
#endif

namespace PTL
{
//--------------------------------------------------------------------------------------//

#if defined(PTL_USE_TBB)

using tbb_global_control_t = ::tbb::global_control;
using tbb_task_group_t     = ::tbb::task_group;
#else

namespace tbb
{
class task_group
{
public:
    // dummy constructor
    task_group() {}
    // dummy wait
    inline void wait() {}
    // run function
    template <typename FuncT>
    inline void run(FuncT f)
    {
        f();
    }
    // run and wait
    template <typename FuncT>
    inline void run_and_wait(FuncT f)
    {
        f();
    }
};

class global_control
{
public:
    enum parameter
    {
        max_allowed_parallelism,
        thread_stack_size
    };

    global_control(parameter p, size_t value);
    ~global_control();
    static size_t active_value(parameter param);
};

}  // namespace tbb

using tbb_global_control_t = tbb::global_control;
using tbb_task_group_t     = tbb::task_group;

#endif

//--------------------------------------------------------------------------------------//

class ThreadPool;
class VUserTaskQueue;

//--------------------------------------------------------------------------------------//

class ThreadData
{
public:
    template <typename Tp>
    using TaskStack = std::deque<Tp>;

    ThreadData(ThreadPool* tp);
    ~ThreadData();

    void update();

public:
    bool                       is_master;
    bool                       within_task;
    intmax_t                   task_depth;
    ThreadPool*                thread_pool;
    VUserTaskQueue*            current_queue;
    TaskStack<VUserTaskQueue*> queue_stack;

public:
    // Public functions
    static ThreadData*& GetInstance();
};

//--------------------------------------------------------------------------------------//

}  // namespace PTL
