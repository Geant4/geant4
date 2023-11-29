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

#ifndef G4GMAKE
#include "PTL/Config.hh"  // IWYU pragma: keep
#endif
#include "PTL/Globals.hh"
#include "PTL/VTask.hh"
#if defined(PTL_USE_LOCKS)
#    include "PTL/AutoLock.hh"
#endif

#include <atomic>
#include <cassert>
#include <list>
#include <memory>
#include <utility>

namespace PTL
{
class TaskSubQueue
{
public:
    template <typename Tp>
    using container = std::list<Tp>;

    using task_pointer   = std::shared_ptr<VTask>;
    using container_type = container<task_pointer>;
    using size_type      = container_type::size_type;

public:
    TaskSubQueue(std::atomic_uintmax_t* _ntasks);
    TaskSubQueue(const TaskSubQueue&);
    ~TaskSubQueue() = default;

    TaskSubQueue& operator=(const TaskSubQueue&) = delete;

public:
    int GetId() const;

    bool AcquireClaim();
    void ReleaseClaim();

    void         PushTask(task_pointer&&) PTL_NO_SANITIZE_THREAD;
    task_pointer PopTask(bool front = true) PTL_NO_SANITIZE_THREAD;

    size_type size() const;
    bool      empty() const;

private:
    // mutex
#if defined(PTL_USE_LOCKS)
    Mutex m_mutex{};
#endif
    // used internally to keep number of tasks
    std::atomic<size_type> m_ntasks;
    // for checking if being modified
    std::atomic_bool m_available;
    // used my master queue to keep track of number of tasks
    std::atomic_uintmax_t* m_all_tasks;
    // queue of tasks
    container_type m_task_queue;
};

//======================================================================================//

inline TaskSubQueue::TaskSubQueue(std::atomic_uintmax_t* _ntasks)
: m_ntasks(0)
, m_available(true)
, m_all_tasks(_ntasks)
{}

//======================================================================================//

inline TaskSubQueue::TaskSubQueue(const TaskSubQueue& rhs)
: m_ntasks(0)
, m_available(true)
, m_all_tasks(rhs.m_all_tasks)
{}

//======================================================================================//

inline bool
TaskSubQueue::AcquireClaim()
{
    bool is_avail = m_available.load(std::memory_order_relaxed);
    if(!is_avail)
        return false;
    return m_available.compare_exchange_strong(is_avail, false,
                                               std::memory_order_relaxed);
}

//======================================================================================//

inline void
TaskSubQueue::ReleaseClaim()
{
    // if(m_available.load(std::memory_order_relaxed))
    //    return;
    m_available.store(true, std::memory_order_release);
}

//======================================================================================//

inline TaskSubQueue::size_type
TaskSubQueue::size() const
{
    return m_ntasks.load();
}

//======================================================================================//

inline bool
TaskSubQueue::empty() const
{
    return (m_ntasks.load() == 0);
}

//======================================================================================//

inline void
TaskSubQueue::PushTask(task_pointer&& task)
{
    // no need to lock these if claim is acquired via atomic
    assert(m_available.load(std::memory_order_relaxed) == false);
    ++m_ntasks;
#if defined(PTL_USE_LOCKS)
    AutoLock lk{ m_mutex };
#endif
    m_task_queue.emplace_front(std::move(task));
}

//======================================================================================//

inline TaskSubQueue::task_pointer
TaskSubQueue::PopTask(bool front)
{
    // no need to lock -- claim is acquired via atomic
    assert(m_available.load(std::memory_order_relaxed) == false);
    if(m_ntasks.load() == 0)
        return nullptr;

    task_pointer _task{ nullptr };
    if(front)
    {
#if defined(PTL_USE_LOCKS)
        AutoLock lk{ m_mutex };
#endif
        _task = std::move(m_task_queue.front());
        m_task_queue.pop_front();
    }
    else
    {
#if defined(PTL_USE_LOCKS)
        AutoLock lk{ m_mutex };
#endif
        _task = std::move(m_task_queue.back());
        m_task_queue.pop_back();
    }
    --m_ntasks;

    return _task;
}

//======================================================================================//
}  // namespace PTL
