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
// class description:
//   This is a class for run control in Tasking for multi-threaded runs
//   It extends RunManager re-implementing multi-threaded behavior in
//   key methods. See documentation for RunManager
//   Users initializes an instance of this class instead of RunManager
//   to start a multi-threaded simulation.

#pragma once

#include "PTL/ThreadPool.hh"
#include "PTL/VUserTaskQueue.hh"

#include <cstddef>
#include <cstdint>
#include <thread>

namespace PTL
{
class TaskManager;

//======================================================================================//

class TaskRunManager
{
public:
    using pointer = TaskRunManager*;

public:
    // Parameters:
    //      m_task_queue: provide a custom task queue
    //      useTBB: only relevant if PTL_USE_TBB defined
    //      grainsize:  0 = auto
    explicit TaskRunManager(bool useTBB = false);
    virtual ~TaskRunManager();

public:
    virtual int GetNumberOfThreads() const
    {
        return (m_thread_pool) ? (int)m_thread_pool->size() : 0;
    }
    virtual size_t GetNumberActiveThreads() const
    {
        return (m_thread_pool) ? m_thread_pool->size() : 0;
    }

public:
    // Inherited methods to re-implement for MT case
    virtual void Initialize(uint64_t n = std::thread::hardware_concurrency());
    virtual void Terminate();
    ThreadPool*  GetThreadPool() const { return m_thread_pool; }
    TaskManager* GetTaskManager() const { return m_task_manager; }
    bool         IsInitialized() const { return m_is_initialized; }
    int          GetVerbose() const { return m_verbose; }
    void         SetVerbose(int val) { m_verbose = val; }

public:  // with description
    // Singleton implementing master thread behavior
    static TaskRunManager* GetInstance(bool useTBB = false);
    static TaskRunManager* GetMasterRunManager(bool useTBB = false);

private:
    static pointer& GetPrivateMasterRunManager();
    static pointer& GetPrivateMasterRunManager(bool init, bool useTBB = false);

protected:
    // Barriers: synch points between master and workers
    bool            m_is_initialized = false;
    int             m_verbose        = 0;
    uint64_t        m_workers        = 0;
    VUserTaskQueue* m_task_queue     = nullptr;
    ThreadPool*     m_thread_pool    = nullptr;
    TaskManager*    m_task_manager   = nullptr;
};

}  // namespace PTL
