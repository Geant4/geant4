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
#include "PTL/TaskRunManager.hh"

#include "PTL/TaskManager.hh"
#include "PTL/ThreadPool.hh"

namespace PTL
{
//======================================================================================//

TaskRunManager::pointer&
TaskRunManager::GetPrivateMasterRunManager()
{
    static pointer _instance = nullptr;
    return _instance;
}

//======================================================================================//

TaskRunManager::pointer&
TaskRunManager::GetPrivateMasterRunManager(bool init, bool useTBB)
{
    auto& _v = GetPrivateMasterRunManager();
    if(!init)
        return _v;
    if(!_v)
        _v = new TaskRunManager(useTBB);
    return _v;
}

//======================================================================================//

TaskRunManager*
TaskRunManager::GetMasterRunManager(bool useTBB)
{
    auto& _v = GetPrivateMasterRunManager(true, useTBB);
    return _v;
}

//======================================================================================//

TaskRunManager*
TaskRunManager::GetInstance(bool useTBB)
{
    return GetMasterRunManager(useTBB);
}

//======================================================================================//

TaskRunManager::TaskRunManager(bool useTBB)
: m_workers(std::thread::hardware_concurrency())
, m_use_tbb(useTBB)
{
    if(!GetPrivateMasterRunManager())
        GetPrivateMasterRunManager() = this;
}

//======================================================================================//

TaskRunManager::~TaskRunManager()
{
    if(GetPrivateMasterRunManager() == this)
        GetPrivateMasterRunManager() = nullptr;
}

//======================================================================================//

void
TaskRunManager::Initialize(uint64_t n)
{
    m_workers = n;

    // create threadpool if needed + task manager
    if(!m_thread_pool)
    {
        ThreadPool::Config cfg;
        cfg.pool_size  = m_workers;
        cfg.task_queue = m_task_queue;
        cfg.use_tbb    = m_use_tbb;
        m_thread_pool  = new ThreadPool(cfg);
        m_task_manager = new TaskManager(m_thread_pool);
    }
    // or resize
    else if(m_workers != m_thread_pool->size())
    {
        m_thread_pool->resize(m_workers);
    }

    m_is_initialized = true;
}

//======================================================================================//

void
TaskRunManager::Terminate()
{
    m_is_initialized = false;
    if(m_thread_pool)
        m_thread_pool->destroy_threadpool();
    delete m_task_manager;
    delete m_thread_pool;
    m_task_manager = nullptr;
    m_thread_pool  = nullptr;
}

//======================================================================================//

}  // namespace PTL
