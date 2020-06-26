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
#include "PTL/AutoLock.hh"
#include "PTL/Task.hh"
#include "PTL/TaskGroup.hh"
#include "PTL/TaskManager.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Threading.hh"
#include "PTL/Utility.hh"

#include <cstdlib>
#include <cstring>
#include <iterator>

using namespace PTL;

//======================================================================================//

TaskRunManager::pointer&
TaskRunManager::GetPrivateMasterRunManager(bool init, bool useTBB)
{
    static pointer _instance = (init) ? new TaskRunManager(useTBB) : nullptr;
    return _instance;
}

//======================================================================================//

TaskRunManager*
TaskRunManager::GetMasterRunManager(bool useTBB)
{
    static pointer& _instance = GetPrivateMasterRunManager(true, useTBB);
    return _instance;
}

//======================================================================================//

TaskRunManager*
TaskRunManager::GetInstance(bool useTBB)
{
    return GetMasterRunManager(useTBB);
}

//======================================================================================//

TaskRunManager::TaskRunManager(bool useTBB)
: m_is_initialized(false)
, m_verbose(0)
, m_workers(std::thread::hardware_concurrency())
, m_task_queue(nullptr)
, m_thread_pool(nullptr)
, m_task_manager(nullptr)
{
    if(!GetPrivateMasterRunManager(false))
    {
        GetPrivateMasterRunManager(false) = this;
    }

#ifdef PTL_USE_TBB
    auto _useTBB = GetEnv<bool>("FORCE_TBB", useTBB);
    if(_useTBB)
        useTBB = true;
#endif

    // handle TBB
    ThreadPool::set_use_tbb(useTBB);
    m_workers = GetEnv<uint64_t>("PTL_NUM_THREADS", m_workers);
}

//======================================================================================//

TaskRunManager::~TaskRunManager() {}

//======================================================================================//

void
TaskRunManager::Initialize(uint64_t n)
{
    m_workers = n;

    // create threadpool if needed + task manager
    if(!m_thread_pool)
    {
        if(m_verbose > 0)
            std::cout << "TaskRunManager :: Creating thread pool..." << std::endl;
        m_thread_pool = new ThreadPool(m_workers, m_task_queue);
        if(m_verbose > 0)
            std::cout << "TaskRunManager :: Creating task manager..." << std::endl;
        m_task_manager = new TaskManager(m_thread_pool);
    }
    // or resize
    else if(m_workers != m_thread_pool->size())
    {
        if(m_verbose > 0)
        {
            std::cout << "TaskRunManager :: Resizing thread pool from "
                      << m_thread_pool->size() << " to " << m_workers << " threads ..."
                      << std::endl;
        }
        m_thread_pool->resize(m_workers);
    }

    // create the joiners
    if(ThreadPool::using_tbb())
    {
        if(m_verbose > 0)
            std::cout << "TaskRunManager :: Using TBB..." << std::endl;
    }
    else
    {
        if(m_verbose > 0)
            std::cout << "TaskRunManager :: Using ThreadPool..." << std::endl;
    }

    m_is_initialized = true;
    if(m_verbose > 0)
        std::cout << "TaskRunManager :: initialized..." << std::endl;
}

//======================================================================================//

void
TaskRunManager::Terminate()
{
    m_is_initialized = false;
    m_thread_pool->destroy_threadpool();
    delete m_task_manager;
    delete m_thread_pool;
    m_task_manager = nullptr;
    m_thread_pool  = nullptr;
}

//======================================================================================//
