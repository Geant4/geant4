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
//
// Class Description:
//
// This file creates an abstract base class for the grouping the thread-pool
// tasking system into independently joinable units
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#include "PTL/VTaskGroup.hh"
#include "PTL/Globals.hh"
#include "PTL/Task.hh"
#include "PTL/TaskRunManager.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/VTask.hh"

using namespace PTL;

//======================================================================================//

std::atomic_uintmax_t&
vtask_group_counter()
{
    static std::atomic_uintmax_t _instance(0);
    return _instance;
}

//======================================================================================//

int VTaskGroup::f_verbose = GetEnv<int>("PTL_VERBOSE", 0);

//======================================================================================//

VTaskGroup::VTaskGroup(ThreadPool* tp)
: m_id(vtask_group_counter()++)
, m_pool(tp)
, m_tot_task_count(std::make_shared<atomic_int>(0))
, m_task_cond(std::make_shared<condition_t>())
, m_task_lock(std::make_shared<lock_t>())
, m_main_tid(std::this_thread::get_id())
{
    if(!m_pool && TaskRunManager::GetMasterRunManager())
        m_pool = TaskRunManager::GetMasterRunManager()->GetThreadPool();

    if(!m_pool)
    {
        std::cerr << __FUNCTION__ << "@" << __LINE__ << " :: Warning! "
                  << "nullptr to thread pool!" << std::endl;
    }
}

//======================================================================================//

VTaskGroup::~VTaskGroup() {}

//======================================================================================//

void
VTaskGroup::wait()
{
    // if no pool was initially present at creation
    if(!m_pool)
    {
        // check for master MT run-manager
        if(TaskRunManager::GetMasterRunManager())
            m_pool = TaskRunManager::GetMasterRunManager()->GetThreadPool();

        // if MTRunManager does not exist or no thread pool created
        if(!m_pool)
        {
            if(f_verbose > 0)
            {
                fprintf(stderr, "%s @ %i :: Warning! nullptr to thread-pool (%p)\n",
                        __FUNCTION__, __LINE__, static_cast<void*>(m_pool));
                std::cerr << __FUNCTION__ << "@" << __LINE__ << " :: Warning! "
                          << "nullptr to thread pool!" << std::endl;
            }
            return;
        }
    }

    ThreadData* data = ThreadData::GetInstance();
    if(!data)
        return;

    ThreadPool*     tpool = (m_pool) ? m_pool : data->thread_pool;
    VUserTaskQueue* taskq = (tpool) ? tpool->get_queue() : data->current_queue;

    bool _is_master   = data->is_master;
    bool _within_task = data->within_task;

    auto is_active_state = [&]() {
        return (tpool->state()->load(std::memory_order_relaxed) !=
                thread_pool::state::STOPPED);
    };

    auto execute_this_threads_tasks = [&]() {
        if(!taskq)
            return;

        // only want to process if within a task
        if((!_is_master || tpool->size() < 2) && _within_task)
        {
            int bin = static_cast<int>(taskq->GetThreadBin());
            // const auto nitr = (tpool) ? tpool->size() : Thread::hardware_concurrency();
            while(this->pending() > 0)
            {
                task_pointer _task = taskq->GetTask(bin);
                if(_task)
                    (*_task)();
            }
        }
    };

    // checks for validity
    if(!is_native_task_group())
    {
        // for external threads
        if(!_is_master || tpool->size() < 2)
            return;
    }
    else if(f_verbose > 0)
    {
        if(!tpool || !taskq)
        {
            // something is wrong, didn't create thread-pool?
            fprintf(
                stderr,
                "%s @ %i :: Warning! nullptr to thread data (%p) or task-queue (%p)\n",
                __FUNCTION__, __LINE__, static_cast<void*>(tpool),
                static_cast<void*>(taskq));
        }
        // return if thread pool isn't built
        else if(is_native_task_group() && !tpool->is_alive())
        {
            fprintf(stderr, "%s @ %i :: Warning! thread-pool is not alive!\n",
                    __FUNCTION__, __LINE__);
        }
        else if(!is_active_state())
        {
            fprintf(stderr, "%s @ %i :: Warning! thread-pool is not active!\n",
                    __FUNCTION__, __LINE__);
        }
    }

    intmax_t wake_size = 2;
    AutoLock _lock(*m_task_lock, std::defer_lock);

    while(is_active_state())
    {
        execute_this_threads_tasks();

        // while loop protects against spurious wake-ups
        while(_is_master && pending() > 0 && is_active_state())
        {
            // auto _wake = [&]() { return (wake_size > pending() || !is_active_state());
            // };

            // lock before sleeping on condition
            if(!_lock.owns_lock())
                _lock.lock();

            // Wait until signaled that a task has been competed
            // Unlock mutex while wait, then lock it back when signaled
            // when true, this wakes the thread
            if(pending() >= wake_size)
            {
                m_task_cond->wait(_lock);
            }
            else
            {
                m_task_cond->wait_for(_lock, std::chrono::microseconds(100));
            }
            // unlock
            if(_lock.owns_lock())
                _lock.unlock();
        }

        // if pending is not greater than zero, we are joined
        if(pending() <= 0)
            break;
    }

    if(_lock.owns_lock())
        _lock.unlock();

    intmax_t ntask = this->task_count().load();
    if(ntask > 0)
    {
        std::stringstream ss;
        ss << "\nWarning! Join operation issue! " << ntask << " tasks still "
           << "are running!" << std::endl;
        std::cerr << ss.str();
        this->wait();
    }
}

//======================================================================================//
