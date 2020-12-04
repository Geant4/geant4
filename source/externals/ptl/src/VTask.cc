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
// This file creates an abstract base class for the thread-pool tasking
// system
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#include "PTL/VTask.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/VTaskGroup.hh"

using namespace PTL;

//======================================================================================//

VTask::VTask()
: m_depth(0)
, m_group(nullptr)
, m_pool(nullptr)
{}

//======================================================================================//

VTask::VTask(VTaskGroup* task_group)
: m_depth(0)
, m_group(task_group)
, m_pool((m_group) ? task_group->pool() : nullptr)
{}

//======================================================================================//

VTask::VTask(ThreadPool* tp)
: m_depth(0)
, m_group(nullptr)
, m_pool(tp)
{}

//======================================================================================//

VTask::~VTask() {}

//======================================================================================//

void
VTask::operator--()
{
    if(m_group)
    {
        intmax_t _count = --(*m_group);
        if(_count < 2)
        {
            try
            {
                m_group->task_cond()->notify_all();
            } catch(std::system_error& e)
            {
                auto     tid = ThreadPool::get_this_thread_id();
                AutoLock l(TypeMutex<decltype(std::cerr)>(), std::defer_lock);
                if(!l.owns_lock())
                    l.lock();
                std::cerr << "[" << tid << "] Caught system error: " << e.what()
                          << std::endl;
            }
        }
    }
}

//======================================================================================//

bool
VTask::is_native_task() const
{
    return (m_group) ? m_group->is_native_task_group() : false;
}

//======================================================================================//

ThreadPool*
VTask::pool() const
{
    return (!m_pool && m_group) ? m_group->pool() : m_pool;
}

//======================================================================================//
