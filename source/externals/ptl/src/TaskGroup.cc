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

#include "PTL/TaskGroup.hh"

#include "PTL/TaskRunManager.hh"
#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"

//======================================================================================//

namespace PTL
{
namespace internal
{
std::atomic_uintmax_t&
task_group_counter()
{
    static std::atomic_uintmax_t _instance(0);
    return _instance;
}

ThreadPool*
get_default_threadpool()
{
    auto* mrm = TaskRunManager::GetMasterRunManager();
    if(mrm)
    {
        if(!mrm->GetThreadPool())
            mrm->Initialize();
        return mrm->GetThreadPool();
    }
    return nullptr;
}

intmax_t
get_task_depth()
{
    return (ThreadData::GetInstance()) ? ThreadData::GetInstance()->task_depth : 0;
}
}  // namespace internal
}  // namespace PTL

//======================================================================================//
