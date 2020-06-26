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

#include "PTL/ThreadData.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/Threading.hh"
#include "PTL/VUserTaskQueue.hh"

#include <iostream>

using namespace PTL;

//======================================================================================//

ThreadData*&
ThreadData::GetInstance()
{
    static thread_local ThreadData* _instance = nullptr;
    return _instance;
}

//======================================================================================//

ThreadData::ThreadData(ThreadPool* tp)
: is_master(tp->is_master())
, within_task(false)
, task_depth(0)
, thread_pool(tp)
, current_queue(tp->get_queue())
, queue_stack({ current_queue })
{}

//======================================================================================//

void
ThreadData::update()
{
    current_queue = thread_pool->get_queue();
    queue_stack.push_back(current_queue);
}

//======================================================================================//

ThreadData::~ThreadData() {}

//======================================================================================//
