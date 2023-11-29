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
// Tasking class header file
//
// Class Description:
//
// This file defines types and macros used to expose Tasking threading model.

#pragma once

#include <array>
#include <cstddef>
#include <future>
#include <mutex>
#include <thread>

namespace PTL
{
// global thread types
using Thread       = std::thread;
using NativeThread = std::thread::native_handle_type;
// std::thread::id does not cast to integer
using Pid_t = std::thread::id;

// Condition
using Condition = std::condition_variable;

// Thread identifier
using ThreadId = Thread::id;

// will be used in the future when migrating threading to task-based style
template <typename Tp>
using Future = std::future<Tp>;
template <typename Tp>
using SharedFuture = std::shared_future<Tp>;
template <typename Tp>
using Promise = std::promise<Tp>;

// global mutex types
using Mutex          = std::mutex;
using RecursiveMutex = std::recursive_mutex;

// static functions: get_id(), sleep_for(...), sleep_until(...), yield(),
namespace ThisThread
{
using namespace std::this_thread;
}

// Helper function for getting a unique static mutex for a specific
// class or type
// Usage example:
//		a template class "Cache<T>" that required a static
//		mutex for specific to type T:
//			AutoLock l(TypeMutex<Cache<T>>());
template <typename Tp, typename MutexTp = Mutex, size_t N = 4>
MutexTp&
TypeMutex(const unsigned int& _n = 0)
{
    static std::array<MutexTp, N> _mutex_array{};
    return _mutex_array[_n % N];
}

//======================================================================================//

namespace Threading
{
enum
{
    SEQUENTIAL_ID    = -2,
    MASTER_ID        = -1,
    WORKER_ID        = 0,
    GENERICTHREAD_ID = -1000
};

Pid_t
GetPidId();

unsigned
GetNumberOfPhysicalCpus();

unsigned
GetNumberOfCores();

int
GetThreadId();

void
SetThreadId(int aNewValue);

bool
SetPinAffinity(int idx);

bool
SetThreadPriority(int _v);

bool
SetPinAffinity(int idx, NativeThread& _t);

bool
SetThreadPriority(int _v, NativeThread& _t);

}  // namespace Threading
}  // namespace PTL
