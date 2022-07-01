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

#include "PTL/Globals.hh"
#include "PTL/Types.hh"

#include <array>
#include <chrono>
#include <condition_variable>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

namespace PTL
{
// Macro to put current thread to sleep
//
#define THREADSLEEP(tick) std::this_thread::sleep_for(std::chrono::seconds(tick))

// will be used in the future when migrating threading to task-based style
template <typename Tp>
using Future = std::future<Tp>;
template <typename Tp>
using SharedFuture = std::shared_future<Tp>;
template <typename Tp>
using Promise = std::promise<Tp>;

//
//          NOTE ON Tasking SERIAL BUILDS AND MUTEX/UNIQUE_LOCK
//          ==================================================
//
// Mutex and RecursiveMutex are always C++11 std::mutex types
// however, in serial mode, using MUTEXLOCK and MUTEXUNLOCK on these
// types has no effect -- i.e. the mutexes are not actually locked or unlocked
//
// Additionally, when a Mutex or RecursiveMutex is used with AutoLock
// and RecursiveAutoLock, respectively, these classes also suppressing
// the locking and unlocking of the mutex. Regardless of the build type,
// AutoLock and RecursiveAutoLock inherit from std::unique_lock<std::mutex>
// and std::unique_lock<std::recursive_mutex>, respectively. This means
// that in situations (such as is needed by the analysis category), the
// AutoLock and RecursiveAutoLock can be passed to functions requesting
// a std::unique_lock. Within these functions, since std::unique_lock
// member functions are not virtual, they will not retain the dummy locking
// and unlocking behavior
// --> An example of this behavior can be found in AutoLock.hh
//
//  Jonathan R. Madsen (February 21, 2018)
//

// global mutex types
typedef std::mutex           Mutex;
typedef std::recursive_mutex RecursiveMutex;

// mutex macros
#define MUTEX_INITIALIZER                                                                \
    {}
#define MUTEXINIT(mutex)                                                                 \
    ;                                                                                    \
    ;
#define MUTEXDESTROY(mutex)                                                              \
    ;                                                                                    \
    ;

// static functions: get_id(), sleep_for(...), sleep_until(...), yield(),
namespace ThisThread
{
using namespace std::this_thread;
}

// will be used in the future when migrating threading to task-based style
// and are currently used in unit tests
template <typename Tp>
using Promise = std::promise<Tp>;
template <typename Tp>
using Future = std::future<Tp>;
template <typename Tp>
using SharedFuture = std::shared_future<Tp>;

// Some useful types
typedef void* ThreadFunReturnType;
typedef void* ThreadFunArgType;
typedef int (*thread_lock)(Mutex*);
typedef int (*thread_unlock)(Mutex*);

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

// global thread types
using Thread       = std::thread;
using NativeThread = std::thread::native_handle_type;
// std::thread::id does not cast to integer
using Pid_t = std::thread::id;

// Condition
using Condition = std::condition_variable;

// Thread identifier
using ThreadId = Thread::id;

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
