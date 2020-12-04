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
template <typename Tp>
Mutex&
TypeMutex(const unsigned int& _n = 0)
{
    static Mutex* _mutex = new Mutex();
    if(_n == 0)
        return *_mutex;

    static std::vector<Mutex*> _mutexes;
    if(_n > _mutexes.size())
        _mutexes.resize(_n, nullptr);
    if(!_mutexes[_n])
        _mutexes[_n] = new Mutex();
    return *(_mutexes[_n - 1]);
}

// Helper function for getting a unique static recursive_mutex for a
// specific class or type
// Usage example:
//		a template class "Cache<T>" that required a static
//		recursive_mutex for specific to type T:
//			RecursiveAutoLock l(TypeRecursiveMutex<Cache<T>>());
template <typename Tp>
RecursiveMutex&
TypeRecursiveMutex(const unsigned int& _n = 0)
{
    static RecursiveMutex* _mutex = new RecursiveMutex();
    if(_n == 0)
        return *(_mutex);

    static std::vector<RecursiveMutex*> _mutexes;
    if(_n > _mutexes.size())
        _mutexes.resize(_n, nullptr);
    if(!_mutexes[_n])
        _mutexes[_n] = new RecursiveMutex();
    return *(_mutexes[_n - 1]);
}

//======================================================================================//

// global thread types
typedef std::thread                     Thread;
typedef std::thread::native_handle_type NativeThread;

// mutex macros
#define MUTEXLOCK(mutex)                                                                 \
    {                                                                                    \
        (mutex)->lock();                                                                 \
    }
#define MUTEXUNLOCK(mutex)                                                               \
    {                                                                                    \
        (mutex)->unlock();                                                               \
    }

// Macro to join thread
#define THREADJOIN(worker) (worker).join()

// std::thread::id does not cast to integer
typedef std::thread::id Pid_t;

// Instead of previous macro taking one argument, define function taking
// unlimited arguments
template <typename WorkerT, typename FuncT, typename... Args>
void
THREADCREATE(WorkerT*& worker, FuncT func, Args... args)
{
    *worker = Thread(func, std::forward<Args>(args)...);
}

// Conditions
//
// See MTRunManager for example on how to use these
//
typedef std::condition_variable Condition;
#define CONDITION_INITIALIZER                                                            \
    {}
#define CONDITIONWAIT(cond, lock) (cond)->wait(*lock);
#define CONDITIONWAITLAMBDA(cond, lock, lambda) (cond)->wait(*lock, lambda);
#define CONDITIONNOTIFY(cond) (cond)->notify_one();
#define CONDITIONBROADCAST(cond) (cond)->notify_all();
//

//======================================================================================//

// Define here after Thread has been typedef
typedef Thread::id ThreadId;

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
GetNumberOfCores();
int
GetThreadId();
bool
IsWorkerThread();
bool
IsMasterThread();
void
SetThreadId(int aNewValue);
bool
SetPinAffinity(int idx, NativeThread& at);
int
WorkerThreadLeavesPool();
int
WorkerThreadJoinsPool();
int
GetNumberOfRunningWorkerThreads();
}  // namespace Threading

}  // namespace PTL
