//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4Threading
//
// Description:
//
// This unit defines types and macros used to expose Geant4 threading model.

// Author: Andrea Dotti, 15 February 2013 - First Implementation
// Revision: Jonathan R. Madsen, 21 February 2018
// --------------------------------------------------------------------
#ifndef G4Threading_hh
#define G4Threading_hh 1

#include "G4Types.hh"
#include "globals.hh"

#include <chrono>
#include <condition_variable>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

// Macro to put current thread to sleep
//
#define G4THREADSLEEP(tick)                                                    \
  std::this_thread::sleep_for(std::chrono::seconds(tick))

// Will be used in the future when migrating threading to task-based style
template <typename _Tp>
using G4Future = std::future<_Tp>;
template <typename _Tp>
using G4SharedFuture = std::shared_future<_Tp>;
template <typename _Tp>
using G4Promise = std::promise<_Tp>;

//          NOTE ON GEANT4 SERIAL BUILDS AND MUTEX/UNIQUE_LOCK
//          ==================================================
//
// G4Mutex and G4RecursiveMutex are always C++11 std::mutex types
// however, in serial mode, using G4MUTEXLOCK and G4MUTEXUNLOCK on these
// types has no effect -- i.e. the mutexes are not actually locked or unlocked
//
// Additionally, when a G4Mutex or G4RecursiveMutex is used with G4AutoLock
// and G4RecursiveAutoLock, respectively, these classes also suppressing
// the locking and unlocking of the mutex. Regardless of the build type,
// G4AutoLock and G4RecursiveAutoLock inherit from std::unique_lock<std::mutex>
// and std::unique_lock<std::recursive_mutex>, respectively. This means
// that in situations (such as is needed by the analysis category), the
// G4AutoLock and G4RecursiveAutoLock can be passed to functions requesting
// a std::unique_lock. Within these functions, since std::unique_lock
// member functions are not virtual, they will not retain the dummy locking
// and unlocking behavior
// --> An example of this behavior can be found in G4AutoLock.hh

// Global mutex types
using G4Mutex          = std::mutex;
using G4RecursiveMutex = std::recursive_mutex;

// Mutex macros
#define G4MUTEX_INITIALIZER                                                    \
  {}
#define G4MUTEXINIT(mutex)                                                     \
  ;                                                                            \
  ;
#define G4MUTEXDESTROY(mutex)                                                  \
  ;                                                                            \
  ;

// Static functions: get_id(), sleep_for(...), sleep_until(...), yield(),
namespace G4ThisThread
{
  using namespace std::this_thread;
}

// Will be used in the future when migrating threading to task-based style
// and are currently used in unit tests
template <typename _Tp>
using G4Promise = std::promise<_Tp>;
template <typename _Tp>
using G4Future = std::future<_Tp>;
template <typename _Tp>
using G4SharedFuture = std::shared_future<_Tp>;

// Some useful types
using G4ThreadFunReturnType = void*;
using G4ThreadFunArgType    = void*;
using thread_lock =
  G4int (*)(G4Mutex*);  // typedef G4int (*thread_lock)(G4Mutex*);
using thread_unlock =
  G4int (*)(G4Mutex*);  // typedef G4int (*thread_unlock)(G4Mutex*);

// Helper function for getting a unique static mutex for a specific
// class or type
// Usage example:
//   a template class "G4Cache<T>" that required a static
//   mutex for specific to type T:
//      G4AutoLock l(G4TypeMutex<G4Cache<T>>());
template <typename _Tp>
G4Mutex& G4TypeMutex(const unsigned int& _n = 0)
{
  static G4Mutex* _mutex = new G4Mutex();
  if(_n == 0)
    return *_mutex;

  static std::vector<G4Mutex*> _mutexes;
  if(_n > _mutexes.size())
    _mutexes.resize(_n, nullptr);
  if(!_mutexes[_n])
    _mutexes[_n] = new G4Mutex();
  return *(_mutexes[_n - 1]);
}

// Helper function for getting a unique static recursive_mutex for a
// specific class or type
// Usage example:
//                a template class "G4Cache<T>" that required a static
//                recursive_mutex for specific to type T:
//                        G4RecursiveAutoLock
//                        l(G4TypeRecursiveMutex<G4Cache<T>>());
template <typename _Tp>
G4RecursiveMutex& G4TypeRecursiveMutex(const unsigned int& _n = 0)
{
  static auto* _mutex = new G4RecursiveMutex();
  if(_n == 0)
    return *(_mutex);

  static std::vector<G4RecursiveMutex*> _mutexes;
  if(_n > _mutexes.size())
    _mutexes.resize(_n, nullptr);
  if(!_mutexes[_n])
    _mutexes[_n] = new G4RecursiveMutex();
  return *(_mutexes[_n - 1]);
}

#if defined(G4MULTITHREADED)
//==========================================
// G4MULTITHREADED is ON - threading enabled
//==========================================

// global thread types
using G4Thread       = std::thread;
using G4NativeThread = std::thread::native_handle_type;

// mutex macros
#  define G4MUTEXLOCK(mutex)                                                   \
    {                                                                          \
      (mutex)->lock();                                                         \
    }
#  define G4MUTEXUNLOCK(mutex)                                                 \
    {                                                                          \
      (mutex)->unlock();                                                       \
    }

// Macro to join thread
#  define G4THREADJOIN(worker) (worker).join()

// std::thread::id does not cast to integer
using G4Pid_t = std::thread::id;

// Instead of previous macro taking one argument, define function taking
// unlimited arguments
template <typename _Worker, typename _Func, typename... _Args>
void G4THREADCREATE(_Worker*& worker, _Func func, _Args... args)
{
  *worker = G4Thread(func, std::forward<_Args>(args)...);
}

// Conditions
//
// See G4MTRunManager for example on how to use these
//
using G4Condition = std::condition_variable;
#  define G4CONDITION_INITIALIZER                                              \
    {}
#  define G4CONDITIONWAIT(cond, lock) (cond)->wait(*lock);
#  define G4CONDITIONWAITLAMBDA(cond, lock, lambda) (cond)->wait(*lock, lambda);
#  define G4CONDITIONNOTIFY(cond) (cond)->notify_one();
#  define G4CONDITIONBROADCAST(cond) (cond)->notify_all();
//
// we don't define above globally so single-threaded code does not get
// caught in condition with no other thread to wake it up
//

#else
//==========================================
// G4MULTITHREADED is OFF - Sequential build
//==========================================

// implement a dummy thread class that acts like a thread
class G4DummyThread
{
 public:
  using native_handle_type = G4int;
  using id                 = std::thread::id;

 public:
  // does nothing
  G4DummyThread() {}
  // a std::thread-like constructor that execute upon construction
  template <typename _Func, typename... _Args>
  G4DummyThread(_Func func, _Args&&... _args)
  {
    func(std::forward<_Args>(_args)...);
  }

 public:
  native_handle_type native_handle() const { return native_handle_type(); }
  G4bool joinable() const { return true; }
  id get_id() const noexcept { return std::this_thread::get_id(); }
  void swap(G4DummyThread&) {}
  void join() {}
  void detach() {}

 public:
  static unsigned int hardware_concurrency() noexcept
  {
    return std::thread::hardware_concurrency();
  }
};

// global thread types
using G4Thread       = G4DummyThread;
using G4NativeThread = G4DummyThread::native_handle_type;

// mutex macros
#  define G4MUTEXLOCK(mutex)                                                   \
    ;                                                                          \
    ;
#  define G4MUTEXUNLOCK(mutex)                                                 \
    ;                                                                          \
    ;

// Macro to join thread
#  define G4THREADJOIN(worker)                                                 \
    ;                                                                          \
    ;

using G4Pid_t = G4int;

// Instead of previous macro taking one argument, define function taking
// unlimited arguments
template <typename _Worker, typename _Func, typename... _Args>
void G4THREADCREATE(_Worker*& worker, _Func func, _Args... args)
{
  *worker = G4Thread(func, std::forward<_Args>(args)...);
}

using G4Condition = G4int;
#  define G4CONDITION_INITIALIZER 1
#  define G4CONDITIONWAIT(cond, mutex) G4ConsumeParameters(cond, mutex);
#  define G4CONDITIONWAITLAMBDA(cond, mutex, lambda)                           \
    G4ConsumeParameters(cond, mutex, lambda);
#  define G4CONDITIONNOTIFY(cond) G4ConsumeParameters(cond);
#  define G4CONDITIONBROADCAST(cond) G4ConsumeParameters(cond);

#endif  // G4MULTITHREADING

//============================================================================//

// Define here after G4Thread has been typedef
using G4ThreadId = G4Thread::id;

//============================================================================//

namespace G4Threading
{
  enum
  {
    SEQUENTIAL_ID    = -2,
    MASTER_ID        = -1,
    WORKER_ID        = 0,
    GENERICTHREAD_ID = -1000
  };

  G4Pid_t G4GetPidId();
  G4int G4GetNumberOfCores();
  G4int G4GetThreadId();
  G4bool IsWorkerThread();
  G4bool IsMasterThread();
  void G4SetThreadId(G4int aNewValue);
  G4bool G4SetPinAffinity(G4int idx, G4NativeThread& at);
  void SetMultithreadedApplication(G4bool value);
  G4bool IsMultithreadedApplication();
  G4int WorkerThreadLeavesPool();
  G4int WorkerThreadJoinsPool();
  G4int GetNumberOfRunningWorkerThreads();
}  // namespace G4Threading

#endif  // G4Threading_hh
