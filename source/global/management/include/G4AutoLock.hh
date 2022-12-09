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
// G4Autolock
//
// Class Description:
//
// This class provides a mechanism to create a mutex and locks/unlocks it.
// Can be used by applications to implement in a portable way a mutexing logic.
// Usage Example:
//
//      #include "G4Threading.hh"
//      #include "G4AutoLock.hh"
//
//      // defined somewhere -- static so all threads see the same mutex
//      static G4Mutex aMutex;
//
//      // somewhere else:
//      // The G4AutoLock instance will automatically unlock the mutex when it
//      // goes out of scope. One typically defines the scope within { } if
//      // there is thread-safe code following the auto-lock
//
//      {
//          G4AutoLock l(&aMutex);
//          ProtectedCode();
//      }
//
//      UnprotectedCode();
//
//      // When ProtectedCode() is calling a function that also tries to lock
//      // a normal G4AutoLock + G4Mutex will "deadlock". In other words, the
//      // the mutex in the ProtectedCode() function will wait forever to
//      // acquire the lock that is being held by the function that called
//      // ProtectedCode(). In this situation, use a G4RecursiveAutoLock +
//      // G4RecursiveMutex, e.g.
//
//      // defined somewhere -- static so all threads see the same mutex
//      static G4RecursiveMutex aRecursiveMutex;
//
//      // this function is sometimes called directly and sometimes called
//      // from SomeFunction_B(), which also locks the mutex
//      void SomeFunction_A()
//      {
//          // when called from SomeFunction_B(), a G4Mutex + G4AutoLock will
//          // deadlock
//          G4RecursiveAutoLock l(&aRecursiveMutex);
//          // do something
//      }
//
//      void SomeFunction_B()
//      {
//
//          {
//              G4RecursiveAutoLock l(&aRecursiveMutex);
//              SomeFunction_A();
//          }
//
//          UnprotectedCode();
//      }
//

// --------------------------------------------------------------------
// Author: Andrea Dotti (15 Feb 2013): First Implementation
//
// Update: Jonathan Madsen (9 Feb 2018): Replaced custom implementation
//      with inheritance from C++11 unique_lock, which inherits the
//      following member functions:
//
//      - unique_lock(unique_lock&& other) noexcept;
//      - explicit unique_lock(mutex_type& m);
//      - unique_lock(mutex_type& m, std::defer_lock_t t) noexcept;
//      - unique_lock(mutex_type& m, std::try_to_lock_t t);
//      - unique_lock(mutex_type& m, std::adopt_lock_t t);
//
//      - template <typename Rep, typename Period>
//        unique_lock(mutex_type& m,
//                   const std::chrono::duration<Rep,Period>& timeout_duration);
//
//      - template<typename Clock, typename Duration>
//        unique_lock(mutex_type& m,
//              const std::chrono::time_point<Clock,Duration>& timeout_time);
//
//      - void lock();
//      - void unlock();
//      - bool try_lock();
//
//      - template <typename Rep, typename Period>
//        bool try_lock_for(const std::chrono::duration<Rep,Period>&);
//
//      - template <typename Rep, typename Period>
//        bool try_lock_until(const std::chrono::time_point<Clock,Duration>&);
//
//      - void swap(unique_lock& other) noexcept;
//      - mutex_type* release() noexcept;
//      - mutex_type* mutex() const noexcept;
//      - bool owns_lock() const noexcept;
//      - explicit operator bool() const noexcept;
//      - unique_lock& operator=(unique_lock&& other);
//
// --------------------------------------------------------------------
//
// Note that G4AutoLock is defined also for a sequential Geant4 build but below
// regarding implementation (also found in G4Threading.hh)
//
//
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
// --> An example of this behavior can be found below
//
//  Jonathan R. Madsen (February 21, 2018)
//
/**

//============================================================================//

void print_threading()
{
#ifdef G4MULTITHREADED
    std::cout << "\nUsing G4MULTITHREADED version..." << std::endl;
#else
    std::cout << "\nUsing G4SERIAL version..." << std::endl;
#endif
}

//============================================================================//

typedef std::unique_lock<std::mutex> unique_lock_t;
// functions for casting G4AutoLock to std::unique_lock to demonstrate
// that G4AutoLock is NOT polymorphic
void as_unique_lock(unique_lock_t* lock) { lock->lock(); }
void as_unique_unlock(unique_lock_t* lock) { lock->unlock(); }

//============================================================================//

void run(const uint64_t& n)
{
    // sync the threads a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // get two mutexes to avoid deadlock when l32 actually locks
    G4AutoLock l32(G4TypeMutex<int32_t>(), std::defer_lock);
    G4AutoLock l64(G4TypeMutex<int64_t>(), std::defer_lock);

    // when serial: will not execute std::unique_lock::lock() because
    // it overrides the member function
    l32.lock();
    // regardless of serial or MT: will execute std::unique_lock::lock()
    // because std::unique_lock::lock() is not virtual
    as_unique_lock(&l64);

    std::cout << "Running iteration " << n << "..." << std::endl;
}

//============================================================================//
// execute some work
template <typename thread_type = std::thread>
void exec(uint64_t n)
{
    // get two mutexes to avoid deadlock when l32 actually locks
    G4AutoLock l32(G4TypeMutex<int32_t>(), std::defer_lock);
    G4AutoLock l64(G4TypeMutex<int64_t>(), std::defer_lock);

    std::vector<thread_type*> threads(n, nullptr);
    for(uint64_t i = 0; i < n; ++i)
    {
        threads[i] = new thread_type();
        *(threads[i]) = std::move(thread_type(run, i));
    }

    // when serial: will not execute std::unique_lock::lock() because
    // it overrides the member function
    l32.lock();
    // regardless of serial or MT: will execute std::unique_lock::lock()
    // because std::unique_lock::lock() is not virtual
    as_unique_lock(&l64);

    std::cout << "Joining..." << std::endl;

    // when serial: will not execute std::unique_lock::unlock() because
    // it overrides the member function
    l32.unlock();
    // regardless of serial or MT: will execute std::unique_lock::unlock()
    // because std::unique_lock::unlock() is not virtual
    as_unique_unlock(&l64);

    // NOTE ABOUT UNLOCKS:
    // in MT, commenting out either
    //      l32.unlock();
    // or
    //      as_unique_unlock(&l64);
    // creates a deadlock; in serial, commenting out
    //      as_unique_unlock(&l64);
    // creates a deadlock but commenting out
    //      l32.unlock();
    // does not

    // clean up and join
    for(uint64_t i = 0; i < n; ++i)
    {
        threads[i]->join();
        delete threads[i];
    }
    threads.clear();
}

//============================================================================//

int main()
{
    print_threading();

    uint64_t n = 30;
    std::cout << "\nRunning with real threads...\n" << std::endl;
    exec<std::thread>(n);
    std::cout << "\nRunning with fake threads...\n" << std::endl;
    exec<G4DummyThread>(n);

}

**/
// --------------------------------------------------------------------
#ifndef G4AUTOLOCK_HH
#define G4AUTOLOCK_HH

#include "G4Threading.hh"

#include <chrono>
#include <iostream>
#include <mutex>
#include <system_error>

// Note: Note that G4TemplateAutoLock by itself is not thread-safe and
//       cannot be shared among threads due to the locked switch
//
template <typename _Mutex_t>
class G4TemplateAutoLock : public std::unique_lock<_Mutex_t>
{
 public:
  //------------------------------------------------------------------------//
  // Some useful typedefs
  //------------------------------------------------------------------------//
  using unique_lock_t = std::unique_lock<_Mutex_t>;
  using this_type = G4TemplateAutoLock<_Mutex_t>;
  using mutex_type = typename unique_lock_t::mutex_type;

 public:
  //------------------------------------------------------------------------//
  // STL-consistent reference form constructors
  //------------------------------------------------------------------------//

  // reference form is consistent with STL lock_guard types
  // Locks the associated mutex by calling m.lock(). The behavior is
  // undefined if the current thread already owns the mutex except when
  // the mutex is recursive
  G4TemplateAutoLock(mutex_type& _mutex)
    : unique_lock_t(_mutex, std::defer_lock)
  {
    // call termination-safe locking. if serial, this call has no effect
    _lock_deferred();
  }

  // Tries to lock the associated mutex by calling
  // m.try_lock_for(_timeout_duration). Blocks until specified
  // _timeout_duration has elapsed or the lock is acquired, whichever comes
  // first. May block for longer than _timeout_duration.
  template <typename Rep, typename Period>
  G4TemplateAutoLock(
    mutex_type& _mutex,
    const std::chrono::duration<Rep, Period>& _timeout_duration)
    : unique_lock_t(_mutex, std::defer_lock)
  {
    // call termination-safe locking. if serial, this call has no effect
    _lock_deferred(_timeout_duration);
  }

  // Tries to lock the associated mutex by calling
  // m.try_lock_until(_timeout_time). Blocks until specified _timeout_time has
  // been reached or the lock is acquired, whichever comes first. May block
  // for longer than until _timeout_time has been reached.
  template <typename Clock, typename Duration>
  G4TemplateAutoLock(
    mutex_type& _mutex,
    const std::chrono::time_point<Clock, Duration>& _timeout_time)
    : unique_lock_t(_mutex, std::defer_lock)
  {
    // call termination-safe locking. if serial, this call has no effect
    _lock_deferred(_timeout_time);
  }

  // Does not lock the associated mutex.
  G4TemplateAutoLock(mutex_type& _mutex, std::defer_lock_t _lock) noexcept
    : unique_lock_t(_mutex, _lock)
  {}

#ifdef G4MULTITHREADED

  // Tries to lock the associated mutex without blocking by calling
  // m.try_lock(). The behavior is undefined if the current thread already
  // owns the mutex except when the mutex is recursive.
  G4TemplateAutoLock(mutex_type& _mutex, std::try_to_lock_t _lock)
    : unique_lock_t(_mutex, _lock)
  {}

  // Assumes the calling thread already owns m
  G4TemplateAutoLock(mutex_type& _mutex, std::adopt_lock_t _lock)
    : unique_lock_t(_mutex, _lock)
  {}

#else

  // serial dummy version (initializes unique_lock but does not lock)
  G4TemplateAutoLock(mutex_type& _mutex, std::try_to_lock_t)
    : unique_lock_t(_mutex, std::defer_lock)
  {}

  // serial dummy version (initializes unique_lock but does not lock)
  G4TemplateAutoLock(mutex_type& _mutex, std::adopt_lock_t)
    : unique_lock_t(_mutex, std::defer_lock)
  {}

#endif  // defined(G4MULTITHREADED)

 public:
  //------------------------------------------------------------------------//
  // Backwards compatibility versions (constructor with pointer to mutex)
  //------------------------------------------------------------------------//
  G4TemplateAutoLock(mutex_type* _mutex)
    : unique_lock_t(*_mutex, std::defer_lock)
  {
    // call termination-safe locking. if serial, this call has no effect
    _lock_deferred();
  }

  G4TemplateAutoLock(mutex_type* _mutex, std::defer_lock_t _lock) noexcept
    : unique_lock_t(*_mutex, _lock)
  {}

#if defined(G4MULTITHREADED)

  G4TemplateAutoLock(mutex_type* _mutex, std::try_to_lock_t _lock)
    : unique_lock_t(*_mutex, _lock)
  {}

  G4TemplateAutoLock(mutex_type* _mutex, std::adopt_lock_t _lock)
    : unique_lock_t(*_mutex, _lock)
  {}

#else  // NOT defined(G4MULTITHREADED) -- i.e. serial

  G4TemplateAutoLock(mutex_type* _mutex, std::try_to_lock_t)
    : unique_lock_t(*_mutex, std::defer_lock)
  {}

  G4TemplateAutoLock(mutex_type* _mutex, std::adopt_lock_t)
    : unique_lock_t(*_mutex, std::defer_lock)
  {}

#endif  // defined(G4MULTITHREADED)

 public:
  //------------------------------------------------------------------------//
  // Non-constructor overloads
  //------------------------------------------------------------------------//

#if defined(G4MULTITHREADED)

  // overload nothing

#else  // NOT defined(G4MULTITHREADED) -- i.e. serial

  // override unique lock member functions to keep from locking/unlocking
  // but does not override in polymorphic usage
  void lock() {}
  void unlock() {}
  bool try_lock() { return true; }

  template <typename Rep, typename Period>
  bool try_lock_for(const std::chrono::duration<Rep, Period>&)
  {
    return true;
  }

  template <typename Clock, typename Duration>
  bool try_lock_until(const std::chrono::time_point<Clock, Duration>&)
  {
    return true;
  }

  void swap(this_type& other) noexcept { std::swap(*this, other); }
  bool owns_lock() const noexcept { return false; }

  // no need to overload
  // explicit operator bool() const noexcept;
  // this_type& operator=(this_type&& other);
  // mutex_type* release() noexcept;
  // mutex_type* mutex() const noexcept;

#endif  // defined(G4MULTITHREADED)

 private:
// helpful macros
#define _is_stand_mutex(_Tp) (std::is_same<_Tp, G4Mutex>::value)
#define _is_recur_mutex(_Tp) (std::is_same<_Tp, G4RecursiveMutex>::value)
#define _is_other_mutex(_Tp) (!_is_stand_mutex(_Tp) && !_is_recur_mutex(_Tp))

  template <typename _Tp                                             = _Mutex_t,
            typename std::enable_if<_is_stand_mutex(_Tp), int>::type = 0>
  std::string GetTypeString()
  {
    return "G4AutoLock<G4Mutex>";
  }

  template <typename _Tp                                             = _Mutex_t,
            typename std::enable_if<_is_recur_mutex(_Tp), int>::type = 0>
  std::string GetTypeString()
  {
    return "G4AutoLock<G4RecursiveMutex>";
  }

  template <typename _Tp                                             = _Mutex_t,
            typename std::enable_if<_is_other_mutex(_Tp), int>::type = 0>
  std::string GetTypeString()
  {
    return "G4AutoLock<UNKNOWN_MUTEX>";
  }

// pollution is bad
#undef _is_stand_mutex
#undef _is_recur_mutex
#undef _is_other_mutex

  // used in _lock_deferred chrono variants to avoid ununsed-variable warning
  template <typename _Tp>
  void suppress_unused_variable(const _Tp&)
  {}

  //========================================================================//
  // NOTE on _lock_deferred(...) variants:
  //      a system_error in lock means that the mutex is unavailable
  //      we want to throw the error that comes from locking an unavailable
  //      mutex so that we know there is a memory leak
  //      if the mutex is valid, this will hold until the other thread
  //      finishes

  // sometimes certain destructors use locks, this isn't an issue unless
  // the object is leaked. When this occurs, the application finalization
  // (i.e. the real or implied "return 0" part of main) will call destructors
  // on Geant4 object after some static mutex variables are deleted, leading
  // to the error code (typically on Clang compilers):
  //      libc++abi.dylib: terminating with uncaught exception of type
  //      std::__1::system_error: mutex lock failed: Invalid argument
  // this function protects against this failure until such a time that
  // these issues have been resolved

  //========================================================================//
  // standard locking
  inline void _lock_deferred()
  {
#if defined(G4MULTITHREADED)
    try
    {
      this->unique_lock_t::lock();
    } catch(std::system_error& e)
    {
      PrintLockErrorMessage(e);
    }
#endif
  }

  //========================================================================//
  // Tries to lock the associated mutex by calling
  // m.try_lock_for(_timeout_duration). Blocks until specified
  // _timeout_duration has elapsed or the lock is acquired, whichever comes
  // first. May block for longer than _timeout_duration.
  template <typename Rep, typename Period>
  void _lock_deferred(
    const std::chrono::duration<Rep, Period>& _timeout_duration)
  {
#if defined(G4MULTITHREADED)
    try
    {
      this->unique_lock_t::try_lock_for(_timeout_duration);
    } catch(std::system_error& e)
    {
      PrintLockErrorMessage(e);
    }
#else
    suppress_unused_variable(_timeout_duration);
#endif
  }

  //========================================================================//
  // Tries to lock the associated mutex by calling
  // m.try_lock_until(_timeout_time). Blocks until specified _timeout_time has
  // been reached or the lock is acquired, whichever comes first. May block
  // for longer than until _timeout_time has been reached.
  template <typename Clock, typename Duration>
  void _lock_deferred(
    const std::chrono::time_point<Clock, Duration>& _timeout_time)
  {
#if defined(G4MULTITHREADED)
    try
    {
      this->unique_lock_t::try_lock_until(_timeout_time);
    } catch(std::system_error& e)
    {
      PrintLockErrorMessage(e);
    }
#else
    suppress_unused_variable(_timeout_time);
#endif
  }

  //========================================================================//
  // the message for what mutex lock fails due to deleted static mutex
  // at termination
  void PrintLockErrorMessage(std::system_error& e)
  {
    // use std::cout/std::endl to avoid include dependencies
    using std::cout;
    using std::endl;
    // the error that comes from locking an unavailable mutex
#if defined(G4VERBOSE)
    cout << "Non-critical error: mutex lock failure in "
         << GetTypeString<mutex_type>() << ". "
         << "If the app is terminating, Geant4 failed to "
         << "delete an allocated resource and a Geant4 destructor is "
         << "being called after the statics were destroyed. \n\t--> "
         << "Exception: [code: " << e.code() << "] caught: " << e.what()
         << endl;
#else
    suppress_unused_variable(e);
#endif
  }
};

// -------------------------------------------------------------------------- //
//
//      Use the non-template types below:
//          - G4AutoLock with G4Mutex
//          - G4RecursiveAutoLock with G4RecursiveMutex
//
// -------------------------------------------------------------------------- //

using G4AutoLock          = G4TemplateAutoLock<G4Mutex>;
using G4RecursiveAutoLock = G4TemplateAutoLock<G4RecursiveMutex>;

// provide abbriviated type if another mutex type is desired to be used
// aside from above
template <typename _Tp>
using G4TAutoLock = G4TemplateAutoLock<_Tp>;

#endif  // G4AUTOLOCK_HH
