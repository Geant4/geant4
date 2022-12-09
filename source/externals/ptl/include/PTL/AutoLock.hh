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
//
// ---------------------------------------------------------------
// Tasking class header file
//
/// Class Description:
///
/// This class provides a mechanism to create a mutex and locks/unlocks it.
/// Can be used by applications to implement in a portable way a mutexing logic.
/// Usage Example:
///
///      #include "Threading.hh"
///      #include "AutoLock.hh"
///
///      /// defined somewhere -- static so all threads see the same mutex
///      static Mutex aMutex;
///
///      /// somewhere else:
///      /// The AutoLock instance will automatically unlock the mutex when it
///      /// goes out of scope. One typically defines the scope within { } if
///      /// there is thread-safe code following the auto-lock
///
///      {
///          AutoLock l(&aMutex);
///          ProtectedCode();
///      }
///
///      UnprotectedCode();
///
///      /// When ProtectedCode() is calling a function that also tries to lock
///      /// a normal AutoLock + Mutex will "deadlock". In other words, the
///      /// the mutex in the ProtectedCode() function will wait forever to
///      /// acquire the lock that is being held by the function that called
///      /// ProtectedCode(). In this situation, use a RecursiveAutoLock +
///      /// RecursiveMutex, e.g.
///
///      /// defined somewhere -- static so all threads see the same mutex
///      static RecursiveMutex aRecursiveMutex;
///
///      /// this function is sometimes called directly and sometimes called
///      /// from SomeFunction_B(), which also locks the mutex
///      void SomeFunction_A()
///      {
///          /// when called from SomeFunction_B(), a Mutex + AutoLock will
///          /// deadlock
///          RecursiveAutoLock l(&aRecursiveMutex);
///          /// do something
///      }
///
///      void SomeFunction_B()
///      {
///
///          {
///              RecursiveAutoLock l(&aRecursiveMutex);
///              SomeFunction_A();
///          }
///
///          UnprotectedCode();
///      }
///
///
/// ---------------------------------------------------------------
/// Author: Andrea Dotti (15 Feb 2013): First Implementation
///
/// Update: Jonathan Madsen (9 Feb 2018): Replaced custom implementation
///      with inheritance from C++11 unique_lock, which inherits the
///      following member functions:
///
///      - unique_lock(unique_lock&& other) noexcept;
///      - explicit unique_lock(mutex_type& m);
///      - unique_lock(mutex_type& m, std::defer_lock_t t) noexcept;
///      - unique_lock(mutex_type& m, std::try_to_lock_t t);
///      - unique_lock(mutex_type& m, std::adopt_lock_t t);
///
///      - template <typename Rep, typename Period>
///        unique_lock(mutex_type& m,
///                   const std::chrono::duration<Rep,Period>&
///                   timeout_duration);
///
///      - template<typename Clock, typename Duration>
///        unique_lock(mutex_type& m,
///              const std::chrono::time_point<Clock,Duration>& timeout_time);
///
///      - void lock();
///      - void unlock();
///      - bool try_lock();
///
///      - template <typename Rep, typename Period>
///        bool try_lock_for(const std::chrono::duration<Rep,Period>&);
///
///      - template <typename Rep, typename Period>
///        bool try_lock_until(const std::chrono::time_point<Clock,Duration>&);
///
///      - void swap(unique_lock& other) noexcept;
///      - mutex_type* release() noexcept;
///      - mutex_type* mutex() const noexcept;
///      - bool owns_lock() const noexcept;
///      - explicit operator bool() const noexcept;
///      - unique_lock& operator=(unique_lock&& other);
///
/// ---------------------------------------------------------------
///
/// Note that AutoLock is defined also for a sequential Tasking build but below
/// regarding implementation (also found in Threading.hh)
///
///
///          NOTE ON Tasking SERIAL BUILDS AND MUTEX/UNIQUE_LOCK
///          ==================================================
///
/// Mutex and RecursiveMutex are always C++11 std::mutex types
/// however, in serial mode, using MUTEXLOCK and MUTEXUNLOCK on these
/// types has no effect -- i.e. the mutexes are not actually locked or unlocked
///
/// Additionally, when a Mutex or RecursiveMutex is used with AutoLock
/// and RecursiveAutoLock, respectively, these classes also suppressing
/// the locking and unlocking of the mutex. Regardless of the build type,
/// AutoLock and RecursiveAutoLock inherit from std::unique_lock<std::mutex>
/// and std::unique_lock<std::recursive_mutex>, respectively. This means
/// that in situations (such as is needed by the analysis category), the
/// AutoLock and RecursiveAutoLock can be passed to functions requesting
/// a std::unique_lock. Within these functions, since std::unique_lock
/// member functions are not virtual, they will not retain the dummy locking
/// and unlocking behavior
/// --> An example of this behavior can be found below
///
///  Jonathan R. Madsen (February 21, 2018)
///
/***

//======================================================================================//

typedef std::unique_lock<std::mutex> unique_lock_t;
// functions for casting AutoLock to std::unique_lock to demonstrate
// that AutoLock is NOT polymorphic
void as_unique_lock(unique_lock_t* lock) { lock->lock(); }
void as_unique_unlock(unique_lock_t* lock) { lock->unlock(); }

//======================================================================================//

void run(const uint64_t& n)
{
    // sync the threads a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    // get two mutexes to avoid deadlock when l32 actually locks
    AutoLock l32(TypeMutex<int32_t>(), std::defer_lock);
    AutoLock l64(TypeMutex<int64_t>(), std::defer_lock);

    // when serial: will not execute std::unique_lock::lock() because
    // it overrides the member function
    l32.lock();
    // regardless of serial or MT: will execute std::unique_lock::lock()
    // because std::unique_lock::lock() is not virtual
    as_unique_lock(&l64);

    std::cout << "Running iteration " << n << "..." << std::endl;
}

//======================================================================================//
// execute some work
template <typename thread_type = std::thread>
void exec(uint64_t n)
{
    // get two mutexes to avoid deadlock when l32 actually locks
    AutoLock l32(TypeMutex<int32_t>(), std::defer_lock);
    AutoLock l64(TypeMutex<int64_t>(), std::defer_lock);

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

//======================================================================================//

int main()
{
    print_threading();

    uint64_t n = 30;
    std::cout << "\nRunning with real threads...\n" << std::endl;
    exec<std::thread>(n);
    std::cout << "\nRunning with fake threads...\n" << std::endl;
    exec<DummyThread>(n);

}

***/

#pragma once

#include "PTL/Threading.hh"
#include "PTL/Utility.hh"

#include <chrono>
#include <iostream>
#include <mutex>
#include <system_error>

namespace PTL
{
// Note: Note that TemplateAutoLock by itself is not thread-safe and
//       cannot be shared among threads due to the locked switch
//
template <typename MutexT>
class TemplateAutoLock : public std::unique_lock<MutexT>
{
public:
    //------------------------------------------------------------------------//
    // Some useful typedefs
    //------------------------------------------------------------------------//
    using unique_lock_t = std::unique_lock<MutexT>;
    using this_type     = TemplateAutoLock<MutexT>;
    using mutex_type    = typename unique_lock_t::mutex_type;

public:
    //------------------------------------------------------------------------//
    // STL-consistent reference form constructors
    //------------------------------------------------------------------------//

    // reference form is consistent with STL lock_guard types
    // Locks the associated mutex by calling m.lock(). The behavior is
    // undefined if the current thread already owns the mutex except when
    // the mutex is recursive
    explicit TemplateAutoLock(mutex_type& _mutex)
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
    TemplateAutoLock(mutex_type&                               _mutex,
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
    TemplateAutoLock(mutex_type&                                     _mutex,
                     const std::chrono::time_point<Clock, Duration>& _timeout_time)
    : unique_lock_t(_mutex, std::defer_lock)
    {
        // call termination-safe locking. if serial, this call has no effect
        _lock_deferred(_timeout_time);
    }

    // Does not lock the associated mutex.
    TemplateAutoLock(mutex_type& _mutex, std::defer_lock_t _lock) noexcept
    : unique_lock_t(_mutex, _lock)
    {}

    // Tries to lock the associated mutex without blocking by calling
    // m.try_lock(). The behavior is undefined if the current thread already
    // owns the mutex except when the mutex is recursive.
    TemplateAutoLock(mutex_type& _mutex, std::try_to_lock_t _lock)
    : unique_lock_t(_mutex, _lock)
    {}

    // Assumes the calling thread already owns m
    TemplateAutoLock(mutex_type& _mutex, std::adopt_lock_t _lock)
    : unique_lock_t(_mutex, _lock)
    {}

public:
    //------------------------------------------------------------------------//
    // Backwards compatibility versions (constructor with pointer to mutex)
    //------------------------------------------------------------------------//
    TemplateAutoLock(mutex_type* _mutex)
    : unique_lock_t(*_mutex, std::defer_lock)
    {
        // call termination-safe locking. if serial, this call has no effect
        _lock_deferred();
    }

    TemplateAutoLock(mutex_type* _mutex, std::defer_lock_t _lock) noexcept
    : unique_lock_t(*_mutex, _lock)
    {}

    TemplateAutoLock(mutex_type* _mutex, std::try_to_lock_t _lock)
    : unique_lock_t(*_mutex, _lock)
    {}

    TemplateAutoLock(mutex_type* _mutex, std::adopt_lock_t _lock)
    : unique_lock_t(*_mutex, _lock)
    {}

private:
// helpful macros
#define _is_stand_mutex(Tp) (std::is_same<Tp, Mutex>::value)
#define _is_recur_mutex(Tp) (std::is_same<Tp, RecursiveMutex>::value)
#define _is_other_mutex(Tp) (!_is_stand_mutex(Tp) && !_is_recur_mutex(Tp))

    template <typename Tp                                             = MutexT,
              typename std::enable_if<_is_stand_mutex(Tp), int>::type = 0>
    std::string GetTypeString()
    {
        return "AutoLock<Mutex>";
    }

    template <typename Tp                                             = MutexT,
              typename std::enable_if<_is_recur_mutex(Tp), int>::type = 0>
    std::string GetTypeString()
    {
        return "AutoLock<RecursiveMutex>";
    }

    template <typename Tp                                             = MutexT,
              typename std::enable_if<_is_other_mutex(Tp), int>::type = 0>
    std::string GetTypeString()
    {
        return "AutoLock<UNKNOWN_MUTEX>";
    }

// pollution is bad
#undef _is_stand_mutex
#undef _is_recur_mutex
#undef _is_other_mutex

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
    // on Tasking object after some static mutex variables are deleted, leading
    // to the error code (typically on Clang compilers):
    //      libc++abi.dylib: terminating with uncaught exception of type
    //      std::__1::system_error: mutex lock failed: Invalid argument
    // this function protects against this failure until such a time that
    // these issues have been resolved

    //========================================================================//
    // standard locking
    inline void _lock_deferred()
    {
        try
        {
            this->unique_lock_t::lock();
        } catch(std::system_error& e)
        {
            PrintLockErrorMessage(e);
        }
    }

    //========================================================================//
    // Tries to lock the associated mutex by calling
    // m.try_lock_for(_timeout_duration). Blocks until specified
    // _timeout_duration has elapsed or the lock is acquired, whichever comes
    // first. May block for longer than _timeout_duration.
    template <typename Rep, typename Period>
    void _lock_deferred(const std::chrono::duration<Rep, Period>& _timeout_duration)
    {
        try
        {
            this->unique_lock_t::try_lock_for(_timeout_duration);
        } catch(std::system_error& e)
        {
            PrintLockErrorMessage(e);
        }
    }

    //========================================================================//
    // Tries to lock the associated mutex by calling
    // m.try_lock_until(_timeout_time). Blocks until specified _timeout_time has
    // been reached or the lock is acquired, whichever comes first. May block
    // for longer than until _timeout_time has been reached.
    template <typename Clock, typename Duration>
    void _lock_deferred(const std::chrono::time_point<Clock, Duration>& _timeout_time)
    {
        try
        {
            this->unique_lock_t::try_lock_until(_timeout_time);
        } catch(std::system_error& e)
        {
            PrintLockErrorMessage(e);
        }
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
#if defined(VERBOSE)
        cout << "Non-critical error: mutex lock failure in "
             << GetTypeString<mutex_type>() << ". "
             << "If the app is terminating, Tasking failed to "
             << "delete an allocated resource and a Tasking destructor is "
             << "being called after the statics were destroyed. \n\t--> "
             << "Exception: [code: " << e.code() << "] caught: " << e.what() << std::endl;
#else
        ConsumeParameters(e);
#endif
    }
};

// -------------------------------------------------------------------------- //
//
//      Use the non-template types below:
//          - AutoLock with Mutex
//          - RecursiveAutoLock with RecursiveMutex
//
// -------------------------------------------------------------------------- //

using AutoLock          = TemplateAutoLock<Mutex>;
using RecursiveAutoLock = TemplateAutoLock<RecursiveMutex>;

}  // namespace PTL
