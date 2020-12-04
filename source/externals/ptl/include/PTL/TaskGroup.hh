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
// Class Description:
//
// This file creates the a class for handling a group of tasks that
// can be independently joined
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Feb 13th 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/Task.hh"
#include "PTL/ThreadPool.hh"
#include "PTL/VTaskGroup.hh"

#include <cstdint>
#include <deque>
#include <future>
#include <list>
#include <vector>

#ifdef PTL_USE_TBB
#    include <tbb/tbb.h>
#endif

namespace PTL
{
class ThreadPool;

//--------------------------------------------------------------------------------------//

#if !defined(PTL_DEFAULT_OBJECT)
#    define PTL_DEFAULT_OBJECT(NAME)                                                     \
        NAME()                = default;                                                 \
        ~NAME()               = default;                                                 \
        NAME(const NAME&)     = default;                                                 \
        NAME(NAME&&) noexcept = default;                                                 \
        NAME& operator=(const NAME&) = default;                                          \
        NAME& operator=(NAME&&) noexcept = default;
#endif

//--------------------------------------------------------------------------------------//

template <typename JoinT, typename JoinArg>
struct JoinFunction
{
public:
    using Type = std::function<JoinT(JoinT&, JoinArg&&)>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    template <typename... Args>
    JoinT& operator()(Args&&... args)
    {
        return std::move(m_func(std::forward<Args>(args)...));
    }

private:
    Type m_func = [](JoinT& lhs, JoinArg&&) { return lhs; };
};

//--------------------------------------------------------------------------------------//

template <typename JoinArg>
struct JoinFunction<void, JoinArg>
{
public:
    using Type = std::function<void(JoinArg)>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    template <typename... Args>
    void operator()(Args&&... args)
    {
        m_func(std::forward<Args>(args)...);
    }

private:
    Type m_func = [](JoinArg) {};
};

//--------------------------------------------------------------------------------------//

template <>
struct JoinFunction<void, void>
{
public:
    using Type = std::function<void()>;

public:
    PTL_DEFAULT_OBJECT(JoinFunction)

    template <typename Func>
    JoinFunction(Func&& func)
    : m_func(std::forward<Func>(func))
    {}

    void operator()() { m_func(); }

private:
    Type m_func = []() {};
};

//--------------------------------------------------------------------------------------//

template <typename Tp, typename Arg = Tp>
class TaskGroup
: public VTaskGroup
, public TaskAllocator<TaskGroup<Tp, Arg>>
{
public:
    //------------------------------------------------------------------------//
    typedef decay_t<Arg>                                 ArgTp;
    typedef Tp                                           result_type;
    typedef TaskGroup<Tp, Arg>                           this_type;
    typedef std::promise<ArgTp>                          promise_type;
    typedef std::future<ArgTp>                           future_type;
    typedef std::packaged_task<ArgTp()>                  packaged_task_type;
    typedef list_type<future_type>                       task_list_t;
    typedef typename JoinFunction<Tp, Arg>::Type         join_type;
    typedef typename task_list_t::iterator               iterator;
    typedef typename task_list_t::reverse_iterator       reverse_iterator;
    typedef typename task_list_t::const_iterator         const_iterator;
    typedef typename task_list_t::const_reverse_iterator const_reverse_iterator;
    //------------------------------------------------------------------------//
    template <typename... Args>
    using task_type = Task<ArgTp, decay_t<Args>...>;
    //------------------------------------------------------------------------//

public:
    // Constructor
    template <typename Func>
    TaskGroup(Func&& _join, ThreadPool* _tp = nullptr)
    : VTaskGroup(_tp)
    , m_join(std::forward<Func>(_join))
    {}
    template <typename Up = Tp, enable_if_t<std::is_same<Up, void>::value, int> = 0>
    explicit TaskGroup(ThreadPool* _tp = nullptr)
    : VTaskGroup(_tp)
    , m_join([]() {})
    {}
    // Destructor
    virtual ~TaskGroup() { this->clear(); }

    // delete copy-construct
    TaskGroup(const this_type&) = delete;
    // define move-construct
    TaskGroup(this_type&& rhs) = default;
    // delete copy-assign
    this_type& operator=(const this_type& rhs) = delete;
    // define move-assign
    this_type& operator=(this_type&& rhs) = default;

public:
    //------------------------------------------------------------------------//
    template <typename Up>
    Up* operator+=(Up* _task)
    {
        // store in list
        vtask_list.push_back(_task);
        // thread-safe increment of tasks in task group
        operator++();
        // add the future
        m_task_set.push_back(std::move(_task->get_future()));
        // return
        return _task;
    }

public:
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args>
    task_type<Args...>* wrap(Func&& func, Args... args)
    {
        return operator+=(
            new task_type<Args...>(this, std::forward<Func>(func), args...));
    }

public:
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args>
    void exec(Func&& func, Args... args)
    {
        m_pool->add_task(wrap(std::forward<Func>(func), args...));
    }
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args>
    void run(Func&& func, Args... args)
    {
        m_pool->add_task(wrap(std::forward<Func>(func), args...));
    }
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args>
    void parallel_for(const intmax_t& nitr, const intmax_t& chunks, Func&& func,
                      Args... args)
    {
        auto nsplit = nitr / chunks;
        auto nmod   = nitr % chunks;
        if(nsplit < 1)
            nsplit = 1;
        for(intmax_t n = 0; n < nsplit; ++n)
        {
            auto _beg = n * chunks;
            auto _end = (n + 1) * chunks + ((n + 1 == nsplit) ? nmod : 0);
            run(std::forward<Func>(func), std::move(_beg), std::move(_end), args...);
        }
    }

protected:
    //------------------------------------------------------------------------//
    // shorter typedefs
    typedef iterator               itr_t;
    typedef const_iterator         citr_t;
    typedef reverse_iterator       ritr_t;
    typedef const_reverse_iterator critr_t;

public:
    //------------------------------------------------------------------------//
    // Get tasks with non-void return types
    //
    task_list_t&       get_tasks() { return m_task_set; }
    const task_list_t& get_tasks() const { return m_task_set; }

    //------------------------------------------------------------------------//
    // iterate over tasks with return type
    //
    itr_t   begin() { return m_task_set.begin(); }
    itr_t   end() { return m_task_set.end(); }
    citr_t  begin() const { return m_task_set.begin(); }
    citr_t  end() const { return m_task_set.end(); }
    citr_t  cbegin() const { return m_task_set.begin(); }
    citr_t  cend() const { return m_task_set.end(); }
    ritr_t  rbegin() { return m_task_set.rbegin(); }
    ritr_t  rend() { return m_task_set.rend(); }
    critr_t rbegin() const { return m_task_set.rbegin(); }
    critr_t rend() const { return m_task_set.rend(); }

    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, enable_if_t<!std::is_void<Up>::value, int> = 0>
    inline Up join(Up accum = {})
    {
        this->wait();
        for(auto& itr : m_task_set)
        {
            using RetT = decay_t<decltype(itr.get())>;
            accum = std::move(m_join(std::ref(accum), std::forward<RetT>(itr.get())));
        }
        this->clear();
        return accum;
    }
    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, typename Rp = Arg,
              enable_if_t<std::is_void<Up>::value && std::is_void<Rp>::value, int> = 0>
    inline void join()
    {
        this->wait();
        for(auto& itr : m_task_set)
            itr.get();
        m_join();
        this->clear();
    }
    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, typename Rp = Arg,
              enable_if_t<std::is_void<Up>::value && !std::is_void<Rp>::value, int> = 0>
    inline void join()
    {
        this->wait();
        for(auto& itr : m_task_set)
        {
            using RetT = decay_t<decltype(itr.get())>;
            m_join(std::forward<RetT>(itr.get()));
        }
        this->clear();
    }
    //------------------------------------------------------------------------//
    // clear the task result history
    void clear()
    {
        m_task_set.clear();
        VTaskGroup::clear();
    }

protected:
    // Protected variables
    task_list_t m_task_set;
    join_type   m_join;
};

}  // namespace PTL
