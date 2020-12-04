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
// This file wraps a TBB task_group into a TaskGroup
//
// ---------------------------------------------------------------
// Author: Jonathan Madsen (Jun 21st 2018)
// ---------------------------------------------------------------

#pragma once

#include "PTL/TaskGroup.hh"
#include <functional>
#include <memory>

#if defined(PTL_USE_TBB)
#    include <tbb/tbb.h>
#endif

namespace PTL
{
class ThreadPool;

//--------------------------------------------------------------------------------------//
#if defined(PTL_USE_TBB)

class ThreadPool;
namespace
{
typedef ::tbb::task_group tbb_task_group_t;
}

//--------------------------------------------------------------------------------------//

template <typename Tp, typename Arg = Tp>
class TBBTaskGroup : public TaskGroup<Tp, Arg>
{
public:
    //------------------------------------------------------------------------//
    typedef decay_t<Arg>                                             ArgTp;
    typedef TBBTaskGroup<Tp, Arg>                                    this_type;
    typedef TaskGroup<Tp, Arg>                                       base_type;
    typedef typename base_type::result_type                          result_type;
    typedef typename base_type::packaged_task_type                   packaged_task_type;
    typedef typename base_type::future_type                          future_type;
    typedef typename base_type::promise_type                         promise_type;
    typedef typename JoinFunction<Tp, Arg>::Type                     join_type;
    typedef tbb::task_group                                          tbb_task_group_t;
    //------------------------------------------------------------------------//
    template <typename... Args>
    using task_type = Task<ArgTp, decay_t<Args>...>;
    //------------------------------------------------------------------------//

public:
    // Constructor
    template <typename Func>
    TBBTaskGroup(Func&& _join, ThreadPool* _tp = nullptr)
    : base_type(std::forward<Func>(_join), _tp)
    , m_tbb_task_group(new tbb_task_group_t)
    {}
    template <typename Up = Tp, enable_if_t<std::is_same<Up, void>::value, int> = 0>
    TBBTaskGroup(ThreadPool* _tp = nullptr)
    : base_type(_tp)
    , m_tbb_task_group(new tbb_task_group_t)
    {}

    // Destructor
    virtual ~TBBTaskGroup()
    {
        delete m_tbb_task_group;
        this->clear();
    }

    // delete copy-construct
    TBBTaskGroup(const this_type&) = delete;
    // define move-construct
    TBBTaskGroup(this_type&& rhs) = default;

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
    void run(Func&& func, Args... args)
    {
        auto _func = [=]() { return func(args...); };
        auto _task = wrap(std::move(_func));
        auto _lamb = [=]() { (*_task)(); };
        m_tbb_task_group->run(_lamb);
    }
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args>
    void exec(Func&& func, Args... args)
    {
        auto _func = [=]() { return func(args...); };
        auto _task = wrap(std::move(_func));
        auto _lamb = [=]() { (*_task)(); };
        m_tbb_task_group->run(_lamb);
    }
    //------------------------------------------------------------------------//
    template <typename Func, typename... Args, typename Up = Tp,
              enable_if_t<std::is_same<Up, void>::value, int> = 0>
    void parallel_for(uintmax_t nitr, uintmax_t, Func&& func, Args... args)
    {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nitr),
                          [&](const tbb::blocked_range<size_t>& range) {
                              for(size_t i = range.begin(); i != range.end(); ++i)
                                  func(args...);
                          });
    }

public:
    //------------------------------------------------------------------------//
    // this is not a native Tasking task group
    virtual bool is_native_task_group() const override { return false; }

    //------------------------------------------------------------------------//
    // wait on tbb::task_group, not internal thread-pool
    virtual void wait() override
    {
        base_type::wait();
        m_tbb_task_group->wait();
    }

public:
    using base_type::begin;
    using base_type::cbegin;
    using base_type::cend;
    using base_type::clear;
    using base_type::end;
    using base_type::get_tasks;

    //------------------------------------------------------------------------//
    // wait to finish
    template <typename Up = Tp, enable_if_t<!std::is_void<Up>::value, int> = 0>
    inline Up join(Up accum = {})
    {
        this->wait();
        for(auto& itr : m_task_set)
        {
            using RetT = decay_t<decltype(itr.get())>;
            accum      = m_join(std::ref(accum), std::forward<RetT>(itr.get()));
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

protected:
    // Protected variables
    tbb_task_group_t* m_tbb_task_group;
    using base_type:: operator++;
    using base_type:: operator--;
    using base_type::m_join;
    using base_type::m_task_set;
    using base_type::vtask_list;
};

//--------------------------------------------------------------------------------------//
#else
//--------------------------------------------------------------------------------------//

template <typename Tp, typename Arg = Tp>
using TBBTaskGroup = TaskGroup<Tp, Arg>;

//--------------------------------------------------------------------------------------//
#endif
//--------------------------------------------------------------------------------------//

}  // namespace PTL
