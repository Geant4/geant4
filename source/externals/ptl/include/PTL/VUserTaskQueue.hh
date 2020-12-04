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
//  ---------------------------------------------------------------
//  Tasking class header
//  Class Description:
//      Abstract base class for creating a task queue used by
//      ThreadPool
//  ---------------------------------------------------------------
//  Author: Jonathan Madsen
//  ---------------------------------------------------------------

#pragma once

#include "PTL/Globals.hh"
#include "PTL/Threading.hh"
#include "PTL/Types.hh"

#include <cstddef>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>

namespace PTL
{
class VTask;
class VTaskGroup;
class ThreadPool;
class ThreadData;

class VUserTaskQueue
{
public:
    typedef VTask*                task_pointer;
    typedef std::atomic<intmax_t> AtomicInt;
    typedef uintmax_t             size_type;
    typedef std::function<void()> function_type;
    typedef std::set<ThreadId>    ThreadIdSet;

public:
    // Constructor - accepting the number of workers
    explicit VUserTaskQueue(intmax_t nworkers = -1);
    // Virtual destructors are required by abstract classes
    // so add it by default, just in case
    virtual ~VUserTaskQueue();

public:
    // Virtual function for getting a task from the queue
    // parameters:
    //      1. int - get from specific sub-queue
    //      2. int - number of iterations
    // returns:
    //      VTask* - a task or nullptr
    virtual task_pointer GetTask(intmax_t subq = -1, intmax_t nitr = -1) = 0;

    // Virtual function for inserting a task into the queue
    // parameters:
    //      1. VTask* - task to insert
    //      2. int - sub-queue to inserting into
    // return:
    //      int - subqueue inserted into
    virtual intmax_t InsertTask(task_pointer, ThreadData* = nullptr,
                                intmax_t subq = -1) = 0;

    // Overload this function to hold threads
    virtual void     Wait()               = 0;
    virtual intmax_t GetThreadBin() const = 0;

    virtual void resize(intmax_t) = 0;

    // these are used for stanard checking
    virtual size_type size() const  = 0;
    virtual bool      empty() const = 0;

    virtual size_type bin_size(size_type bin) const  = 0;
    virtual bool      bin_empty(size_type bin) const = 0;

    // these are for slower checking, default to returning normal size()/empty
    virtual size_type true_size() const { return size(); }
    virtual bool      true_empty() const { return empty(); }

    // a method of executing a specific function on all threads
    virtual void ExecuteOnAllThreads(ThreadPool* tp, function_type f) = 0;

    virtual void ExecuteOnSpecificThreads(ThreadIdSet tid_set, ThreadPool* tp,
                                          function_type f) = 0;

    intmax_t workers() const { return m_workers; }

    virtual VUserTaskQueue* clone() = 0;

    // operator for number of tasks
    //      prefix versions
    // virtual uintmax_t operator++() = 0;
    // virtual uintmax_t operator--() = 0;
    //      postfix versions
    // virtual uintmax_t operator++(int) = 0;
    // virtual uintmax_t operator--(int) = 0;

public:
    template <typename ContainerT, size_t... Idx>
    static auto ContainerToTupleImpl(ContainerT&& container, mpl::index_sequence<Idx...>)
        -> decltype(std::make_tuple(std::forward<ContainerT>(container)[Idx]...))
    {
        return std::make_tuple(std::forward<ContainerT>(container)[Idx]...);
    }

    template <std::size_t N, typename ContainerT>
    static auto ContainerToTuple(ContainerT&& container)
        -> decltype(ContainerToTupleImpl(std::forward<ContainerT>(container),
                                         mpl::make_index_sequence<N>{}))
    {
        return ContainerToTupleImpl(std::forward<ContainerT>(container),
                                    mpl::make_index_sequence<N>{});
    }

    template <std::size_t N, std::size_t Nt, typename TupleT,
              enable_if_t<(N == Nt), int> = 0>
    static void TExecutor(TupleT&& _t)
    {
        if(std::get<N>(_t).get())
            (*(std::get<N>(_t)))();
    }

    template <std::size_t N, std::size_t Nt, typename TupleT,
              enable_if_t<(N < Nt), int> = 0>
    static void TExecutor(TupleT&& _t)
    {
        if(std::get<N>(_t).get())
            (*(std::get<N>(_t)))();
        TExecutor<N + 1, Nt, TupleT>(std::forward<TupleT>(_t));
    }

    template <typename TupleT, std::size_t N = std::tuple_size<decay_t<TupleT>>::value>
    static void Executor(TupleT&& __t)
    {
        TExecutor<0, N - 1, TupleT>(std::forward<TupleT>(__t));
    }

    template <typename Container,
              typename std::enable_if<std::is_same<Container, task_pointer>::value,
                                      int>::type = 0>
    static void Execute(Container& obj)
    {
        if(obj.get())
            (*obj)();
    }

    template <typename Container,
              typename std::enable_if<!std::is_same<Container, task_pointer>::value,
                                      int>::type = 0>
    static void Execute(Container& tasks)
    {
        /*
        for(auto& itr : tasks)
        {
            if(itr.get())
                (*itr)();
        }*/

        size_type n     = tasks.size();
        size_type max_n = 4;
        while(n > 0)
        {
            auto compute = (n > max_n) ? max_n : n;
            switch(compute)
            {
                case 4: {
                    auto t = ContainerToTuple<4>(tasks);
                    Executor(t);
                    break;
                }
                case 3: {
                    auto t = ContainerToTuple<3>(tasks);
                    Executor(t);
                    break;
                }
                case 2: {
                    auto t = ContainerToTuple<2>(tasks);
                    Executor(t);
                    break;
                }
                case 1: {
                    auto t = ContainerToTuple<1>(tasks);
                    Executor(t);
                    break;
                }
                case 0: break;
            }
            // tasks.erase(tasks.begin(), tasks.begin() + compute);
            n -= compute;
        }
    }

protected:
    intmax_t m_workers = 0;
};

}  // namespace PTL
