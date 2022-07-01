# Parallel Tasking Library (PTL)
Lightweight C++11 multithreading tasking system featuring thread-pool, task-groups, and lock-free task queue

## Basic Interface

```cpp
#include "PTL/PTL.hh"

#include <cassert>

inline long
fibonacci(long n)
{
    return (n < 2) ? n : (fibonacci(n - 1) + fibonacci(n - 2));
}

int main()
{
    bool use_tbb     = false;
    auto num_threads = 4;
    auto run_manager = PTL::TaskRunManager(use_tbb);

    run_manager.Initialize(num_threads);

    auto* task_manager = run_manager.GetTaskManager();

    // add a task via the task manager
    auto baz = task_manager->async<long>(fibonacci, 40);

    // functor to combine results
    auto join = [](long& lhs, long rhs) { return lhs += rhs; };

    // create a task group for 10 fibonacci calculations
    PTL::TaskGroup<long> foo(join);
    for(uint64_t i = 0; i < 10; ++i)
        foo.exec(fibonacci, 40);

    // create a task group for 10 fibonacci calculations
    PTL::TaskGroup<void> bar{};

    long ret_bar = 0;
    auto run     = [&ret_bar](long n) { ret_bar += fibonacci(n); };
    for(uint64_t i = 0; i < 10; ++i)
        bar.exec(run, 40);

    auto ret_baz = baz->get();
    auto ret_foo = foo.join();
    bar.join();

    assert(ret_baz * 10 == ret_foo);
    assert(ret_baz * 10 == ret_bar);
    assert(ret_foo == ret_bar);
}
```

## Explicit Thread-Pool

Using `PTL::TaskRunManager` is not necessary with task-groups.
You can create new thread-pools directly and pass them to task-groups:

```cpp
long example()
{
    // create a new thread-pool explicitly
    PTL::ThreadPool tp(4);

    // combines results
    auto join = [](long& lhs, long rhs) { return lhs += rhs; };

    // specify thread-pool explicitly
    PTL::TaskGroup<long> foo(join, &tp);

    for(int i = 0; i < 10; ++i)
        foo.exec(fibonacci, 40);

    // blocks until tasks in group are completed
    // thread-pool is destroyed after function returns
    return foo.get();
}
```
