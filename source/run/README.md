# Geant4 Tasking

This directory contains a Geant4 run manager which uses a tasking system for the G4Event loop.
This tasking system is fully compatible with TBB if `GEANT4_USE_TBB=ON` is specified when
configuring CMake. The default behavior, however, is to submit the tasks to an internal
thread-pool and task-queue.

## G4TaskRunManager

`G4TaskRunManager` multiply inherits from `G4MTRunManager` and `PTL::TaskRunManager`.
`PTL::TaskRunManager` holds the thread-pool instance, the size of the thread-pool,
and the default task-queue. The constructor of `G4TaskRunManager` takes a `G4VUserTaskQueue`
pointer (can be nullptr), a boolean for whether to use TBB if available, and a grainsize.

### Concepts

#### Grainsize

> Environment Variable: `G4FORCE_GRAINSIZE=N`

The grainsize is essentially the number of tasks. If set to 0, the default grainsize
will be `poolSize` and each thread will get `numEvents / poolSize` events.
If the grainsize is set to 1, then _all the events_ will be submitted as one task (i.e. be
processed serially by one thread in the pool). If the grainsize is set to 50 and there are 500 events,
then 50 tasks of 10 events will be submitted.

#### Events Per Tasks

> Environment Variable: `G4FORCE_EVENTS_PER_TASK=N`

Sometimes is easier to specify the number of events in a task instead of the grainsize.
If the events-per-task is set to 10 and there are 500 events,
then 50 tasks of 10 events will be submitted.

### Default Constructor

```cpp
    G4TaskRunManager(G4VUserTaskQueue* = nullptr, bool useTBB = false, G4int grainsize = 0);
```

## G4RunManagerFactory

An enumeration `G4RunManagerType` and a function `G4RunManagerFactory::CreateRunManager(...)`
was added to `"G4RunManagerFactory.hh"` to simplify the selection of the various run managers.
The first parameter is either one of the enumerated `G4RunManagerType` or a string identifier

| Enumeration                     | String ID   | Class               |
| ------------------------------- | ----------- | ------------------- |
| `G4RunManagerType::Serial`      | `"Serial"`  | `G4RunManager`      |
| `G4RunManagerType::MT`          | `"MT"`      | `G4MTRunManager`    |
| `G4RunManagerType::Tasking`     | `"Tasking"` | `G4TaskRunManager`  |
| `G4RunManagerType::TBB`         | `"TBB"`     | `G4TaskRunManager`  |
| `G4RunManagerType::Default`     | `"Default"` | Environment setting |
| `G4RunManagerType::SerialOnly`  | `"Serial"`  | `G4RunManager`      |
| `G4RunManagerType::MTOnly`      | `"MT"`      | `G4MTRunManager`    |
| `G4RunManagerType::TaskingOnly` | `"Tasking"` | `G4TaskRunManager`  |
| `G4RunManagerType::TBBOnly`     | `"TBB"`     | `G4TaskRunManager`  |


The `Default` enumeration value will defer to the following environment variable `G4RUN_MANAGER_TYPE`
if specified and will default to `"MT"` if MT is supported and serial if MT is not supported.
If the `G4FORCE_RUN_MANAGER_TYPE` environment variable is set, this variable will override the
value passed to the `CreateRunManager` function unless `G4RunManagerType` matches one of the `<TYPE>Only`
values. In this case, the environment variable is ignored and the run manager will be `<TYPE>`.

| Environment Variable       | Options                                  | Description                                                                            |
| -------------------------- | ---------------------------------------- | -------------------------------------------------------------------------------------- |
| `G4RUN_MANAGER_TYPE`       | `"Serial"`, `"MT"`, `"Tasking"`, `"TBB"` | Only applicable when `G4RunManagerType::Default` is used                               |
| `G4FORCE_RUN_MANAGER_TYPE` | `"Serial"`, `"MT"`, `"Tasking"`, `"TBB"` | Will override explicitly specifed `G4RunManagerType` if application allows and fail if type is not available |

## Creating the G4RunManager

- The `G4RunManagerFactory::CreateRunManager(...)` function takes either `G4RunManagerType` enumerated type or string to specify the desired G4RunManager
  - If a string is used, regex matching is used which is case-insensitive
  - Returns a `G4RunManager*`
  - Various overloads exist which just reorder passing in:
    - `int numberOfThreads` - executes `G4MTRunManager::SetNumberOfThreads(numberOfThreads)` before returning if > 0
      - default: `0`
    - `bool fail_if_unavail` - will cause a runtime failure if requested type is not available with Geant4 build
      - default: `true`
    - `G4VTaskQueue*` - a task-queue manager
      - default: `nullptr`

```cpp
#include "G4RunManagerFactory.hh"

int main()
{
    // specify {Serial, MT, Tasking, TBB} as the default, can be overridden
    // with "G4FORCE_RUN_MANAGER_TYPE" env variable
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Tasking);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::TBB);

    // specify {Serial, MT, Tasking, TBB} as the required type, cannot be overridden
    // with "G4FORCE_RUN_MANAGER_TYPE" env variable
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::TaskingOnly);
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::TBBOnly);

    // defer to "G4RUN_MANAGER_TYPE" env variable and default to MT if
    // env variable is not set
    auto* runmanager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // same as above
    auto* runmanager = G4RunManagerFactory::CreateRunManager();
}
```

## Using the Tasking System

With G4TaskRunManager, Geant4 events will be launched asynchronously as tasks. These tasks are
placed into a queue until one of the thread in the pool is available to execute the task. Users can
take advantage of this system to load-balance expensive sub-event calculations which might have
previously resulted in serial bottlenecks. For example, if an application needs to do extensive event
analysis on electrons and thread #1 ends up with 10x as many of these events, the other threads
might finish their G4Run significantly eariler and be idle while thread #1 has a lot of work.
Tasking allows these analysis calculations to be offload back into the queue so that other
threads can contribute to their completion.

### Option 1 - Submit Directly to Thread-Pool

- To execute the function `foo(int, double)` asynchronously:

```cpp
// get the task manager
auto* task_manager = G4TaskRunManager::GetTaskManager();

// submit task to thread-pool and receive a future for when the result is need
std::future<void> _fvoid = task_manager->async<void>(foo, 1, 1.0);
std::future<int>  _fint  = task_manager->async<int>(bar, 1.0);

// wait for task to execute
_fvoid.wait();
_fint.wait();

// get the result (if non-void)
auto result = _fint.get();
```

### Option 2 - Submit to task-group

- Obtain a pointer to the thread-pool instance
- Create a `task_group<T>` object where `T` is the return type of all the functions in the group
  - If `T` is non-void, you must provide a join functor who return type and first argument are both references
  to the joined type and the second argument is type `T`, e.g. `task_group<int>` can provide a join functor
  with `vector<int>&` as the return type and `T` as the second argument or `int&` as the return and first argument
  and `int` as the second argument
  - If `T` is void, the join functor is optional and can be treated as a final synchronization operation after
  all the tasks have been completed.

> NOTE: The join functor for task-groups are called sequentially on the thread that is
> waiting on `task_group<T>::join()` member function.

#### Global Definitions for Examples

```cpp
// obtain thread-pool instance from task manager
static auto* thread_pool = G4TaskRunManager::GetThreadPool();

// trivial int function which just returns value passed
int foo(int v) { return v; }

// function which launches CUDA kernel
void bar(int v)
{
    cuda_bar<<<512, 1>>>(v);
}
```

#### Example with non-void return types from tasks

```cpp
// put all return values from tasks into an array
auto join_vec = [](std::vector<int>& lhs, int rhs) { lhs.push_back(rhs); return lhs; };

// sum the values returned by tasks
auto sum_int = [](int& lhs, int rhs) { return lhs += rhs; };

// task group which applies 'join_vec' to all task return values
task_group<int>  vec_tg(join_vec, thread_pool);
// task group with applies 'sum_int' to all task return values
task_group<int>  sum_tg(sum_int, thread_pool);

// submit work to task-groups
vec_tg.exec(foo, 1);
vec_tg.exec(foo, 2);
sum_tg.exec(foo, 1);
sum_tg.exec(foo, 2);

// produces std::vector{ 1, 2 };
auto vec_result = vec_tg.join();

// produces 1 + 2 = 3
auto sum_result = sum_tg.join();
```

#### Example with void return type from tasks

```cpp
// wait for the GPU to finish
auto sync = []() { cudaDeviceSynchronize(); };

// task group which applies 'sync' after all tasks have been executed
task_group<void> gpu_tg(sync, thread_pool);
// generic task group w/o a join functor
task_group<void> general_tg(thread_pool);

// submit work to task-groups
gpu_tg.exec(bar, 1);
gpu_tg.exec(bar, 2);
general_tg.exec(bar, 1);
general_tg.exec(bar, 2);

// 'sync()' will get called after all tasks in group have executed
// (i.e. return from 'bar' function). 'sync' will then block until
// all GPU work has been completed
gpu_tg.join();

// will block only until all tasks in group have been executed
// (i.e. returned from 'bar' function)
generic_tg.join();

```
