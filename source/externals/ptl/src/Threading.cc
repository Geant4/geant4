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
//  Tasking class implementation
//
// Threading.cc
//

#include "PTL/Threading.hh"
#include "PTL/Types.hh"
#include "PTL/Utility.hh"

#if defined(PTL_WINDOWS)
#    include <Windows.h>
#endif

#if defined(PTL_MACOS)
#    include <sys/sysctl.h>
#endif

#if defined(PTL_LINUX)
#    include <fstream>
#endif

using namespace PTL;

//======================================================================================//

namespace
{
thread_local int ThreadID = Threading::MASTER_ID;

}  // namespace

//======================================================================================//

Pid_t
Threading::GetPidId()
{
    // In multithreaded mode return Thread ID
    return std::this_thread::get_id();
}

//======================================================================================//

unsigned
Threading::GetNumberOfCores()
{
    return std::thread::hardware_concurrency();
}

//======================================================================================//

unsigned
Threading::GetNumberOfPhysicalCpus()
{
#if defined(PTL_MACOS)
    int    count;
    size_t count_len = sizeof(count);
    sysctlbyname("hw.physicalcpu", &count, &count_len, nullptr, 0);
    return static_cast<unsigned>(count);
#elif defined(PTL_LINUX)
    unsigned      core_id_count = 0;
    std::ifstream ifs("/proc/cpuinfo");
    if(ifs)
    {
        std::set<std::string> core_ids;
        while(true)
        {
            std::string line = {};
            getline(ifs, line);
            if(!ifs.good())
                break;
            if(line.find("core id") != std::string::npos)
            {
                for(std::string itr : { "core id", ":", " ", "\t" })
                {
                    static auto _npos = std::string::npos;
                    auto        _pos  = _npos;
                    while((_pos = line.find(itr)) != _npos)
                        line = line.replace(_pos, itr.length(), "");
                }
                core_ids.insert(line);
            }
        }
        core_id_count = static_cast<unsigned>(core_ids.size());
        if(core_id_count > 0)
            return core_id_count;
    }
    return GetNumberOfCores();
#else
    return GetNumberOfCores();
#endif
}

//======================================================================================//

void
Threading::SetThreadId(int value)
{
    ThreadID = value;
}

int
Threading::GetThreadId()
{
    return ThreadID;
}

//======================================================================================//

bool
Threading::SetPinAffinity(int _cpu)
{
#if defined(__linux__) || defined(_AIX)
    cpu_set_t _cpu_set{};
    CPU_ZERO(&_cpu_set);
    if(_cpu >= 0)
        CPU_SET(_cpu, &_cpu_set);
    return (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &_cpu_set) == 0);
#else  // Not available for Mac, WIN,...
    ConsumeParameters(_cpu);
    return true;
#endif
}

//======================================================================================//

bool
Threading::SetThreadPriority(int _prio)
{
#if defined(__linux__) || defined(_AIX)
    return (pthread_setschedprio(pthread_self(), _prio) == 0);
#else  // Not available for Mac, WIN,...
    ConsumeParameters(_prio);
    return true;
#endif
}

//======================================================================================//

bool
Threading::SetPinAffinity(int _cpu, NativeThread& _t)
{
#if defined(__linux__) || defined(_AIX)
    cpu_set_t _cpu_set{};
    CPU_ZERO(&_cpu_set);
    CPU_SET(_cpu, &_cpu_set);
    return (pthread_setaffinity_np(static_cast<pthread_t>(_t), sizeof(cpu_set_t),
                                   &_cpu_set) == 0);
#else  // Not available for Mac, WIN,...
    ConsumeParameters(_cpu, _t);
    return true;
#endif
}

//======================================================================================//

bool
Threading::SetThreadPriority(int _prio, NativeThread& _t)
{
#if defined(__linux__) || defined(_AIX)
    return (pthread_setschedprio(static_cast<pthread_t>(_t), _prio) == 0);
#else
    ConsumeParameters(_prio, _t);
    return true;
#endif
}

//======================================================================================//
