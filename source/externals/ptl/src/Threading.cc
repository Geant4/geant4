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
#include "PTL/AutoLock.hh"
#include "PTL/Globals.hh"

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#    include <Windows.h>
#else
#    include <sys/syscall.h>
#    include <sys/types.h>
#    include <unistd.h>
#endif

#include <atomic>

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
Threading::SetPinAffinity(int cpu, NativeThread& aT)
{
#if defined(__linux__) || defined(_AIX)
    cpu_set_t* aset = new cpu_set_t;
    CPU_ZERO(aset);
    CPU_SET(cpu, aset);
    pthread_t& _aT = static_cast<pthread_t&>(aT);
    return (pthread_setaffinity_np(_aT, sizeof(cpu_set_t), aset) == 0);
#else  // Not available for Mac, WIN,...
    ConsumeParameters(cpu, aT);
    return true;
#endif
}

//======================================================================================//
