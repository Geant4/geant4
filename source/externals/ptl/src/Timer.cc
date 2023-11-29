//
// MIT License
// Copyright (c) 2019 Jonathan R. Madsen
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
// ----------------------------------------------------------------------
// class Timer
//
// Implementation
// 29.04.97 G.Cosmo Added timings for Windows systems

#include "PTL/Timer.hh"

#include <stdexcept>

using namespace PTL;

#if !(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
#    include <unistd.h>
#else
#    include <sys/types.h>
#    include <windows.h>

//======================================================================================//

// extract milliseconds time unit
int
sysconf(int a)
{
    if(a == _SC_CLK_TCK)
        return 1000;
    else
        return 0;
}

//======================================================================================//

static clock_t
filetime2msec(FILETIME* t)
{
    return (clock_t)((((float) t->dwHighDateTime) * 429496.7296) +
                     (((float) t->dwLowDateTime) * .0001));
}

//======================================================================================//

clock_t
times(struct tms* t)
{
    FILETIME   ct = { 0, 0 }, et = { 0, 0 }, st = { 0, 0 }, ut = { 0, 0 }, rt = { 0, 0 };
    SYSTEMTIME realtime;

    GetSystemTime(&realtime);
    SystemTimeToFileTime(&realtime, &rt);  // get real time in 10^-9 sec
    if(t != 0)
    {
        GetProcessTimes(GetCurrentProcess(), &ct, &et, &st,
                        &ut);  // get process time in 10^-9 sec
        t->tms_utime = t->tms_cutime = filetime2msec(&ut);
        t->tms_stime = t->tms_cstime = filetime2msec(&st);
    }
    return filetime2msec(&rt);
}

//======================================================================================//

#endif /* WIN32 */

//======================================================================================//

double
Timer::GetRealElapsed() const
{
    if(!fValidTimes)
    {
        throw std::runtime_error("Timer::GetRealElapsed() - "
                                 "Timer not stopped or times not recorded!");
    }
    std::chrono::duration<double> diff = fEndRealTime - fStartRealTime;
    return diff.count();
}

//======================================================================================//

double
Timer::GetSystemElapsed() const
{
    if(!fValidTimes)
    {
        throw std::runtime_error("Timer::GetSystemElapsed() - "
                                 "Timer not stopped or times not recorded!");
    }
    double diff = fEndTimes.tms_stime - fStartTimes.tms_stime;
    return diff / sysconf(_SC_CLK_TCK);
}

//======================================================================================//

double
Timer::GetUserElapsed() const
{
    if(!fValidTimes)
    {
        throw std::runtime_error("Timer::GetUserElapsed() - "
                                 "Timer not stopped or times not recorded!");
    }
    double diff = fEndTimes.tms_utime - fStartTimes.tms_utime;
    return diff / sysconf(_SC_CLK_TCK);
}

//======================================================================================//
