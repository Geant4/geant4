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
// Class Timer
//
// Class description:
//
// Class for timer objects, able to measure elasped user/system process time.
//
// Note: Uses <sys/times.h> & <unistd.h> - POSIX.1 defined
//       If used, this header must be included in the source (.cc) file and it
//       must be the first header file to be included!
//
// Member functions:
//
// Timer()
//   Construct a timer object
// Start()
//   Start timing
// Stop()
//   Stop timing
// bool IsValid()
//   Return true if have a valid time (ie start() and stop() called)
// double GetRealElapsed()
//   Return the elapsed real time between last calling start() and stop()
// double GetSystemElapsed()
//   Return the elapsed system time between last calling start() and stop()
// double GetUserElapsed()
//   Return the elapsed user time between last calling start() and stop()
//
// Operators:
//
// std::ostream& operator << (std::ostream& os, const Timer& t);
//   Print the elapsed real,system and usertimes on os. Prints **s for times
//   if !IsValid
//
// Member data:
//
// bool fValidTimes
//   True after start and stop have both been called more than once and
//   an equal number of times
// clock_t fStartRealTime,fEndRealTime
//   Real times (arbitrary time 0)
// tms fStartTimes,fEndTimes
//   Timing structures (see times(2)) for start and end times

// History:
// 23.08.96 P.Kent Updated to also computed real elapsed time
// 21.08.95 P.Kent
// 29.04.97 G.Cosmo Added timings for Windows/NT

#pragma once

#include "PTL/Types.hh"

#include <chrono>
#include <iomanip>
#include <sstream>

#if !(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
#    include <sys/times.h>
#    include <time.h>
#else
#    include <ctime>
#    define _SC_CLK_TCK 1

extern "C"
{
    int sysconf(int);
};

// Structure returned by times()

struct tms
{
    clock_t tms_utime;  /* user time */
    clock_t tms_stime;  /* system time */
    clock_t tms_cutime; /* user time, children */
    clock_t tms_cstime; /* system time, children */
};

extern "C"
{
    extern clock_t times(struct tms*);
};
#endif /* WIN32 */

namespace PTL
{
class Timer
{
public:
    PTL_DEFAULT_OBJECT(Timer)

public:
    void   Start();
    void   Stop();
    bool   IsValid() const;
    double GetRealElapsed() const;
    double GetSystemElapsed() const;
    double GetUserElapsed() const;

    static const char* GetClockTime();

private:
    bool fValidTimes{ false };
    using clock_type = std::chrono::high_resolution_clock;
    std::chrono::time_point<clock_type> fStartRealTime, fEndRealTime;
    tms                                 fStartTimes, fEndTimes;
};

// ------------------------------------------------------------
//      class inline implementation
// ------------------------------------------------------------

inline void
Timer::Start()
{
    fValidTimes = false;
    times(&fStartTimes);
    fStartRealTime = clock_type::now();
}

inline void
Timer::Stop()
{
    times(&fEndTimes);
    fEndRealTime = clock_type::now();
    fValidTimes  = true;
}

inline bool
Timer::IsValid() const
{
    return fValidTimes;
}

inline const char*
Timer::GetClockTime()
{
    time_t     rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    return asctime(timeinfo);
}

inline std::ostream&
operator<<(std::ostream& os, const Timer& t)
{
    // so fixed doesn't propagate
    std::stringstream ss;
    ss << std::fixed;
    if(t.IsValid())
    {
        ss << "Real=" << t.GetRealElapsed() << "s User=" << t.GetUserElapsed()
           << "s Sys=" << t.GetSystemElapsed() << "s";

        // avoid possible FPE error
        if(t.GetRealElapsed() > 1.0e-6)
        {
            double cpu_util =
                (t.GetUserElapsed() + t.GetSystemElapsed()) / t.GetRealElapsed() * 100.0;
            ss << std::setprecision(1);
            ss << " [Cpu=" << std::setprecision(1) << cpu_util << "%]";
        }
    }
    else
    {
        ss << "Real=****s User=****s Sys=****s";
    }
    os << ss.str();

    return os;
}

}  // namespace PTL