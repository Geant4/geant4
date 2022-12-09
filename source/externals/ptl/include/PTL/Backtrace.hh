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
// Backtrace
//
// Description:
//
//  Prints backtraces after signals are caught. Available on Unix.
//
// Usage:
//  A standard set of signals are enabled by default:
//
//     SIGQUIT, SIGILL, SIGABRT, SIGKILL, SIGBUS, SIGSEGV
//
//  These should not interfere with debuggers and/or FPEDetection.
//  In order to turn off handling for one or more signals, one can do:
//
//    BackTrace::DefaultSignals() = std::set<int>{};
//    BackTrace::DefaultSignals() = std::set<int>{ SIGSEGV };
//
//  and so on, *before* creating the run-manager. After the run-manager
//  has been created, one should disable the signals:
//
//    BackTrace::Disable(BackTrace::DefaultSignals());
//
//  Additionally, at runtime, the environment variable "BACKTRACE" can
//  be set to select a specific set of signals or none, e.g. in bash:
//
//    export BACKTRACE="SIGQUIT,SIGSEGV"
//    export BACKTRACE="none"
//
//  The environment variable is case-insensitive and can use any of the
//  following delimiters: space, comma, semi-colon, colon
//
// Author: J.Madsen, 19 October 2020
// --------------------------------------------------------------------

#ifndef PTL_Backtrace_hh
#define PTL_Backtrace_hh 1

#include "Threading.hh"
#include "Types.hh"

#if defined(PTL_UNIX)
#    include <cxxabi.h>
#    include <execinfo.h>
#    include <unistd.h>
#endif

#if defined(PTL_LINUX)
#    include <features.h>
#endif

#include <cfenv>
#include <cmath>
#include <csignal>
#include <cstring>
#include <type_traits>

namespace PTL
{
template <typename FuncT, typename... ArgTypes>
#if __cpp_lib_is_invocable >= 201703
using ResultOf_t = std::invoke_result_t<FuncT, ArgTypes...>;
#else
using ResultOf_t = typename std::result_of<FuncT(ArgTypes...)>::type;
#endif
}  // namespace PTL

// compatible OS and compiler
#if defined(PTL_UNIX) &&                                                                 \
    (defined(__GNUC__) || defined(__clang__) || defined(_INTEL_COMPILER))
#    if !defined(PTL_SIGNAL_AVAILABLE)
#        define PTL_SIGNAL_AVAILABLE
#    endif
#    if !defined(PTL_DEMANGLE_AVAILABLE)
#        define PTL_DEMANGLE_AVAILABLE
#    endif
#endif

#if !defined(PTL_PSIGINFO_AVAILABLE)
#    if _XOPEN_SOURCE >= 700 || _POSIX_C_SOURCE >= 200809L
#        define PTL_PSIGINFO_AVAILABLE 1
#    else
#        define PTL_PSIGINFO_AVAILABLE 0
#    endif
#endif

//----------------------------------------------------------------------------//

namespace PTL
{
inline std::string
Demangle(const char* _str)
{
#if defined(PTL_DEMANGLE_AVAILABLE)
    // demangling a string when delimiting
    int   _status = 0;
    char* _ret    = ::abi::__cxa_demangle(_str, nullptr, nullptr, &_status);
    if(_ret && _status == 0)
        return std::string(const_cast<const char*>(_ret));
    return _str;
#else
    return _str;
#endif
}

//----------------------------------------------------------------------------//

inline std::string
Demangle(const std::string& _str)
{
    return Demangle(_str.c_str());
}

//----------------------------------------------------------------------------//

template <typename Tp>
inline std::string
Demangle()
{
    return Demangle(typeid(Tp).name());
}
}  // namespace PTL

//----------------------------------------------------------------------------//
//
//      ONLY IF SIGNAL_AVAILABLE
//
//----------------------------------------------------------------------------//
//
#if defined(PTL_SIGNAL_AVAILABLE)
//
//  these are not in the original POSIX.1-1990 standard so we are defining
//  them in case the OS hasn't
//  POSIX-1.2001
#    ifndef SIGTRAP
#        define SIGTRAP 5
#    endif
//  not specified in POSIX.1-2001, but nevertheless appears on most other
//  UNIX systems, where its default action is typically to terminate the
//  process with a core dump.
#    ifndef SIGEMT
#        define SIGEMT 7
#    endif
//  POSIX-1.2001
#    ifndef SIGURG
#        define SIGURG 16
#    endif
//  POSIX-1.2001
#    ifndef SIGXCPU
#        define SIGXCPU 24
#    endif
//  POSIX-1.2001
#    ifndef SIGXFSZ
#        define SIGXFSZ 25
#    endif
//  POSIX-1.2001
#    ifndef SIGVTALRM
#        define SIGVTALRM 26
#    endif
//  POSIX-1.2001
#    ifndef SIGPROF
#        define SIGPROF 27
#    endif
//  POSIX-1.2001
#    ifndef SIGINFO
#        define SIGINFO 29
#    endif

//----------------------------------------------------------------------------//

#    include <algorithm>
#    include <array>
#    include <cstdio>
#    include <cstdlib>
#    include <functional>
#    include <iomanip>
#    include <iostream>
#    include <map>
#    include <regex>
#    include <set>
#    include <sstream>
#    include <string>
#    include <tuple>
#    include <vector>

// PTL header
#    include "Threading.hh"

//----------------------------------------------------------------------------//

namespace PTL
{
class Backtrace
{
public:
    using sigaction_t   = struct sigaction;
    using exit_action_t = std::function<void(int)>;
    using frame_func_t  = std::function<std::string(const char*)>;
    using signal_set_t  = std::set<int>;

public:
    struct actions
    {
        using id_entry_t = std::tuple<std::string, int, std::string>;
        using id_list_t  = std::vector<id_entry_t>;

        std::map<int, bool>        is_active    = {};
        std::map<int, sigaction_t> current      = {};
        std::map<int, sigaction_t> previous     = {};
        std::vector<exit_action_t> exit_actions = {};
        const id_list_t            identifiers  = {
            id_entry_t("SIGHUP", SIGHUP, "terminal line hangup"),
            id_entry_t("SIGINT", SIGINT, "interrupt program"),
            id_entry_t("SIGQUIT", SIGQUIT, "quit program"),
            id_entry_t("SIGILL", SIGILL, "illegal instruction"),
            id_entry_t("SIGTRAP", SIGTRAP, "trace trap"),
            id_entry_t("SIGABRT", SIGABRT, "abort program (formerly SIGIOT)"),
            id_entry_t("SIGEMT", SIGEMT, "emulate instruction executed"),
            id_entry_t("SIGFPE", SIGFPE, "floating-point exception"),
            id_entry_t("SIGKILL", SIGKILL, "kill program"),
            id_entry_t("SIGBUS", SIGBUS, "bus error"),
            id_entry_t("SIGSEGV", SIGSEGV, "segmentation violation"),
            id_entry_t("SIGSYS", SIGSYS, "non-existent system call invoked"),
            id_entry_t("SIGPIPE", SIGPIPE, "write on a pipe with no reader"),
            id_entry_t("SIGALRM", SIGALRM, "real-time timer expired"),
            id_entry_t("SIGTERM", SIGTERM, "software termination signal"),
            id_entry_t("SIGURG", SIGURG, "urgent condition present on socket"),
            id_entry_t("SIGSTOP", SIGSTOP, "stop (cannot be caught or ignored)"),
            id_entry_t("SIGTSTP", SIGTSTP, "stop signal generated from keyboard"),
            id_entry_t("SIGCONT", SIGCONT, "continue after stop"),
            id_entry_t("SIGCHLD", SIGCHLD, "child status has changed"),
            id_entry_t("SIGTTIN", SIGTTIN,
                       "background read attempted from control terminal"),
            id_entry_t("SIGTTOU", SIGTTOU,
                       "background write attempted to control terminal"),
            id_entry_t("SIGIO ", SIGIO, "I/O is possible on a descriptor"),
            id_entry_t("SIGXCPU", SIGXCPU, "cpu time limit exceeded"),
            id_entry_t("SIGXFSZ", SIGXFSZ, "file size limit exceeded"),
            id_entry_t("SIGVTALRM", SIGVTALRM, "virtual time alarm"),
            id_entry_t("SIGPROF", SIGPROF, "profiling timer alarm"),
            id_entry_t("SIGWINCH", SIGWINCH, "Window size change"),
            id_entry_t("SIGINFO", SIGINFO, "status request from keyboard"),
            id_entry_t("SIGUSR1", SIGUSR1, "User defined signal 1"),
            id_entry_t("SIGUSR2", SIGUSR2, "User defined signal 2")
        };
    };

public:
    // a functor called for each frame in the backtrace
    static frame_func_t& FrameFunctor();
    // default set of signals
    static signal_set_t& DefaultSignals();
    // the signal handler
    static void Handler(int sig, siginfo_t* sinfo, void* context);
    // information message about the signal, performs exit-actions
    // and prints back-trace
    static void Message(int sig, siginfo_t* sinfo, std::ostream&);
    // calls user-provided functions after signal is caught but before abort
    static void ExitAction(int sig);
    // enable signals via a string (which is tokenized)
    static int Enable(const std::string&);
    // enable signals via set of integers, anything less than zero is ignored
    static int Enable(const signal_set_t& _signals = DefaultSignals());
    // disable signals
    static int Disable(signal_set_t _signals = {});
    // gets the numeric value for a signal name
    static int GetSignal(const std::string&);
    // provides a description of the signal
    static std::string Description(int sig);

    // adds an exit action
    template <typename FuncT>
    static void AddExitAction(FuncT&& func);

    // gets a backtrace of "Depth" frames. The offset parameter is used
    // to ignore initial frames (such as this function). A callback
    // can be provided to inspect and/or tweak the frame string
    template <size_t Depth, size_t Offset = 0, typename FuncT = frame_func_t>
    static std::array<ResultOf_t<FuncT, const char*>, Depth> GetMangled(
        FuncT&& func = FrameFunctor());

    // gets a demangled backtrace of "Depth" frames. The offset parameter is
    // used to ignore initial frames (such as this function). A callback
    // can be provided to inspect and/or tweak the frame string
    template <size_t Depth, size_t Offset = 0, typename FuncT = frame_func_t>
    static std::array<ResultOf_t<FuncT, const char*>, Depth> GetDemangled(
        FuncT&& func = FrameFunctor());

private:
    static actions& GetData()
    {
        static auto _instance = actions{};
        return _instance;
    }
};

//----------------------------------------------------------------------------//

// a functor called for each frame in the backtrace
inline Backtrace::frame_func_t&
Backtrace::FrameFunctor()
{
    static frame_func_t _instance = [](const char* inp) { return std::string(inp); };
    return _instance;
}

//----------------------------------------------------------------------------//

// default set of signals
inline Backtrace::signal_set_t&
Backtrace::DefaultSignals()
{
    static signal_set_t _instance = {
        SIGQUIT, SIGILL, SIGABRT, SIGKILL, SIGBUS, SIGSEGV
    };
    return _instance;
}

//----------------------------------------------------------------------------//

template <typename FuncT>
inline void
Backtrace::AddExitAction(FuncT&& func)
{
    GetData().exit_actions.emplace_back(std::forward<FuncT>(func));
}

//----------------------------------------------------------------------------//

inline void
Backtrace::ExitAction(int sig)
{
    for(auto& itr : GetData().exit_actions)
        itr(sig);
}

//----------------------------------------------------------------------------//

template <size_t Depth, size_t Offset, typename FuncT>
inline std::array<ResultOf_t<FuncT, const char*>, Depth>
Backtrace::GetMangled(FuncT&& func)
{
    static_assert((Depth - Offset) >= 1, "Error Depth - Offset should be >= 1");

    using type = ResultOf_t<FuncT, const char*>;
    // destination
    std::array<type, Depth> btrace;
    btrace.fill((std::is_pointer<type>::value) ? nullptr : type{});

    // plus one for this stack-frame
    std::array<void*, Depth + Offset> buffer;
    // size of returned buffer
    auto sz = backtrace(buffer.data(), Depth + Offset);
    // size of relevant data
    auto n = sz - Offset;

    // skip ahead (Offset + 1) stack frames
    char** bsym = backtrace_symbols(buffer.data() + Offset, n);

    // report errors
    if(bsym == nullptr)
        perror("backtrace_symbols");
    else
    {
        for(decltype(n) i = 0; i < n; ++i)
            btrace[i] = func(bsym[i]);
        free(bsym);
    }
    return btrace;
}

//----------------------------------------------------------------------------//

template <size_t Depth, size_t Offset, typename FuncT>
inline std::array<ResultOf_t<FuncT, const char*>, Depth>
Backtrace::GetDemangled(FuncT&& func)
{
    auto demangle_bt = [&](const char* cstr) {
        auto _trim = [](std::string& _sub, size_t& _len) {
            size_t _pos = 0;
            while((_pos = _sub.find_first_of(' ')) == 0)
            {
                _sub = _sub.erase(_pos, 1);
                --_len;
            }
            while((_pos = _sub.find_last_of(' ')) == _sub.length() - 1)
            {
                _sub = _sub.substr(0, _sub.length() - 1);
                --_len;
            }
            return _sub;
        };

        auto str = Demangle(std::string(cstr));
        auto beg = str.find("(");
        if(beg == std::string::npos)
        {
            beg = str.find("_Z");
            if(beg != std::string::npos)
                beg -= 1;
        }
        auto end = str.find("+", beg);
        if(beg != std::string::npos && end != std::string::npos)
        {
            auto len = end - (beg + 1);
            auto sub = str.substr(beg + 1, len);
            auto dem = Demangle(_trim(sub, len));
            str      = str.replace(beg + 1, len, dem);
        }
        else if(beg != std::string::npos)
        {
            auto len = str.length() - (beg + 1);
            auto sub = str.substr(beg + 1, len);
            auto dem = Demangle(_trim(sub, len));
            str      = str.replace(beg + 1, len, dem);
        }
        else if(end != std::string::npos)
        {
            auto len = end;
            auto sub = str.substr(beg, len);
            auto dem = Demangle(_trim(sub, len));
            str      = str.replace(beg, len, dem);
        }
        return func(str.c_str());
    };
    return GetMangled<Depth, Offset>(demangle_bt);
}

//----------------------------------------------------------------------------//

inline void
Backtrace::Message(int sig, siginfo_t* sinfo, std::ostream& os)
{
    // try to avoid as many dynamic allocations as possible here to avoid
    // overflowing the signal stack

    // ignore future signals of this type
    signal(sig, SIG_IGN);

    os << "\n### CAUGHT SIGNAL: " << sig << " ### ";
    if(sinfo)
        os << "address: " << sinfo->si_addr << ", ";
    os << Description(sig) << ". ";

    if(sig == SIGSEGV)
    {
        if(sinfo)
        {
            switch(sinfo->si_code)
            {
                case SEGV_MAPERR: os << "Address not mapped to object."; break;
                case SEGV_ACCERR: os << "Invalid permissions for mapped object."; break;
                default:
                    os << "Unknown segmentation fault error: " << sinfo->si_code << ".";
                    break;
            }
        }
        else
        {
            os << "Segmentation fault (unknown).";
        }
    }
    else if(sig == SIGFPE)
    {
        if(sinfo)
        {
            switch(sinfo->si_code)
            {
                case FE_DIVBYZERO: os << "Floating point divide by zero."; break;
                case FE_OVERFLOW: os << "Floating point overflow."; break;
                case FE_UNDERFLOW: os << "Floating point underflow."; break;
                case FE_INEXACT: os << "Floating point inexact result."; break;
                case FE_INVALID: os << "Floating point invalid operation."; break;
                default:
                    os << "Unknown floating point exception error: " << sinfo->si_code
                       << ".";
                    break;
            }
        }
        else
        {
            os << "Unknown floating point exception";
            if(sinfo)
                os << ": " << sinfo->si_code;
            os << ". ";
        }
    }

    os << '\n';

    auto bt = GetMangled<256, 3>([](const char* _s) { return _s; });
    char prefix[64];
    snprintf(prefix, 64, "[PID=%i, TID=%i]", (int) getpid(),
             (int) Threading::GetThreadId());
    size_t sz = 0;
    for(auto& itr : bt)
    {
        if(!itr)
            break;
        if(strlen(itr) == 0)
            break;
        ++sz;
    }
    os << "\nBacktrace:\n";
    auto _w = std::log10(sz) + 1;
    for(size_t i = 0; i < sz; ++i)
    {
        os << prefix << "[" << std::setw(_w) << std::right << i << '/' << std::setw(_w)
           << std::right << sz << "]> " << std::left << bt.at(i) << '\n';
    }
    os << std::flush;

    // exit action could cause more signals to be raise so make sure this is done
    // after the message has been printed
    try
    {
        ExitAction(sig);
    } catch(std::exception& e)
    {
        std::cerr << "ExitAction(" << sig << ") threw an exception" << std::endl;
        std::cerr << e.what() << std::endl;
    }
}

//----------------------------------------------------------------------------//

inline void
Backtrace::Handler(int sig, siginfo_t* sinfo, void*)
{
    Message(sig, sinfo, std::cerr);

    char msg[1024];
    snprintf(msg, 1024, "%s", "\n");

    if(sinfo && PTL_PSIGINFO_AVAILABLE > 0)
    {
#    if PTL_PSIGINFO_AVAILABLE > 0
        psiginfo(sinfo, msg);
        fflush(stdout);
        fflush(stderr);
#    endif
    }
    else
    {
        std::cerr << msg << std::flush;
    }

    // ignore any termination signals
    for(auto itr : { SIGKILL, SIGTERM, SIGABRT })
        signal(itr, SIG_IGN);
    abort();
}

//----------------------------------------------------------------------------//

inline int
Backtrace::Enable(const signal_set_t& _signals)
{
    static bool _first = true;
    if(_first)
    {
        std::string       _msg = "!!! Backtrace is activated !!!";
        std::stringstream _filler;
        std::stringstream _spacer;
        _filler.fill('#');
        _filler << std::setw(_msg.length()) << "";
        _spacer << std::setw(10) << "";
        std::cout << "\n\n"
                  << _spacer.str() << _filler.str() << "\n"
                  << _spacer.str() << _msg << "\n"
                  << _spacer.str() << _filler.str() << "\n\n"
                  << std::flush;
    }
    _first  = false;
    int cnt = 0;
    for(const auto& itr : _signals)
    {
        if(itr < 0)
            continue;
        if(GetData().is_active[itr])
            continue;
        ++cnt;
        sigfillset(&(GetData().current[itr].sa_mask));
        sigdelset(&(GetData().current[itr].sa_mask), itr);
        GetData().current[itr].sa_sigaction = &Handler;
        GetData().current[itr].sa_flags     = SA_SIGINFO;
        sigaction(itr, &(GetData().current[itr]), &(GetData().previous[itr]));
    }
    return cnt;
}

//----------------------------------------------------------------------------//

inline int
Backtrace::Enable(const std::string& _signals)
{
    if(_signals.empty())
        return 0;

    auto _add_signal = [](std::string sig, signal_set_t& _targ) {
        if(!sig.empty())
        {
            for(auto& itr : sig)
                itr = toupper(itr);
            _targ.insert(Backtrace::GetSignal(sig));
        }
    };

    const std::regex wsp_re("[ ,;:\t\n]+");
    auto             _maxid  = GetData().identifiers.size();
    auto             _result = std::vector<std::string>(_maxid, "");
    std::copy(std::sregex_token_iterator(_signals.begin(), _signals.end(), wsp_re, -1),
              std::sregex_token_iterator(), _result.begin());
    signal_set_t _sigset{};
    for(auto& itr : _result)
        _add_signal(itr, _sigset);
    return Enable(_sigset);
}

//----------------------------------------------------------------------------//

inline int
Backtrace::Disable(signal_set_t _signals)
{
    if(_signals.empty())
    {
        for(auto& itr : GetData().is_active)
            _signals.insert(itr.first);
    }

    int cnt = 0;
    for(const auto& itr : _signals)
    {
        if(itr < 0)
            continue;
        if(!GetData().is_active[itr])
            continue;
        ++cnt;
        sigaction(itr, &(GetData().previous[itr]), nullptr);
        GetData().current.erase(itr);
        GetData().is_active[itr] = false;
    }
    return cnt;
}

//----------------------------------------------------------------------------//

inline int
Backtrace::GetSignal(const std::string& sid)
{
    for(auto&& itr : GetData().identifiers)
    {
        if(std::get<0>(itr) == sid)
            return std::get<1>(itr);
    }
    return -1;
}

//----------------------------------------------------------------------------//

inline std::string
Backtrace::Description(int sig)
{
    for(auto&& itr : GetData().identifiers)
    {
        if(std::get<1>(itr) == sig)
        {
            std::stringstream ss;
            ss << " signal = " << std::setw(8) << std::get<0>(itr)
               << ", value = " << std::setw(4) << std::get<1>(itr)
               << ", description = " << std::get<2>(itr);
            return ss.str();
        }
    }
    std::stringstream ss;
    ss << " signal = " << std::setw(8) << "unknown"
       << ", value = " << std::setw(4) << sig;
    return ss.str();
}

//----------------------------------------------------------------------------//

}  // namespace PTL

#else

#    include <array>
#    include <functional>
#    include <map>
#    include <set>
#    include <string>
#    include <tuple>
#    include <vector>

namespace PTL
{
// dummy implementation
class Backtrace
{
public:
    struct fake_siginfo
    {};
    struct fake_sigaction
    {};

    using siginfo_t     = fake_siginfo;
    using sigaction_t   = fake_sigaction;
    using exit_action_t = std::function<void(int)>;
    using frame_func_t  = std::function<std::string(const char*)>;
    using signal_set_t  = std::set<int>;

public:
    struct actions
    {
        using id_entry_t = std::tuple<std::string, int, std::string>;
        using id_list_t  = std::vector<id_entry_t>;

        std::map<int, bool>        is_active    = {};
        std::map<int, sigaction_t> current      = {};
        std::map<int, sigaction_t> previous     = {};
        std::vector<exit_action_t> exit_actions = {};
        const id_list_t            identifiers  = {};
    };

public:
    static void        Handler(int, siginfo_t*, void*) {}
    static void        Message(int, siginfo_t*, std::ostream&) {}
    static void        ExitAction(int) {}
    static int         Enable(const std::string&) { return 0; }
    static int         Enable(const signal_set_t& = DefaultSignals()) { return 0; }
    static int         Disable(signal_set_t = {}) { return 0; }
    static int         GetSignal(const std::string&) { return -1; }
    static std::string Description(int) { return std::string{}; }

    template <typename FuncT>
    static void AddExitAction(FuncT&&)
    {}

    template <size_t Depth, size_t Offset = 0, typename FuncT = frame_func_t>
    static std::array<ResultOf_t<FuncT, const char*>, Depth> GetMangled(
        FuncT&& func = FrameFunctor())
    {
        using type = ResultOf_t<FuncT, const char*>;
        auto ret   = std::array<type, Depth>{};
        ret.fill(func(""));
        return ret;
    }

    template <size_t Depth, size_t Offset = 0, typename FuncT = frame_func_t>
    static std::array<ResultOf_t<FuncT, const char*>, Depth> GetDemangled(
        FuncT&& func = FrameFunctor())
    {
        using type = ResultOf_t<FuncT, const char*>;
        auto ret   = std::array<type, Depth>{};
        ret.fill(func(""));
        return ret;
    }

    // a functor called for each frame in the backtrace
    static frame_func_t& FrameFunctor()
    {
        static frame_func_t _instance = [](const char* _s) { return std::string(_s); };
        return _instance;
    }

    // default set of signals
    static signal_set_t& DefaultSignals()
    {
        static signal_set_t _instance = {};
        return _instance;
    }

    static actions& GetData()
    {
        static auto _instance = actions{};
        return _instance;
    }
};

//----------------------------------------------------------------------------//

}  // namespace PTL

#endif  // SIGNAL_AVAILABLE
#endif  // Backtrace_hh
