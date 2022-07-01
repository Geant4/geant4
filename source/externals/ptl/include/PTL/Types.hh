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
// Tasking native types
//

#pragma once

#ifndef G4GMAKE
#    include "PTL/Config.hh"
#endif

#if defined(__APPLE__) || defined(__MACH__)
#    if !defined(PTL_MACOS)
#        define PTL_MACOS 1
#    endif
#    if !defined(PTL_UNIX)
#        define PTL_UNIX 1
#    endif
#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#    if !defined(PTL_WINDOWS)
#        define PTL_WINDOWS 1
#    endif
#endif

#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#    if !defined(PTL_LINUX)
#        define PTL_LINUX 1
#    endif
#    if !defined(PTL_UNIX)
#        define PTL_UNIX 1
#    endif
#endif

#if defined(__unix__) || defined(__unix) || defined(unix)
#    if !defined(PTL_UNIX)
#        define PTL_UNIX 1
#    endif
#endif

#if defined(PTL_WINDOWS)
// Disable warning C4786 on WIN32 architectures:
// identifier was truncated to '255' characters
// in the debug information
//
#    pragma warning(disable : 4786)
//
// Define DLL export macro for WIN32 systems for
// importing/exporting external symbols to DLLs
//
#    if defined PTL_BUILD_DLL
#        define DLLEXPORT __declspec(dllexport)
#        define DLLIMPORT __declspec(dllimport)
#    else
#        define DLLEXPORT
#        define DLLIMPORT
#    endif
//
// Unique identifier for global module
//
#    if defined PTL_ALLOC_EXPORT
#        define PTL_DLL DLLEXPORT
#    else
#        define PTL_DLL DLLIMPORT
#    endif
#else
#    define DLLEXPORT
#    define DLLIMPORT
#    define PTL_DLL
#endif

#if !defined(PTL_DEFAULT_OBJECT)
#    define PTL_DEFAULT_OBJECT(NAME)                                                     \
        NAME()                       = default;                                          \
        ~NAME()                      = default;                                          \
        NAME(const NAME&)            = default;                                          \
        NAME(NAME&&)                 = default;                                          \
        NAME& operator=(const NAME&) = default;                                          \
        NAME& operator=(NAME&&)      = default;
#endif

#include <atomic>
#include <complex>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

namespace PTL
{
//--------------------------------------------------------------------------------------//
//
namespace api
{
struct native
{};
struct tbb
{};
}  // namespace api
//
//--------------------------------------------------------------------------------------//
//
template <typename Tp, typename Tag = api::native, typename Ptr = std::shared_ptr<Tp>,
          typename Pair = std::pair<Ptr, Ptr>>
Pair&
GetSharedPointerPair()
{
    static auto                 _master = std::make_shared<Tp>();
    static std::atomic<int64_t> _count(0);
    static thread_local auto    _inst =
        Pair(_master, Ptr((_count++ == 0) ? nullptr : new Tp()));
    return _inst;
}
//
//--------------------------------------------------------------------------------------//
//
template <typename Tp, typename Tag = api::native, typename Ptr = std::shared_ptr<Tp>,
          typename Pair = std::pair<Ptr, Ptr>>
Ptr
GetSharedPointerPairInstance()
{
    static thread_local auto& _pinst = GetSharedPointerPair<Tp, Tag>();
    static thread_local auto& _inst  = _pinst.second.get() ? _pinst.second : _pinst.first;
    return _inst;
}
//
//--------------------------------------------------------------------------------------//
//
template <typename Tp, typename Tag = api::native, typename Ptr = std::shared_ptr<Tp>,
          typename Pair = std::pair<Ptr, Ptr>>
Ptr
GetSharedPointerPairMasterInstance()
{
    static auto& _pinst = GetSharedPointerPair<Tp, Tag>();
    static auto  _inst  = _pinst.first;
    return _inst;
}

//======================================================================================//

}  // namespace PTL

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
class __void__;
