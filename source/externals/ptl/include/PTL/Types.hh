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

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
// Disable warning C4786 on WIN32 architectures:
// identifier was truncated to '255' characters
// in the debug information
//
#    pragma warning(disable : 4786)
//
// Define DLL export macro for WIN32 systems for
// importing/exporting external symbols to DLLs
//
#    if defined G4LIB_BUILD_DLL
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

#include <atomic>
#include <complex>
#include <limits>
#include <memory>

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

//--------------------------------------------------------------------------------------//

template <typename CountedType>
class CountedObject
{
public:
    typedef CountedObject<CountedType> this_type;
    typedef CountedObject<void>        void_type;

public:
    // return number of existing objects:
    static int64_t           live() { return count(); }
    static constexpr int64_t zero() { return static_cast<int64_t>(0); }
    static int64_t           max_depth() { return fmax_depth; }

    static void enable(const bool& val) { fenabled = val; }
    static void set_max_depth(const int64_t& val) { fmax_depth = val; }
    static bool is_enabled() { return fenabled; }

    template <typename Tp                                                   = CountedType,
              typename std::enable_if<std::is_same<Tp, void>::value>::type* = nullptr>
    static bool enable()
    {
        return fenabled && fmax_depth > count();
    }
    // the void type is consider the global setting
    template <typename Tp = CountedType,
              typename std::enable_if<!std::is_same<Tp, void>::value>::type* = nullptr>
    static bool enable()
    {
        return void_type::is_enabled() && void_type::max_depth() > count() && fenabled &&
               fmax_depth > count();
    }

protected:
    // default constructor
    CountedObject() { ++count(); }
    ~CountedObject() { --count(); }
    CountedObject(const this_type&) { ++count(); }
    explicit CountedObject(this_type&&) { ++count(); }

private:
    // number of existing objects
    static int64_t& thread_number();
    static int64_t& master_count();
    static int64_t& count();
    static int64_t  fmax_depth;
    static bool     fenabled;
};

//--------------------------------------------------------------------------------------//

template <typename CountedType>
int64_t&
CountedObject<CountedType>::thread_number()
{
    static std::atomic<int64_t> _all_instance;
    static thread_local int64_t _instance = _all_instance++;
    return _instance;
}

//--------------------------------------------------------------------------------------//

template <typename CountedType>
int64_t&
CountedObject<CountedType>::master_count()
{
    static int64_t _instance = 0;
    return _instance;
}

//--------------------------------------------------------------------------------------//

template <typename CountedType>
int64_t&
CountedObject<CountedType>::count()
{
    if(thread_number() == 0)
        return master_count();
    static thread_local int64_t _instance = master_count();
    return _instance;
}

//--------------------------------------------------------------------------------------//

template <typename CountedType>
int64_t CountedObject<CountedType>::fmax_depth = std::numeric_limits<int64_t>::max();

//--------------------------------------------------------------------------------------//

template <typename CountedType>
bool CountedObject<CountedType>::fenabled = true;

//======================================================================================//

}  // namespace PTL

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
class __void__;
