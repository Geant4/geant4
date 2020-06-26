// MIT License
//
// Copyright (c) 2020 Jonathan R. Madsen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#pragma once

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <set>
#include <thread>
#include <type_traits>

namespace PTL
{
/// \class PTL::Singleton
/// \brief Singleton object that allows a deleter class to be specified
///
template <typename Type,
          typename PointerT = std::unique_ptr<Type, std::default_delete<Type>>>
class Singleton
{
public:
    using this_type     = Singleton<Type, PointerT>;
    using thread_id_t   = std::thread::id;
    using mutex_t       = std::recursive_mutex;
    using auto_lock_t   = std::unique_lock<mutex_t>;
    using pointer       = Type*;
    using list_t        = std::set<pointer>;
    using smart_pointer = PointerT;
    using deleter_t     = std::function<void(PointerT&)>;

    template <bool B, typename T = int>
    using enable_if_t = typename std::enable_if<B, T>::type;

public:
    // Constructor and Destructors
    Singleton();
    Singleton(pointer);
    ~Singleton();

    Singleton(const Singleton&) = delete;
    Singleton(Singleton&&)      = delete;
    Singleton& operator=(const Singleton&) = delete;
    Singleton& operator=(Singleton&&) = delete;

public:
    // public static functions
    static pointer     GetInstance();
    static pointer     GetMasterInstance();
    static thread_id_t GetMasterThreadID() { return f_master_thread(); }
    static list_t      Children() { return f_children(); }
    static bool        IsMaster(pointer ptr) { return ptr == GetRawMasterInstance(); }
    static bool        IsMasterThread();
    static void        Insert(pointer);
    static void        Remove(pointer);
    static mutex_t&    GetMutex() { return f_mutex(); }

public:
    // public member function
    void Initialize();
    void Initialize(pointer);
    void Destroy();
    void Reset(pointer);
    void Reset();

    // since we are overloading delete we overload new
    void* operator new(size_t)
    {
        this_type* ptr = ::new this_type();
        return static_cast<void*>(ptr);
    }

    // overload delete so that f_master_instance is guaranteed to be
    // a nullptr after deletion
    void operator delete(void* ptr)
    {
        this_type* _instance = (this_type*) (ptr);
        ::delete _instance;
        if(std::this_thread::get_id() == f_master_thread())
            f_master_instance() = nullptr;
    }

protected:
    friend class Type;

    // instance functions that do not Initialize
    smart_pointer&        GetSmartInstance() { return _local_instance(); }
    static smart_pointer& GetSmartMasterInstance() { return _master_instance(); }

    // for checking but not allocating
    pointer GetRawInstance()
    {
        return IsMasterThread() ? f_master_instance() : _local_instance().get();
    }
    static pointer GetRawMasterInstance() { return f_master_instance(); }

private:
    // Private functions
    static smart_pointer& _local_instance()
    {
        static thread_local smart_pointer _instance = smart_pointer();
        return _instance;
    }

    static smart_pointer& _master_instance()
    {
        static smart_pointer _instance = smart_pointer();
        return _instance;
    }

    void* operator new[](std::size_t) noexcept { return nullptr; }
    void  operator delete[](void*) noexcept {}

    template <typename Tp = Type, typename PtrT = PointerT,
              enable_if_t<(std::is_same<PtrT, std::shared_ptr<Tp>>::value)> = 0>
    deleter_t& GetDeleter()
    {
        static deleter_t _instance = [](PointerT&) {};
        return _instance;
    }

    template <typename Tp = Type, typename PtrT = PointerT,
              enable_if_t<!(std::is_same<PtrT, std::shared_ptr<Tp>>::value)> = 0>
    deleter_t& GetDeleter()
    {
        static deleter_t _instance = [](PointerT& _master) {
            auto& del = _master.get_deleter();
            del(_master.get());
            _master.reset(nullptr);
        };
        return _instance;
    }

private:
    // Private variables
    struct persistent_data
    {
        thread_id_t m_master_thread = std::this_thread::get_id();
        mutex_t     m_mutex;
        pointer     m_master_instance = nullptr;
        list_t      m_children        = {};

        persistent_data()                       = default;
        ~persistent_data()                      = default;
        persistent_data(const persistent_data&) = delete;
        persistent_data(persistent_data&&)      = delete;
        persistent_data& operator=(const persistent_data&) = delete;
        persistent_data& operator=(persistent_data&&) = delete;

        persistent_data(pointer _master, std::thread::id _tid)
        : m_master_thread(_tid)
        , m_master_instance(_master)
        {}

        void reset()
        {
            m_master_instance = nullptr;
            m_children.clear();
        }
    };

    bool                m_IsMaster = false;
    static thread_id_t& f_master_thread();
    static mutex_t&     f_mutex();
    static pointer&     f_master_instance();
    static list_t&      f_children();

    static persistent_data& f_persistent_data()
    {
        static persistent_data _instance;
        return _instance;
    }
};

//======================================================================================//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::thread_id_t&
Singleton<Type, PointerT>::f_master_thread()
{
    return f_persistent_data().m_master_thread;
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::pointer&
Singleton<Type, PointerT>::f_master_instance()
{
    return f_persistent_data().m_master_instance;
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::mutex_t&
Singleton<Type, PointerT>::f_mutex()
{
    return f_persistent_data().m_mutex;
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::list_t&
Singleton<Type, PointerT>::f_children()
{
    return f_persistent_data().m_children;
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
Singleton<Type, PointerT>::Singleton()
{
    Initialize();
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
Singleton<Type, PointerT>::Singleton(pointer ptr)
{
    Initialize(ptr);
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
Singleton<Type, PointerT>::~Singleton()
{
    auto& del = GetDeleter();
    del(_master_instance());
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Initialize()
{
    if(!f_master_instance())
    {
        f_master_thread()   = std::this_thread::get_id();
        f_master_instance() = new Type();
    }
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Initialize(pointer ptr)
{
    if(!f_master_instance())
    {
        f_master_thread()   = std::this_thread::get_id();
        f_master_instance() = ptr;
    }
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Destroy()
{
    if(std::this_thread::get_id() == f_master_thread() && f_master_instance())
    {
        delete f_master_instance();
        f_master_instance() = nullptr;
    }
    else
    {
        remove(_local_instance().get());
    }
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::pointer
Singleton<Type, PointerT>::GetInstance()
{
    if(std::this_thread::get_id() == f_master_thread())
        return GetMasterInstance();
    else if(!_local_instance().get())
    {
        _local_instance().reset(new Type());
        Insert(_local_instance().get());
    }
    return _local_instance().get();
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
typename Singleton<Type, PointerT>::pointer
Singleton<Type, PointerT>::GetMasterInstance()
{
    if(!f_master_instance())
    {
        f_master_thread()   = std::this_thread::get_id();
        f_master_instance() = new Type();
    }
    return f_master_instance();
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Reset(pointer ptr)
{
    if(IsMaster(ptr))
    {
        if(_master_instance().get())
            _master_instance().reset();
        else if(f_master_instance())
        {
            auto& del = GetDeleter();
            del(_master_instance());
            f_master_instance() = nullptr;
        }
        f_persistent_data().reset();
    }
    else
    {
        _local_instance().reset();
    }
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Reset()
{
    if(IsMasterThread())
        _master_instance().reset();
    _local_instance().reset();
    f_persistent_data().reset();
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
bool
Singleton<Type, PointerT>::IsMasterThread()
{
    return std::this_thread::get_id() == f_master_thread();
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Insert(pointer itr)
{
    auto_lock_t l(f_mutex());
    f_children().insert(itr);
}

//--------------------------------------------------------------------------------------//

template <typename Type, typename PointerT>
void
Singleton<Type, PointerT>::Remove(pointer itr)
{
    auto_lock_t l(f_mutex());
    for(auto litr = f_children().begin(); litr != f_children().end(); ++litr)
    {
        if(*litr == itr)
        {
            f_children().erase(litr);
            break;
        }
    }
}

//--------------------------------------------------------------------------------------//

}  // namespace PTL
