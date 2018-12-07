//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file parallel/ThreadsafeScorers/include/G4atomic_defines.hh
/// \brief Definition of the G4atomic_defines class
//
//
//
//
/// This is a functional class for G4atomic. The functions in this
///     file are not intended to be used outside of their implementation
///     in G4atomic.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef G4atomic_defines_hh_
#define G4atomic_defines_hh_

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <functional>
#include <atomic>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace atomics
{
    //------------------------------------------------------------------------//
    namespace details
    {
        //--------------------------------------------------------------------//
        template <typename _Tp>
        using OpFunction = std::function<_Tp(const _Tp&, const _Tp&)>;
        //--------------------------------------------------------------------//
        template <typename _Tp>
        inline void do_fetch_and_store(std::atomic<_Tp>* _atomic,
                                       const _Tp& _value,
                                       std::memory_order mem_odr)
        {
            _atomic->store(_value, mem_odr);
        }
        //--------------------------------------------------------------------//
        template <typename _Tp>
        inline void do_fetch_and_store(std::atomic<_Tp>* _atomic,
                                       const std::atomic<_Tp>& _value,
                                       std::memory_order mem_odr)
        {
            _atomic->store(_value.load(), mem_odr);
        }
        //--------------------------------------------------------------------//
        template <typename _Tp>
        inline void do_compare_and_swap(std::atomic<_Tp>* _atomic,
                          const _Tp& _value,
                          const OpFunction<_Tp>& _operator,
                          std::memory_order mem_odr)
        {
            _Tp _expected = _Tp();
            do {
                _expected = _atomic->load();
            } while (!(_atomic->compare_exchange_weak(_expected,
                                                      _operator(_expected,
                                                                _value),
                                                      mem_odr)));
        }
        //--------------------------------------------------------------------//
        template <typename _Tp>
        inline void do_compare_and_swap(std::atomic<_Tp>* _atomic,
                                        const std::atomic<_Tp>&
                                        _atomic_value,
                                        const OpFunction<_Tp>& _operator,
                        std::memory_order mem_odr)
        {
            _Tp _expected = _Tp();
            do {
                _expected = _atomic->load();
            } while (!(_atomic->compare_exchange_weak(_expected,
                                                      _operator(_expected,
                                                          _atomic_value.load()),
                                                      mem_odr)));
        }
        //--------------------------------------------------------------------//
    }
    //------------------------------------------------------------------------//
    //  WITH ATOMIC TEMPLATE BASE TYPE AS SECOND PARAMETER
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(std::atomic<T>* _atomic,
                    const T& _desired,
                    std::memory_order mem_odr
                    = std::memory_order_seq_cst)
    {
        details::do_compare_and_swap(_atomic,
                                     _desired,
                                     details::OpFunction<T>
                                     ([](const T&, const T& y){return y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(std::atomic<T>& _atomic,
                    const T& _desired,
                    std::memory_order mem_odr
                    = std::memory_order_seq_cst)
    {
        set(&_atomic, _desired, mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void increment(std::atomic<T>* _atomic,
                          const T& _increment,
                          std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _increment,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x+y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void decrement(std::atomic<T>* _atomic, const T& _decrement,
                          std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _decrement,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x-y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void multiply(std::atomic<T>* _atomic, const T& _factor,
                         std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _factor,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x*y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void divide(std::atomic<T>* _atomic, const T& _factor,
                       std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _factor,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x/y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    //  WITH ATOMICS AS SECOND PARAMETER
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(std::atomic<T>* _atomic,
                    const std::atomic<T>& _atomic_desired,
                    std::memory_order mem_odr
                    = std::memory_order_seq_cst)
    {
        //details::do_fetch_and_store(_atomic, _desired);
        details::do_compare_and_swap(_atomic, _atomic_desired,
                                     details::OpFunction<T>
                                     ([](const T&, const T& y){return y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(std::atomic<T>& _atomic,
                    const std::atomic<T>& _atomic_desired,
                    std::memory_order mem_odr)
    {
        set(&_atomic, _atomic_desired, mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void increment(std::atomic<T>* _atomic,
                          const std::atomic<T>& _atomic_increment,
                          std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _atomic_increment,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x+y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void decrement(std::atomic<T>* _atomic,
                          const std::atomic<T>& _atomic_decrement,
                          std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _atomic_decrement,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x-y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void multiply(std::atomic<T>* _atomic,
                         const std::atomic<T>& _atomic_factor,
                         std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _atomic_factor,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x*y;}),
                                     mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void divide(std::atomic<T>* _atomic,
                       const std::atomic<T>& _atomic_factor,
                       std::memory_order mem_odr)
    {
        details::do_compare_and_swap(_atomic, _atomic_factor,
                                     details::OpFunction<T>
                                     ([](const T& x, const T& y){return x/y;}),
                                     mem_odr);
    }

    //------------------------------------------------------------------------//
    //               STANDARD OVERLOAD                                        //
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(T* _non_atomic, const T& _desired)
    {
      *_non_atomic = _desired;
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline void set(T& _non_atomic, const T& _desired)
    {
      set(&_non_atomic, _desired);
    }
    //------------------------------------------------------------------------//
    //               STL PAIR OVERLOAD                                        //
    //------------------------------------------------------------------------//
    //
    //------------------------------------------------------------------------//
    //  WITH ATOMIC TEMPLATE TYPE AS SECOND PARAMETER
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void set(std::pair<std::atomic<T>,
                              std::atomic<U> >* _atomic,
                    const std::pair<T, U>& _desired)
    {
        set(&_atomic->first, _desired.first);
        set(&_atomic->second, _desired.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void set(std::pair<std::atomic<T>,
                              std::atomic<U> >& _atomic,
                    const std::pair<T, U>& _desired)
    {
        set(&_atomic, _desired);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void increment(std::pair<std::atomic<T>,
                                    std::atomic<U> >* _atomic,
                          const std::pair<T, U>& _increment)
    {
        increment(&_atomic->first, _increment.first);
        increment(&_atomic->second, _increment.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void decrement(std::pair<std::atomic<T>,
                                    std::atomic<U> >* _atomic,
                          const std::pair<T, U>& _decrement)
    {
        decrement(&_atomic->first, _decrement.first);
        decrement(&_atomic->second, _decrement.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void multiply(std::pair<std::atomic<T>,
                                   std::atomic<U> >* _atomic,
                         const std::pair<T, U>& _factor)
    {
        multiply(&_atomic->first, _factor.first);
        multiply(&_atomic->second, _factor.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void divide(std::pair<std::atomic<T>,
                                 std::atomic<U> >* _atomic,
                       const std::pair<T, U>& _factor)
    {
        divide(&_atomic->first, _factor.first);
        divide(&_atomic->second, _factor.second);
    }
    //------------------------------------------------------------------------//
    //  WITH ATOMICS AS SECOND PARAMETER
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void set(std::pair<std::atomic<T>,
                              std::atomic<U> >* _atomic,
                    const std::pair<std::atomic<T>,
                                    std::atomic<U> >& _desired)
    {
        set(&_atomic->first, _desired.first);
        set(&_atomic->second, _desired.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void set(std::pair<std::atomic<T>,
                              std::atomic<U> >& _atomic,
                    const std::pair<std::atomic<T>,
                                    std::atomic<U> >& _desired)
    {
        set(&_atomic, _desired);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void increment(std::pair<std::atomic<T>,
                                    std::atomic<U> >* _atomic,
                          const std::pair<std::atomic<T>,
                                          std::atomic<U> >& _increment)
    {
        increment(&_atomic->first, _increment.first);
        increment(&_atomic->second, _increment.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void decrement(std::pair<std::atomic<T>,
                                    std::atomic<U> >* _atomic,
                          const std::pair<std::atomic<T>,
                                          std::atomic<U> >& _decrement)
    {
        decrement(&_atomic->first, _decrement.first);
        decrement(&_atomic->second, _decrement.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void multiply(std::pair<std::atomic<T>,
                                   std::atomic<U> >* _atomic,
                         const std::pair<std::atomic<T>,
                                         std::atomic<U> >& _factor)
    {
        multiply(&_atomic->first, _factor.first);
        multiply(&_atomic->second, _factor.second);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline void divide(std::pair<std::atomic<T>,
                                 std::atomic<U> >* _atomic,
                       const std::pair<std::atomic<T>,
                                       std::atomic<U> >& _factor)
    {
        divide(&_atomic->first, _factor.first);
        divide(&_atomic->second, _factor.second);
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    template <typename T>
    inline T get(const T& _non_atomic)
    {
        return _non_atomic;
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline T get(const T& _non_atomic, std::memory_order)
    {
        return _non_atomic;
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline T get(const std::atomic<T>& _atomic)
    {
        return _atomic.load();
    }
    //------------------------------------------------------------------------//
    template <typename T>
    inline T get(const std::atomic<T>& _atomic,
                 std::memory_order mem_odr)
    {
        return _atomic.load(mem_odr);
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline std::pair<T, U> get(const std::pair<std::atomic<T>,
                               std::atomic<U> >& _atomic)
    {
        return std::pair<T, U>(get(_atomic.first), get(_atomic.second));
    }
    //------------------------------------------------------------------------//
    template <typename T, typename U>
    inline std::pair<T, U> get(const std::pair<std::atomic<T>,
                               std::atomic<U> >& _atomic,
                               std::memory_order mem_odr)
    {
        return std::pair<T, U>(get(_atomic.first, mem_odr),
                               get(_atomic.second, mem_odr));
    }
    //------------------------------------------------------------------------//



    //------------------------------------------------------------------------//
    // for plain old data (POD) and pairs (e.g. std::pair<atomic<T>, atomic<U>>)
    template <typename _Tp_base, typename _Tp_atom>
    inline _Tp_base base(const _Tp_atom& _atomic)
    {
        return get(_atomic);
    }
    //------------------------------------------------------------------------//

} // namespace atomics

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // G4MULTITHREADED

#endif  // atomic_typedefs_hh_
