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
/// \file parallel/ThreadsafeScorers/include/G4atomic.hh
/// \brief Definition of the G4atomic class
//
//
//
//
/// This is an friendly implementation of the STL atomic class.
///     This class has the same interface as the STL atomic but can be used
///     in an extremely similar fashion to plain old data (POD) types.
///     In other words, the copy constructor and assignment operators are
///     defined, and a load() does not need to be called to get the POD
///     value.
///
/// IMPORTANT: Care must be taken when using this class as a RHS term.
///     The best use case scenario for this class is as a LHS term that is
///     only used as a RHS term outside of the multithreaded operations on it.
///
/// FOR EXAMPLE:
///     Proper use:
///         Goal: sum energy deposited in run
///         Impl: Is a member variable of derived
///                   G4VUserActionInitialization (of which there will
///                   always be just one instance). Accumulate
///                   thread-local energy deposit into EventAction,
///                   add to ActionInit at end of event, print sum
///                   on master EndOfRunAction
///         Why: the sum is only a LHS term
///     Improper use:
///         Goal: compute error during event processing
///         Impl: sum += x; sum_sq += x*x; counts++;
///               error = sqrt(sum_sq/(sum*sum) - 1/counts;
///               (where sum, sum_sq, counts are G4atomics)
///         Why: This will work but error can/will be wrong when
///              sum, sum_sq, and counts are updated by another thread
///              while error is being calculated, i.e. they are used as
///              RHS terms
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef G4atomic_hh_
#define G4atomic_hh_

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED

#include "G4atomic_defines.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename _Tp>
class G4atomic
{
  public:

    typedef typename std::atomic<_Tp>   base_type;
    typedef _Tp                         value_type;

  private:

    using mem_ord = std::memory_order;

  public:

    // constructors
    explicit
    G4atomic(mem_ord mo = std::memory_order_acq_rel) : fMemOrder(mo)
    { atomics::set(&fvalue, value_type()); }

    explicit
    G4atomic(const value_type& _init,
             mem_ord mo = std::memory_order_acq_rel) : fMemOrder(mo)
    { atomics::set(&fvalue, _init); }

    // copy-constructor from pure C++11 atomic
    explicit
    G4atomic(const base_type& rhs,
             mem_ord mo = std::memory_order_acq_rel) : fMemOrder(mo)
    { atomics::set(&fvalue, rhs); }

    // copy-constructor
    explicit
    G4atomic(const G4atomic& rhs) : fMemOrder(rhs.fMemOrder)
    { atomics::set(&fvalue, rhs.base()); }

    // assignment operators
    G4atomic& operator=(const G4atomic& rhs)
    {
        if(this != &rhs)
            atomics::set(&fvalue, rhs.fvalue);
        return *this;
    }

    G4atomic& operator=(const value_type& rhs)
    {
        atomics::set(&fvalue, rhs);
        return *this;
    }

    G4atomic& operator=(const base_type& rhs)
    {
        atomics::set(&fvalue, rhs);
        return *this;
    }

    // destructor
    ~G4atomic() { fvalue.~base_type(); }

    // base version
    base_type& base() { return fvalue; }
    const base_type& base() const { return fvalue; }
    base_type& base() volatile { return fvalue; }
    const base_type& base() const volatile { return fvalue; }

    // check if atomic is lock-free
    bool is_lock_free() const { return fvalue.is_lock_free(); }
    bool is_lock_free() const volatile { return fvalue.is_lock_free(); }

    // store functions
    void store(_Tp _desired, mem_ord mo = std::memory_order_seq_cst)
    { atomics::set(fvalue, _desired, mo); }
    void store(_Tp _desired, mem_ord mo = std::memory_order_seq_cst) volatile
    { atomics::set(fvalue, _desired, mo); }

    // load functions
    _Tp load(mem_ord mo = std::memory_order_seq_cst) const
    { return atomics::get(fvalue, mo); }
    _Tp load(mem_ord mo = std::memory_order_seq_cst) const volatile
    { return atomics::get(fvalue, mo); }

    // implicit conversion functions
    operator _Tp() const { return this->load(); }
    operator _Tp() const volatile { return this->load(); }

    operator base_type&() const { return fvalue; }

    // compare-and-swap functions
    bool compare_exchange_weak(_Tp& _expected, _Tp _desired,
                               mem_ord _success, mem_ord _failure)
    { return fvalue.compare_exchange_weak(_expected, _desired,
                                          _success, _failure); }
    bool compare_exchange_weak(_Tp& _expected, _Tp _desired,
                               mem_ord _success, mem_ord _failure) volatile
    { return fvalue.compare_exchange_weak(_expected, _desired,
                                          _success, _failure); }

    bool compare_exchange_weak(_Tp& _expected, _Tp _desired, mem_ord _order)
    { return fvalue.compare_exchange_weak(_expected, _desired, _order); }
    bool compare_exchange_weak(_Tp& _expected, _Tp _desired,
                               mem_ord _order) volatile
    { return fvalue.compare_exchange_weak(_expected, _desired, _order); }

    bool compare_exchange_strong(_Tp& _expected, _Tp _desired,
                                 mem_ord _success, mem_ord _failure)
    { return fvalue.compare_exchange_weak(_expected, _desired,
                                          _success, _failure); }
    bool compare_exchange_strong(_Tp& _expected, _Tp _desired,
                                 mem_ord _success, mem_ord _failure) volatile
    { return fvalue.compare_exchange_weak(_expected, _desired,
                                          _success, _failure); }

    bool compare_exchange_strong(_Tp& _expected, _Tp _desired, mem_ord _order)
    { return fvalue.compare_exchange_weak(_expected, _desired, _order); }
    bool compare_exchange_strong(_Tp& _expected, _Tp _desired,
                                 mem_ord _order) volatile
    { return fvalue.compare_exchange_weak(_expected, _desired, _order); }

    // value_type operators
    G4atomic& operator+=(const value_type& rhs)
    { atomics::increment(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator-=(const value_type& rhs)
    { atomics::decrement(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator*=(const value_type& rhs)
    { atomics::multiply(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator/=(const value_type& rhs)
    { atomics::divide(&fvalue, rhs, fMemOrder); return *this; }

    // atomic operators
    G4atomic& operator+=(const G4atomic& rhs)
    { atomics::increment(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator-=(const G4atomic& rhs)
    { atomics::decrement(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator*=(const G4atomic& rhs)
    { atomics::multiply(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator/=(const G4atomic& rhs)
    { atomics::divide(&fvalue, rhs.fvalue); return *this; }

    G4atomic& operator+=(const G4atomic& rhs) volatile
    { atomics::increment(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator-=(const G4atomic& rhs) volatile
    { atomics::decrement(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator*=(const G4atomic& rhs) volatile
    { atomics::multiply(&fvalue, rhs.fvalue); return *this; }
    G4atomic& operator/=(const G4atomic& rhs) volatile
    { atomics::divide(&fvalue, rhs.fvalue); return *this; }

    // STL atomic operators
    G4atomic& operator+=(const std::atomic<_Tp>& rhs)
    { atomics::increment(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator-=(const std::atomic<_Tp>& rhs)
    { atomics::decrement(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator*=(const std::atomic<_Tp>& rhs)
    { atomics::multiply(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator/=(const std::atomic<_Tp>& rhs)
    { atomics::divide(&fvalue, rhs, fMemOrder); return *this; }

    G4atomic& operator+=(const std::atomic<_Tp>& rhs) volatile
    { atomics::increment(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator-=(const std::atomic<_Tp>& rhs) volatile
    { atomics::decrement(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator*=(const std::atomic<_Tp>& rhs) volatile
    { atomics::multiply(&fvalue, rhs, fMemOrder); return *this; }
    G4atomic& operator/=(const std::atomic<_Tp>& rhs) volatile
    { atomics::divide(&fvalue, rhs, fMemOrder); return *this; }

    // increment operators
    value_type operator++() { value_type _tmp = ++fvalue; return _tmp; }
    value_type operator++(int)
    { value_type _tmp = fvalue++; return _tmp; }

    value_type operator--() { value_type _tmp = --fvalue; return _tmp; }
    value_type operator--(int)
    { value_type _tmp = fvalue--; return _tmp; }

  protected:

    base_type fvalue;
    mem_ord fMemOrder;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#else // ! G4MULTITHREADED

template <typename _Tp> using G4atomic = _Tp;

#endif // G4MULTITHREADED

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // G4atomic_hh_
