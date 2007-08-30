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
// $Id: G4GPRTypeTraits.hh,v 1.2 2007-08-30 19:33:45 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. Creation. Loki style type traits:
//                       "Modern C++ Design, Andrei Alexandrescu"   
//
#ifndef G4GPRTYPETRAITS_HH
#define G4GPRTYPETRAITS_HH

#include "G4GPRTypeList.hh"

template <bool flag, typename T, typename U>
struct Select
{
  typedef T Result;
};

template <typename T, typename U>
struct Select<false, T, U>
{
  typedef U Result;
};

namespace Private {
  
  template <typename U> struct AddPointer
  {
    typedef U* Result;
  };
  
  template <typename U> struct AddPointer<U&>
  {
    typedef U* Result;
  };
  
  template <class U> struct AddReference
  {
    typedef U & Result;
  };
  
  template <class U> struct AddReference<U &>
  {
    typedef U & Result;
  };
  
  template <> struct AddReference<void>
  {
    typedef G4GPRNullType Result;
  };
  
  template <class U> struct AddParameterType
  {
    typedef const U & Result;
  };

  template <class U> struct AddParameterType<U &>
  {
    typedef U & Result;
  };
  
  template <> struct AddParameterType<void>
  {
    typedef G4GPRNullType Result;
  };
  
  template <typename T>
  struct IsFunctionPointerRaw
  {
    enum{result = 0};
  };
  
  template <typename T>
  struct IsFunctionPointerRaw<T(*)()> 
  {enum 
      {result = 1};
  };
  
  template <typename T, 
            typename P01>
  struct IsFunctionPointerRaw<T(*)(P01)> 
  {
    enum {result = 1};
  };

  template <typename T, 
	    typename P01, typename P02>
  struct IsFunctionPointerRaw<T(*)(P01, P02)> 
  {
  enum {result = 1};
  };

   template <typename T, 
	     typename P01, typename P02, typename P03>
  struct IsFunctionPointerRaw<T(*)(P01, P02, P03)> 
  {
    enum {result = 1};
  };

  template <typename T, 
	    typename P01, typename P02, typename P03, typename P04>
  struct IsFunctionPointerRaw<T(*)(P01, P02, P03, P04)> 
  {
    enum {result = 1};
  };

  template <typename T>
  struct IsMemberFunctionPointerRaw
  {
    enum{result = 0};
  };

  template <typename T, typename S>
  struct IsMemberFunctionPointerRaw<T (S::*)()> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01)> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01, typename P02>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02)> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01, typename P02, typename P03>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02, P03)> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
	    typename P01, typename P02, typename P03, typename P04>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02, P03, P04)> 
  {
    enum {result = 1};
  };

  // Const versions
  template <typename T, typename S>
  struct IsMemberFunctionPointerRaw<T (S::*)() const> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01) const> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01, typename P02>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02) const> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
            typename P01, typename P02, typename P03>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02, P03) const> 
  {
    enum {result = 1};
  };

  template <typename T, typename S, 
	    typename P01, typename P02, typename P03, typename P04>
  struct IsMemberFunctionPointerRaw<T (S::*)(P01, P02, P03, P04) const> 
  {
    enum {result = 1};
  };
               
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
class G4GPRTypeTraits {
  
private:
  
  template <class U> struct ReferenceTraits
  {
    enum { result = false };
    typedef U ReferredType;
  };
  
  template <class U> struct ReferenceTraits<U&>
  {
    enum { result = true };
    typedef U ReferredType;
  };
  
  template <class U> struct PointerTraits
  {
    enum { result = false };
    typedef G4GPRNullType PointeeType;
  };
  
  template <class U> struct PointerTraits<U*>
  {
    enum { result = true };
    typedef U PointeeType;
  };
  
  template <class U> struct PointerTraits<U*&>
  {
    enum { result = true };
    typedef U PointeeType;
  };
  
  template <class U> struct PToMTraits
  {
    enum { result = false };
  };
  
  template <class U, class V> struct PToMTraits<U V::*>
  {
    enum { result = true };
  };
  
  template <class U, class V> struct PToMTraits<U V::*&>
  {
    enum { result = true };
  };
  
  template <class U> struct FunctionPointerTraits
  {
    enum{ result = Private::IsFunctionPointerRaw<U>::result };
  };
  
  template <typename U> struct PToMFunctionTraits
  {
    enum{ result = Private::IsMemberFunctionPointerRaw<U>::result };
  };
  
  template <class U> struct UnConst
  {
    typedef U Result;
    enum { isConst = 0 };
  };
  
  template <class U> struct UnConst<const U>
  {
    typedef U Result;
    enum { isConst = 1 };
  };

  template <class U> struct UnConst<const U&>
  {
    typedef U& Result;
    enum { isConst = 1 };
  };
  
  template <class U> struct UnVolatile
  {
    typedef U Result;
    enum { isVolatile = 0 };
  };
  
public:
  typedef typename UnConst<T>::Result 
  NonConstType;
  typedef typename UnVolatile<T>::Result 
  NonVolatileType;
  typedef typename UnVolatile<typename UnConst<T>::Result>::Result 
  UnqualifiedType;
  typedef typename PointerTraits<UnqualifiedType>::PointeeType 
  PointeeType;
  typedef typename ReferenceTraits<T>::ReferredType 
  ReferredType;
  
  enum { isConst          = UnConst<T>::isConst };
  enum { isReference      = ReferenceTraits<UnqualifiedType>::result };
  enum { isFunction       = FunctionPointerTraits<typename Private::AddPointer<T>::Result >::result };
  enum { isFunctionPointer= FunctionPointerTraits<
	 typename ReferenceTraits<UnqualifiedType>::ReferredType >::result };
  enum { isMemberFunctionPointer= PToMFunctionTraits<
	 typename ReferenceTraits<UnqualifiedType>::ReferredType >::result };
  enum { isMemberPointer  = PToMTraits<
	 typename ReferenceTraits<UnqualifiedType>::ReferredType >::result ||
	 isMemberFunctionPointer };
  enum { isPointer        = PointerTraits<
	 typename ReferenceTraits<UnqualifiedType>::ReferredType >::result ||
	 isFunctionPointer };
  
  typedef typename Select<isPointer || isMemberPointer, T, typename Private::AddParameterType<T>::Result>::Result ParameterType;

};

// J. Tinslay - add DetermineFlavour
struct TypeTraitType {
  struct PtrToFunction{};
  struct PtrToMfn{};
  struct RegularPtr{};
  struct Object{};
};

template <typename Func>
struct DetermineFlavour {
  typedef TypeTraitType::Object Result;
};

template <typename T>
struct DetermineFlavour<T(*)()> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01>
struct DetermineFlavour<T(*)(P01)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01, typename P02>
struct DetermineFlavour<T(*)(P01, P02)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01, typename P02, typename P03>
struct DetermineFlavour<T(*)(P01, P02, P03)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01, typename P02, typename P03, typename P04>
struct DetermineFlavour<T(*)(P01, P02, P03, P04)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01, typename P02, typename P03>
struct DetermineFlavour<T(*const)(P01, P02, P03)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T, typename P01, typename P02, typename P03, typename P04>
struct DetermineFlavour<T(*const)(P01, P02, P03, P04)> 
{
  typedef TypeTraitType::PtrToFunction Result;
};

template <typename T>
struct DetermineFlavour<T*> 
{
    typedef TypeTraitType::RegularPtr Result;
};

template <typename T, typename S>
struct DetermineFlavour<T (S::*)()> 
{
  typedef TypeTraitType::PtrToMfn Result;
};

template <typename T, typename S, typename P01>
struct DetermineFlavour<T (S::*)(P01)> 
{
  typedef TypeTraitType::PtrToMfn Result;
};

template <typename T, typename S, 
	  typename P01, typename P02>
struct DetermineFlavour<T (S::*)(P01, P02)> 
{
  typedef TypeTraitType::PtrToMfn Result;
};

template <typename T, typename S, 
	  typename P01, typename P02, typename P03>
struct DetermineFlavour<T (S::*)(P01, P02, P03)> 
{
  typedef TypeTraitType::PtrToMfn Result;
};


template <typename T, typename S, 
	  typename P01, typename P02, typename P03, typename P04>
struct DetermineFlavour<T (S::*)(P01, P02, P03, P04)> 
{
  typedef TypeTraitType::PtrToMfn Result;
};

#endif 

