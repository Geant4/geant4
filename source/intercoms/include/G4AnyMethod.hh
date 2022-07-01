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
// G4AnyMethod
//
// Class description:
//
// The G4AnyMethod class represents any object method.
// The class only holds a member pointer.

// See http://www.boost.org/libs/any for Documentation.
// Copyright Kevlin Henney, 2000, 2001, 2002. All rights reserved.
//
// Permission to use, copy, modify, and distribute this software for any
// purpose is hereby granted without fee, provided that this copyright and
// permissions notice appear in all copies and derivatives.
//
// This software is provided "as is" without express or implied warranty.
// What:  variant At boost::any
// who:   contributed by Kevlin Henney,
//        with features contributed and bugs found by
//        Ed Brey, Mark Rodgers, Peter Dimov, and James Curran
// when:  July 2001
// --------------------------------------------------------------------
#ifndef G4AnyMethod_hh
#define G4AnyMethod_hh 1

#include <functional>

/** Bad Argument exception */
class G4BadArgument : public std::bad_cast
{
 public:
  G4BadArgument() = default;
  const char* what() const throw() override
  {
    return "G4BadArgument: failed operator()";
  }
};

#include <type_traits>
using std::remove_const;
using std::remove_reference;

class G4AnyMethod
{
  public:

    /** contructors */

   G4AnyMethod() = default;

   template <class S, class T>
   G4AnyMethod(S (T::*f)())
   {
     fContent = new FuncRef<S, T>(f);
    }

    template <class S, class T, class A0>
    G4AnyMethod(S (T::*f)(A0))
      : narg(1)
    {
      fContent = new FuncRef1<S, T, A0>(f);
    }

    template <class S, class T, class A0, class A1>
    G4AnyMethod(S (T::*f)(A0, A1))
      : narg(2)
    {
      fContent = new FuncRef2<S, T, A0, A1>(f);
    }

    G4AnyMethod(const G4AnyMethod& other)
      : fContent(other.fContent != nullptr ? other.fContent->Clone() : nullptr)
      , narg(other.narg)
    {}

    /** destructor */

    ~G4AnyMethod() { delete fContent; }

    G4AnyMethod& Swap(G4AnyMethod& rhs)
    {
      std::swap(fContent, rhs.fContent);
      std::swap(narg, rhs.narg);
      return *this;
    }

    /** Assignment operators */

    template <class S, class T>
    G4AnyMethod& operator=(S (T::*f)())
    {
      G4AnyMethod(f).Swap(*this);
      narg = 0;
      return *this;
    }

    template <class S, class T, class A0>
    G4AnyMethod& operator=(S (T::*f)(A0))
    {
      G4AnyMethod(f).Swap(*this);
      narg = 1;
      return *this;
    }
    template <class S, class T, class A0, class A1>
    G4AnyMethod& operator=(S (T::*f)(A0, A1))
    {
      G4AnyMethod(f).Swap(*this);
      narg = 1;
      return *this;
    }

    G4AnyMethod& operator=(const G4AnyMethod& rhs)
    {
      G4AnyMethod(rhs).Swap(*this);
      narg = rhs.narg;
      return *this;
    }

    /** Query */

    G4bool Empty() const { return fContent == nullptr; }

    /** call operators */

    void operator()(void* obj) { fContent->operator()(obj); }
    void operator()(void* obj, const std::string& a0)
    {
      fContent->operator()(obj, a0);
    }

    /** Number of arguments */

    std::size_t NArg() const { return narg; }

    const std::type_info& ArgType(size_t n = 0) const
    {
      return fContent != nullptr ? fContent->ArgType(n) : typeid(void);
    }

  private:

    class Placeholder
    {
      public:
       Placeholder()                                       = default;
       virtual ~Placeholder()                              = default;
       virtual Placeholder* Clone() const                  = 0;
       virtual void operator()(void*)                      = 0;
       virtual void operator()(void*, const std::string&)  = 0;
       virtual const std::type_info& ArgType(size_t) const = 0;
    };

    template <class S, class T>
    struct FuncRef : public Placeholder
    {
      FuncRef(S (T::*f)())
        : fRef(f)
      {}

      void operator()(void* obj) override { ((T*) obj->*fRef)(); }
      void operator()(void*, const std::string&) override
      {
        throw G4BadArgument();
      }
      Placeholder* Clone() const override { return new FuncRef(fRef); }
      const std::type_info& ArgType(std::size_t) const override
      {
        return typeid(void);
      }
      S (T::*fRef)();
    };

    template <class S, class T, class A0>
    struct FuncRef1 : public Placeholder
    {
      using nakedA0 =
        typename remove_const<typename remove_reference<A0>::type>::type;

      FuncRef1(S (T::*f)(A0))
        : fRef(f)
      {}

      void operator()(void*) override { throw G4BadArgument(); }
      void operator()(void* obj, const std::string& s0) override
      {
        nakedA0 a0;
        std::stringstream strs(s0);
        strs >> a0;
        ((T*) obj->*fRef)(a0);
      }
      Placeholder* Clone() const override { return new FuncRef1(fRef); }
      const std::type_info& ArgType(size_t) const override
      {
        return typeid(A0);
      }
      S (T::*fRef)(A0);
    };

    template <class S, class T, class A0, class A1>
    struct FuncRef2 : public Placeholder
    {
      using nakedA0 =
        typename remove_const<typename remove_reference<A0>::type>::type;
      using nakedA1 =
        typename remove_const<typename remove_reference<A1>::type>::type;

      FuncRef2(S (T::*f)(A0, A1))
        : fRef(f)
      {}

      void operator()(void*) override { throw G4BadArgument(); }
      void operator()(void* obj, const std::string& s0) override
      {
        nakedA0 a0;
        nakedA1 a1;
        std::stringstream strs(s0);
        strs >> a0 >> a1;
        ((T*) obj->*fRef)(a0, a1);
      }
      Placeholder* Clone() const override { return new FuncRef2(fRef); }
      const std::type_info& ArgType(size_t i) const override
      {
        return i == 0 ? typeid(A0) : typeid(A1);
      }
      S (T::*fRef)(A0, A1);
    };

    Placeholder* fContent = nullptr;
    std::size_t narg = 0;
};

#endif
