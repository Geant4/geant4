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
//
// Author: Mathieu Karamitros

/*
 * Compile time counter to use flags with switch and add new flags
 *
 *   - Usage 1
 *
 * class A {
 *  G4CT_COUNT_INIT(0)
 *  CT_COUNT(a)
 *  CT_COUNT(b)
 *  CT_COUNT(c)
 * };
 *
 * class B1 : public A{
 *  CT_COUNT(d)
 * };
 *
 * class B2 : public A{
 *  CT_COUNT(e)
 * };
 * class C : public B1{
 *  CT_COUNT(f)
 * };
 *
 * Result
 *  a = 0
 *  b = 1
 *  c = 2
 *  d & e = 3
 *  f = 4
 *
 *   - Usage 2
 *  namespace CounterNSpace{
 *    G4CT_COUNT_INIT(0)
 *    CT_COUNT(a)
 *    CT_COUNT(b)
 *    CT_COUNT(c)
 *  }
 *
 *  // extend
 *  namespace CounterNSpace{
 *    CT_COUNT(d)
 *    CT_COUNT(e)
 *    CT_COUNT(f)
 *  }
 */

#ifndef G4CTCounter_h
#define G4CTCounter_h

template<int N>
struct G4Number: public G4Number<N-1>{
  static constexpr int value=N;
};

template<>
struct G4Number<0>{
  static constexpr int value=0;
};

//------------------------------------------------------------------------------

#define G4CT_COUNT_INIT(init_value)     \
  static constexpr G4Number<init_value> \
   Counter(G4Number<init_value>) {      \
     return G4Number<init_value>();     \
   }

// Allow overridable maximum enum depth by setting G4CT_MAX_COUNT
#ifndef G4CT_MAX_COUNT
#ifndef __NVCOMPILER
// Maximum in Geant4 V11.0.0 is hardcoded to 255
#define G4CT_MAX_COUNT 255
#else
// NVC++ compiler default template instantiation recursion limit is 64.
// This can be changed by setting, e.g., "-Wc,--pending_instantiations=256"
#define G4CT_MAX_COUNT 63
#endif
#endif

#define G4CT_COUNT(flagName)                                                                  \
  static constexpr const int flagName = decltype(Counter(G4Number<G4CT_MAX_COUNT>{}))::value; \
  static constexpr G4Number<flagName + 1> Counter(G4Number<flagName + 1>)                     \
  {                                                                                           \
    static_assert(flagName + 1 < G4CT_MAX_COUNT, "Maximum enumeration count exeeded");        \
    return G4Number<flagName + 1>{};                                                          \
  }

//------------------------------------------------------------------------------
// On Win platforms, static functions must not be inlined, use the following
// macros to separate definition and implementation of static function
// to avoid potential crashes

#define G4CT_COUNT_INIT_DEF(init_value) \
  static constexpr G4Number<init_value> \
   Counter(G4Number<init_value>);

#define G4CT_COUNT_INIT_IMPL(enumName, init_value) \
  constexpr G4Number<init_value> \
   enumName::Counter(G4Number<init_value>){ \
     return G4Number<init_value>{}; \
   }

#define G4CT_COUNT_DEF(flagName) \
  static constexpr const int flagName = \
   decltype(Counter(G4Number<G4CT_MAX_COUNT>{}))::value; \
  static constexpr G4Number<flagName + 1> \
   Counter(G4Number<flagName + 1>);

#define G4CT_COUNT_IMPL(enumName, flagName) \
  constexpr G4Number<enumName::flagName + 1> \
   enumName::Counter(G4Number<enumName::flagName + 1>){ \
     return G4Number<enumName::flagName +1 >{}; \
   }

//------------------------------------------------------------------------------

#endif
