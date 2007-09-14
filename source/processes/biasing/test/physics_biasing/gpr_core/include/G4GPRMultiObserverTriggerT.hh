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
// $Id: G4GPRMultiObserverTriggerT.hh,v 1.6 2007-09-14 16:44:29 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRMULTIOBSERVERTRIGGERT_HH
#define G4GPRMULTIOBSERVERTRIGGERT_HH

#include "G4GPRObserverCollectionT.hh"

template <typename Scope>
class G4GPRMultiObserverTriggerT {

public:
    
  template <typename Ptr, typename PtrToMfn>
  G4GPRMultiObserverTriggerT(Ptr* ptr, PtrToMfn mfn, G4bool initState = true)
    :fWrapper("tst", ptr, mfn)
    ,fCached(initState)
  {}

  template <typename Input>
  G4GPRMultiObserverTriggerT(Input* input, G4bool initState = true)
    :fWrapper("tst", input)
    ,fCached(initState)
  {}

  ~G4GPRMultiObserverTriggerT() {}
  
  template <typename Pointer, typename PointerToMfn>
  void AddObserver(Pointer* pointer, PointerToMfn mfn) 
  {
    fObserverCollection.RegisterObserver("tmp", pointer, mfn);
  }
  
  void Fire() 
  {
    G4bool result = fWrapper();
    if (result == fCached) return;

    fCached = result;
    fObserverCollection();    
  }

  template <typename Arg1>
  void Fire(const Arg1& a1) 
  {
    G4bool result = fWrapper(a1);
    G4cout<<"jane here multitrigger::fire 1 arg,  #observers: "<<fObserverCollection.GetNumberOfObservers()<<", Result:"<<result<<", Cached value: "<<fCached<<G4endl;

    if (result == fCached) return;

    fCached = result;
    fObserverCollection();     
  }

  template <typename Arg1, typename Arg2>
  void Fire(const Arg1& a1, const Arg2& a2) 
  {
    G4bool result = fWrapper(a1, a2);
    //    G4cout<<"jane here multitrigger::fire 2 args size "<<fObserverCollection.GetNumberOfObservers()<<" "<<result<<" "<<fCached<<G4endl;
    if (result == fCached) return;

    fCached = result;
    fObserverCollection();     
  }

private:

  G4GPRObserverCollectionT<G4GPRNullType> fObserverCollection;
  typename Scope::TriggerWrapper fWrapper;
  G4bool fCached;
  
};

#endif
