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
// $Id: G4GPRTriggerManagerT.hh,v 1.1 2007-07-27 22:13:08 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRTRIGGERMANAGERT_HH
#define G4GPRTRIGGERMANAGERT_HH

#include "G4GPRLinearHierarchyT.hh"
#include "G4GPRTypeList.hh"
#include "G4GPRMultiObserverTriggerT.hh"
#include "G4GPRAssocT.hh"

template <typename Scope>
class G4GPRTriggerManagerT {

public:

  typedef G4GPRMultiObserverTriggerT<Scope> Trigger;
  typedef typename Scope::TriggerFunc TriggerFunc;
  typedef std::vector<Trigger*> TriggerList;
  typedef G4GPRAssocT<void*, Trigger*> PointerAssoc;
  typedef G4GPRAssocT<TriggerFunc, Trigger*> FuncAssoc;
  typedef G4GPRLinearHierarchyT< G4GPRTypeList_2(PointerAssoc, FuncAssoc) > Registry;

  template <typename ObserverPtr, typename PointerToMfn>
  void Register(const TriggerFunc& func, ObserverPtr* observer, PointerToMfn mfn) 
  {
    
    Trigger* trigger(0);

    if (!fRegistry.FuncAssoc::Retrieve(func, trigger)) {
      G4cout<<"jane didn't find func trigger"<<G4endl;
      trigger = new Trigger(func);
      fTriggerList.push_back(trigger);
      fRegistry.FuncAssoc::Register(func, trigger);
    }
    else {
      G4cout<<"jane did find func trigger"<<G4endl;
    }
    trigger->AddObserver(observer, mfn);
  }

  template <typename TriggerPtr, typename ObserverPtr, typename PointerToMfn>
  void Register(TriggerPtr*const & triggerPtr, ObserverPtr* observer, PointerToMfn mfn) 
  {
    Trigger* trigger(0);

    if (!fRegistry.PointerAssoc::Retrieve(triggerPtr, trigger)) {
      G4cout<<"jane didn't find ptr trigger"<<G4endl;
      typename Scope::template TriggerMfn<TriggerPtr>::PtrToMfn mfn = &TriggerPtr::operator();
      trigger = new Trigger(triggerPtr, mfn);
      fTriggerList.push_back(trigger);
      fRegistry.PointerAssoc::Register(triggerPtr, trigger);
    }
    else {
      G4cout<<"jane did find ptr trigger"<<G4endl;
    }
    trigger->AddObserver(observer, mfn);    
  }

  void Fire() 
  {
    for (typename TriggerList::iterator iter = fTriggerList.begin(); iter != fTriggerList.end(); ++iter) {
      (*iter)->Fire();
    }
  }

  template <typename Arg1>
  void Fire(const Arg1& a1) 
  {
    G4cout<<"jane triggermgr::fire 1 arg size "<<fTriggerList.size()<<G4endl;
    for (typename TriggerList::iterator iter = fTriggerList.begin(); iter != fTriggerList.end(); ++iter) {
      (*iter)->Fire(a1);
    }
  }

  template <typename Arg1, typename Arg2>
  void Fire(const Arg1& a1, const Arg2& a2) 
  {
    G4cout<<"jane triggermgr::fire 2 args size "<<fTriggerList.size()<<G4endl;
    for (typename TriggerList::iterator iter = fTriggerList.begin(); iter != fTriggerList.end(); ++iter) {
      (*iter)->Fire(a1, a2);
    }
  }

private:

  TriggerList fTriggerList;
  Registry fRegistry;

};



#endif
