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
// $Id: G4ObserverCollectionT.hh,v 1.1 2007-06-11 19:25:47 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Responsible for storing a list of observers and notifying the 
// observers when prompted. Observers need not be of the same type -
// they can be registered as a pointer to a function, or as a
// pointer to an object along with a pointer to a member function.
// The function or member function signature must be of the form:
//
// void MyFunc(Arguments)
//
// where "Arguments" are those specified in the G4TypeList which is used
// as the first template parameter.
// 
// Jane Tinslay, June 2007. Creation - basically rename G4Signal.  
//
#ifndef G4OBSERVERCOLLECTIONT_HH
#define G4OBSERVERCOLLECTIONT_HH

#include "G4Functor.hh"
#include "G4FunctorIdentifier.hh"
#include <vector>

template <typename TList, typename Identifier=G4FunctorIdentifier>
class G4ObserverCollectionT {
  
  // Wrapper definition - void return value since we're just notifying
  // observers that "something" happened which caused them to be notified.
  typedef G4Functor<void, Identifier, TList> Wrapper;
  
public:
  
  // Register pointer to member function observer
  template <typename Pointer>
  void RegisterObserver(const Identifier& id, const Pointer& pointer) 
  {
    fObserverCollection.push_back(Wrapper(id, pointer));
  }
 
  // Register pointer to object + pointer to member function observer
  template <typename Pointer, typename MemberFunc>
  void RegisterObserver(const Identifier& id, const Pointer& pointer, MemberFunc memberFunc) 
  {
    fObserverCollection.push_back(Wrapper(id, pointer, memberFunc));
  }
  
  // Notification operators. Multiple methods provided to allow the compiler
  // to select the operator corresponding to the number of arguments declared in the G4TypeList
  // (first template parameter). 
  void operator()()
  {
    for (iter = fObserverCollection.begin(); iter != fObserverCollection.end(); ++iter) {
      (*iter)();
    }
  }
  
  template <typename Arg1>
  void operator()(const Arg1& a1)
  {
    for (iter = fObserverCollection.begin(); iter != fObserverCollection.end(); ++iter) {
      iter->operator()(a1);
    }
      
  }
  
  template <typename Arg1, typename Arg2>
  void operator()(const Arg1& a1, const Arg2& a2)
  {
    for (iter = fObserverCollection.begin(); iter != fObserverCollection.end(); ++iter) {
      (*iter)(a1, a2);
    }
  }
  
  // Accessor
  unsigned GetNumberOfObservers() {return fObserverCollection.size();}

private:

  // Data members
  std::vector<Wrapper> fObserverCollection;
  typename std::vector<Wrapper>::iterator iter;

};

#endif

// jane fixme to do list:
// 1) Add accessor to observer collection
// 2) Add notify methods corresponding to operators to make things clearer ? Not sure.
