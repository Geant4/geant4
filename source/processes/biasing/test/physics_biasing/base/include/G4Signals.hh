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
// $Id: G4Signals.hh,v 1.1 2007-05-25 19:14:37 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation. Signal/slot idiom.   
//
#ifndef G4SIGNALS_HH
#define G4SIGNALS_HH

#include "G4Functor.hh"
#include "G4FunctorIdentifier.hh"
#include "G4TypeList.hh"
#include <vector>

template <typename TList>
class G4Signal {
  
  typedef G4Functor<void, G4FunctorIdentifier, TList> Slot;
  typedef G4FunctorIdentifier Identifier;
  
public:
  
  template <typename Id, typename PtrObj, typename MemFn>
  void Connect(const Id& id, const PtrObj&ptr, MemFn memFn) 
  {
    fSlots.push_back(Slot(id, ptr, memFn));
  }
  
  void AddSlot(const Slot& slot) 
  {
    fSlots.push_back(slot);
  }
  
  void operator()()
  {
    for (iter = fSlots.begin(); iter != fSlots.end(); ++iter) {
      (*iter)();
    }
  }
  
  template <typename Parm1>
  void operator()(const Parm1& p1)
  {
    for (iter = fSlots.begin(); iter != fSlots.end(); ++iter) {
      iter->operator()(p1);
    }
      
  }
  
  template <typename Parm1, typename Parm2>
  void operator()(const Parm1& p1, const Parm2& p2)
  {
    for (iter = fSlots.begin(); iter != fSlots.end(); ++iter) {
      (*iter)(p1, p2);
    }
  }
  
  unsigned GetNumberOfSlots() {return fSlots.size();}

private:

  typename std::vector<Slot>::iterator iter;
  std::vector<Slot> fSlots;

};

#endif
