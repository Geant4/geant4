//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//   G4MCTGenEvent.hh
//
// ====================================================================
#ifndef MCT_GEN_EVENT_H
#define MCT_GEN_EVENT_H

#include "G4Types.hh"
#include "g4std/iostream"
#include "g4std/vector"
#include "CLHEP/HepMC/GenEvent.h"
 
// ====================================================================
//
// class definition
//
// ====================================================================

class G4MCTGenEvent {
protected:
  G4std::vector<HepMC::GenEvent*> eventList;

public:
  G4MCTGenEvent();
  virtual ~G4MCTGenEvent();
 
  // copy constructor and assignment operator
  G4MCTGenEvent(const G4MCTGenEvent& right);
  const G4MCTGenEvent& operator=(const G4MCTGenEvent& right);
  
  // methods...  
  int AddGenEvent(const HepMC::GenEvent* genevent);
  int GetNofEvents() const;
  const HepMC::GenEvent* GetGenEvent(int i);

  void ClearEvent();

  void Print(G4std::ostream& ostr= G4std::cout) const;
};

// ====================================================================
// inline functions
// ====================================================================

inline G4MCTGenEvent::G4MCTGenEvent(const G4MCTGenEvent& right)
{
  *this= right;
}
 
inline const G4MCTGenEvent& G4MCTGenEvent::operator=(const G4MCTGenEvent& right)
{
  eventList= right.eventList;  // shallow copy

  return *this;
}

#endif
