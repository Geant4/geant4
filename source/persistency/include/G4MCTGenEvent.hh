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
//   G4MCTGenEvent.hh
//
// ====================================================================
#ifndef MCT_GEN_EVENT_H
#define MCT_GEN_EVENT_H

#include "G4Types.hh"
#include <iostream>
#include <vector>
#include "CLHEP/HepMC/GenEvent.h"
 
// ====================================================================
//
// class definition
//
// ====================================================================

class G4MCTGenEvent {
protected:
  std::vector<HepMC::GenEvent*> eventList;

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

  void Print(std::ostream& ostr= std::cout) const;
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
