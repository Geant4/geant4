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
// G4MCTGenEvent

// Author: Youhei Morita, 12.09.2001
// --------------------------------------------------------------------
#ifndef G4MCTGENEVENT_HH
#define G4MCTGENEVENT_HH 1

#include <iostream>
#include <vector>

#include "G4Types.hh"

class G4MCTGenEvent
{
  public:

    G4MCTGenEvent();
    virtual ~G4MCTGenEvent();

    inline G4MCTGenEvent(const G4MCTGenEvent& right);
    inline G4MCTGenEvent& operator=(const G4MCTGenEvent& right);
      // copy constructor and assignment operator

    G4int AddGenEvent(const void* genevent);
    G4int GetNofEvents() const;
    const void* GetGenEvent(G4int i);

    void ClearEvent();

  protected:

    std::vector<void*> eventList;
};

// ====================================================================
// inline methods
// ====================================================================

inline G4MCTGenEvent::G4MCTGenEvent(const G4MCTGenEvent& right)
{
  *this = right;
}

inline G4MCTGenEvent& G4MCTGenEvent::operator=(const G4MCTGenEvent& right)
{
  eventList = right.eventList;  // shallow copy

  return *this;
}

#endif
