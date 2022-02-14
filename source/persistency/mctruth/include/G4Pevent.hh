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
// G4Pevent
//
// Class Description:
//
// Geant4 event object for store.
//
// This class has pointers to MCTruth and G4Event.
//
// In the event store operation, this object will be created in the concrete
// class of G4VPEventIO::Store() method, and will be deleted immediately after
// creating persistent Geant4 event.
//
// In the event retrieve operation, this object will be created in the
// concrete Persistency::Retrieve() method.  The retrieved G4Pevent
// has to be deleted by the user method which called the Retrieve().

// Author: Youhei Morita, 18.11.2001
// --------------------------------------------------------------------
#ifndef G4PEVENT_HH
#define G4PEVENT_HH 1

#include "G4Event.hh"
#include "G4MCTEvent.hh"

class G4Pevent
{
  public:

    G4Pevent(G4MCTEvent* mctevt, G4Event* g4evt);
      // Constructor

    ~G4Pevent();
      // Destructor

    G4int GetEventID() { return m_id; }
      // Returns the event ID

    G4Event* GetEvent() { return f_g4evt; }
      // Returns the G4Event

    G4MCTEvent* GetMCTEvent() { return f_mctevt; }
      // Returns the MCTruth event

    G4int GetGenEventID() const { return genEventID; }
      // Returns the GenEvent ID

    void SetGenEventID(G4int id) { genEventID = id; }
      // Sets the GenEvent ID

  private:

    G4MCTEvent* f_mctevt = nullptr;
    G4Event* f_g4evt = nullptr;
    G4int genEventID = -1;
    G4int m_id = -1;
};

#endif
