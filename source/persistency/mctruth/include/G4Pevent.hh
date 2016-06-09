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
// File: G4Pevent.hh
//
// History:
//   '01.11.18  Youhei Morita  Initial creation

#ifndef G4PEVENT_HH
#define G4PEVENT_HH 1

#include "G4Event.hh"
#include "G4MCTEvent.hh"

// Class Description:
//   Geant4 event object for store.
// 
//   This class has pointers to MCTruth and G4Event.
// 
//   In the event store operation, this object will be created in the concrete
//   class of G4VPEventIO::Store() method, and will be deleted immediately after
//   creating persistent Geant4 event. 
// 
//   In the event retrieve operation, this object will be created in the
//   concrete Persistency::Retrieve() method.  The retrieved G4Pevent
//   has to be deleted by the user method which called the Retrieve().

class G4Pevent
{
    public: // With description
      G4Pevent( G4MCTEvent* mctevt, G4Event* g4evt );
      // Constructor

      ~G4Pevent();
      // Destructor

    public: // With description
      int GetEventID() { return m_id; };
      // returns the event ID.

      G4Event* GetEvent() { return f_g4evt; };
      // returns the G4Event.

      G4MCTEvent* GetMCTEvent() { return f_mctevt; };
      // returns the MCTruth event.

      int GetGenEventID() const { return genEventID; };
      // returns the GenEvent ID.

      void SetGenEventID(int id) { genEventID=id; };
      // set the GenEvent ID.

    private:
      int       genEventID;
      G4MCTEvent* f_mctevt;
      G4Event*  f_g4evt;
      int       m_id;

}; // End of class G4Pevent

#endif

