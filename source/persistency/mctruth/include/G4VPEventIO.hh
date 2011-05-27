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
// File: G4VPEventIO.hh
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#ifndef V_P_EVENT_I_O_HH
#define V_P_EVENT_I_O_HH 1

#include "G4Event.hh"
#include "G4Pevent.hh"

// Class Description:
//   Abstract base class for storing and retrieving an event object

class G4VPEventIO
{
    public: // With description
      G4VPEventIO();
      // Constructor

      virtual ~G4VPEventIO() {};
      // Destructor

    public: // With description
      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      inline G4int CurrentEventID() { return m_currentEvtID; };
      // Returns the current event id.

      virtual G4bool Store( const G4Event* anEvent ) =0;
      // Store a Geant4 event.

      virtual G4bool Retrieve( G4Pevent*& anEvent ) =0;
      // Retrieve a Geant4 event.

      virtual G4bool Retrieve( G4Event*& anEvent ) =0;
      // Retrieve a Geant4 event.

    protected:
      G4int m_verbose;
      G4int m_currentEvtID;

}; // End of class G4VPEventIO

#endif

