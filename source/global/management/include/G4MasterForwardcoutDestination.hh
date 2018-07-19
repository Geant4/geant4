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
//
// $Id: G4MasterForwardcoutDestination.hh 103582 2017-04-18 17:24:45Z adotti $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// Class description:
//
// This class has the simple role to forward the stream of messages
// to the master thread in a MT application.
// Master thread should implment a G4coutDestination to receive the messages.
// Messages are serialized via a Mutex.

//      ---------------- G4LockcoutDestination ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4MASTERFORWARDCOUTDESTINATION_HH_
#define G4MASTERFORWARDCOUTDESTINATION_HH_

#include <G4coutDestination.hh>

class G4MasterForwardcoutDestination : public G4coutDestination
{
  public:

    G4MasterForwardcoutDestination() = default;

    virtual ~G4MasterForwardcoutDestination();
    virtual G4int ReceiveG4cout(const G4String& msg) override;
    virtual G4int ReceiveG4cerr(const G4String& msg) override;
};

#endif /* G4MASTERFORWARDCOUTDESTINATION_HH_ */
