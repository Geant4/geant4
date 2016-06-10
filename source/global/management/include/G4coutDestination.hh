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
// $Id: G4coutDestination.hh 82279 2014-06-13 14:44:00Z gcosmo $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4coutDestination.hh
//
// ---------------------------------------------------------------
#ifndef G4COUTDESTINATION_HH
#define G4COUTDESTINATION_HH

#include "globals.hh"

class G4coutDestination
{
  public:

    G4coutDestination();
    virtual ~G4coutDestination();

    virtual G4int ReceiveG4cout(const G4String&);
    virtual G4int ReceiveG4cerr(const G4String&);
protected:
    //For MT: if master G4coutDestination derived
    //class wants to intercept the thread outputs
    //derived class should set this pointer.
    //Needed for some G4UIsession like GUIs
    static G4coutDestination* masterG4coutDestination;
};

#endif
