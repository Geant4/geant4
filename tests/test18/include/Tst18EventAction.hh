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
//  File:        Tst18EventAction.hh
//  Description: Event action for radioactive decay system test 
//  Author:      Dennis Wright (SLAC)
//                 (original by F. Lei DERA UK)
//  Date:        14 August 2013
//

#ifndef Tst18EventAction_h
#define Tst18EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Tst18EventActionMessenger;


class Tst18EventAction : public G4UserEventAction
{
  public:
    Tst18EventAction();
   ~Tst18EventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);
    void IncrementParticleNumber();
    
  private:
    Tst18EventActionMessenger* eventMessenger;
    G4int numberOfSecondariesPerEvent;
};

#endif

