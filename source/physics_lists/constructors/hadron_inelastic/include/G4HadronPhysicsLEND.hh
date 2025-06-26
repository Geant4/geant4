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
//---------------------------------------------------------------------------
//
// Description: Create all hadronic (neutron and photon) processes using LEND which is valid up to 20 MeV
//
// Author: Douglas M Wright, LLNL 2022-04-25
//
//----------------------------------------------------------------------------
//

#ifndef G4HadronPhysicsLEND_h
#define G4HadronPhysicsLEND_h 1

#include "G4VPhysicsConstructor.hh"

const G4double maxLEND_Energy = 20*CLHEP::MeV;
const G4double overlapLEND_Energy = 0.1*CLHEP::MeV;

class G4HadronPhysicsLEND : public G4VPhysicsConstructor
{
public: 
    G4HadronPhysicsLEND(G4int verbose =1, const G4String& eval="");
    ~G4HadronPhysicsLEND() override = default;

    void ConstructProcess() override;
    void ConstructParticle() override;

    // copy constructor and hide assignment operator
    G4HadronPhysicsLEND(G4HadronPhysicsLEND &) = delete;
    G4HadronPhysicsLEND & operator =
    (const G4HadronPhysicsLEND &right) = delete;

private:
  G4String evaluation;
  G4int verbose;
};

#endif
