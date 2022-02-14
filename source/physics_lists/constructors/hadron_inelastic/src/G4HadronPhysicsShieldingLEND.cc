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
// ClassName:   
//
// Author: 7 Nov 2017 Tatsumi Koi
//   created from G4HadronPhysicsShielding
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4HadronPhysicsShieldingLEND.hh"
#include "G4HadronicParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsShieldingLEND);

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(G4int verb)
  :  G4HadronPhysicsShieldingLEND()
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
} 

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(const G4String& name)
  :  G4HadronPhysicsShieldingLEND(name, false)
{} 

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(
                              const G4String& name, G4bool qe)
  : G4HadronPhysicsShielding(name, qe)
{
  useLEND_ = true;
}

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(
                              const G4String& name, G4int verb) 
  :  G4HadronPhysicsShieldingLEND(name, false)
{
  G4HadronicParameters::Instance()->SetVerboseLevel(verb);
} 

G4HadronPhysicsShieldingLEND::G4HadronPhysicsShieldingLEND(
                              const G4String& name, G4int verb,
                              G4double minFTFPEnergy, G4double maxBertiniEnergy)
  :  G4HadronPhysicsShielding(name, verb, minFTFPEnergy, maxBertiniEnergy)
{
  useLEND_ = true;
}

G4HadronPhysicsShieldingLEND::~G4HadronPhysicsShieldingLEND()
{}

