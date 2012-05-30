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
// ClassName:   HadronPhysicsCHIPS_HP
//
// Author: 2012  M.Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsCHIPS_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsCHIPS_HP::HadronPhysicsCHIPS_HP(G4int)
:  G4VPhysicsConstructor( "CHIPS_HP hadronic")
,theInelasticCHIPS_HP(0)
//  , verbosity(verbose)
{}

HadronPhysicsCHIPS_HP::HadronPhysicsCHIPS_HP(const G4String& name)
:  G4VPhysicsConstructor(name)
, theInelasticCHIPS_HP(0)
{}

HadronPhysicsCHIPS_HP::~HadronPhysicsCHIPS_HP() 
{
   delete theInelasticCHIPS_HP;
}

void HadronPhysicsCHIPS_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsCHIPS_HP::ConstructProcess()
{
  theInelasticCHIPS_HP = new G4QInelasticCHIPS_HPBuilder(0); // No verbose (@@)
  theInelasticCHIPS_HP->Build();
}

