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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsCHIPS
//
// Author: 2009  M.Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsCHIPS.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsCHIPS::HadronPhysicsCHIPS(G4int)
:  G4VPhysicsConstructor( "CHIPS hadronic")
    , theInelasticCHIPS(0)
//  , verbosity(verbose)
{}

HadronPhysicsCHIPS::HadronPhysicsCHIPS(const G4String& name):  G4VPhysicsConstructor(name), theInelasticCHIPS(0)
{}

HadronPhysicsCHIPS::~HadronPhysicsCHIPS() 
{
   delete theInelasticCHIPS;
}

void HadronPhysicsCHIPS::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsCHIPS::ConstructProcess()
{
  theInelasticCHIPS = new G4QInelasticCHIPSBuilder(0); // No verbose (@@ to be developed)
  theInelasticCHIPS->Build();
}

