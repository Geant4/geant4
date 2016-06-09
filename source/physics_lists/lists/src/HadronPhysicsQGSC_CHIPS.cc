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
// ClassName:   HadronPhysicsQGSC_CHIPS
//
// Author: 2009  M. Kosov
//
// Modified:
//
//----------------------------------------------------------------------------
// Short description: In fact this is the definition of the Hadronic Inelastic
// physics. The definition of the Hadronic Elastic physics one can find in the
// G4HadronQElasticPhysics, which is stable (the same for all physics lists).
// The only "unstable" part of the physics is the Hadronic Inelastic physics,
// which is usually composed of the wixing of the High Energy Inelastic Model
// (HEIM) and the Low Energy Inelastic Model (LEIM), which are applied only for
// some hadrons (mostly nucleons and pi-mesons), above the LHEP model, which
// usually covers all particles (but for Sigma_0 ?) and sometimes covers the
// "hole" between the LEIM and HIME at intermediate energies. The name of the
// Physics list is usually have a form HEIM_LEIM and the inelastic interactions
// are defined in the HadronicPhysicsHEIM_LEIM class. So in this particular
// physics list the low energy model is CHIPS (G4QCollision process) and the
// high energy model is QGSC (QGS with the Energy Flow interface to CHIPS),
// which are in terms of the energy boundary are mixed not on the model level,
// but on the process level (G4DiscProcessMixer class). The LHEP is completely
// excluded from this physics list, because the MiscLHEP is substituted by the
// MiscQGSC class (QGS with the Energy Flow interface to CHIPS), covering all
// particles, which are not N, pi, or K, defined by the separate builders. 
//---------------------------------------------------------------------------
#include <iomanip>   

#include "HadronPhysicsQGSC_CHIPS.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSC_CHIPS::HadronPhysicsQGSC_CHIPS(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSC_CHIPS")
    , theNeut(0)
    , theQGSCNeut(0)
    , thePiK(0)
    , theQGSCPiK(0)
    , theProt(0)
    , theQGSCProt(0)
    , theMiscQGSC(0)
    , QuasiElastic(true)
{}

HadronPhysicsQGSC_CHIPS::HadronPhysicsQGSC_CHIPS(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , theNeut(0)
    , theQGSCNeut(0)
    , thePiK(0)
    , theQGSCPiK(0)
    , theProt(0)
    , theQGSCProt(0)
    , theMiscQGSC(0)
    , QuasiElastic(quasiElastic)
{}

void HadronPhysicsQGSC_CHIPS::CreateModels()
{
  theNeut = new G4QNeutronBuilder;
  theNeut->RegisterMe(theQGSCNeut=new G4QGSC_CHIPSNeutronBuilder(QuasiElastic));
  //theQGSCNeut = new G4QGSC_CHIPSNeutronBuilder(QuasiElastic));

  theQGSCNeut->SetMinEnergy(0.0*GeV);

  theProt = new G4QProtonBuilder;
  theProt->RegisterMe(theQGSCProt = new G4QGSC_CHIPSProtonBuilder(QuasiElastic));
  //theQGSCProt = new G4QGSC_CHIPSProtonBuilder(QuasiElastic);

  theQGSCProt->SetMinEnergy(0.0*GeV);

  thePiK = new G4PiKBuilder;
  thePiK->RegisterMe(theQGSCPiK=new G4QGSC_CHIPSPiKBuilder(QuasiElastic));
  
  theQGSCPiK->SetMinEnergy(0.0*GeV);
   
  theMiscQGSC=new G4MiscQGSCBuilder(0);           // No verbose (@@ to be developed)
}

HadronPhysicsQGSC_CHIPS::~HadronPhysicsQGSC_CHIPS() 
{
  delete theQGSCNeut;
  delete theNeut;

  delete theQGSCProt;
  delete theProt;

  delete theQGSCPiK;
  delete thePiK;

  delete theMiscQGSC;
}

void HadronPhysicsQGSC_CHIPS::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSC_CHIPS::ConstructProcess()
{
  CreateModels();
  theNeut->Build();
  theProt->Build();
  thePiK->Build();
  theMiscQGSC->Build();
}

