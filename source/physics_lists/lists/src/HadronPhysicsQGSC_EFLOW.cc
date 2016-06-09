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
// $Id: HadronPhysicsQGSC_EFLOW.cc,v 1.2 2007/04/26 14:47:11 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_EFLOW
//
// Author: 2006 Gunter Folger
//
// Created from HadronPhysicsQGSC
// Modified:
// 25.04.2007 G.Folger: Add code for quasielastic
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSC_EFLOW.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSC_EFLOW::HadronPhysicsQGSC_EFLOW(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name)  , QuasiElastic(quasiElastic)
{}

void HadronPhysicsQGSC_EFLOW::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSCEflowNeutron=new G4QGSCEflowNeutronBuilder(QuasiElastic));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSCEflowPro=new G4QGSCEflowProtonBuilder(QuasiElastic));
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMaxEnergy(25*GeV);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSCEflowPiK=new G4QGSCEflowPiKBuilder(QuasiElastic));
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(25*GeV);
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsQGSC_EFLOW::~HadronPhysicsQGSC_EFLOW() 
{
   delete theMiscLHEP;
   delete theQGSCEflowNeutron;
   delete theLEPNeutron;
   delete theNeutrons;
   delete theQGSCEflowPro;
   delete theLEPPro;
   delete thePro;
   delete theQGSCEflowPiK;
   delete theLEPPiK;
   delete thePiK;
}

void HadronPhysicsQGSC_EFLOW::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSC_EFLOW::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

