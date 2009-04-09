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
// $Id: HadronPhysicsQGSC_CHIPS.cc,v 1.3 2009-04-09 08:22:53 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_CHIPS
//
// Author: 2007  G.Folger
//           created from HadronPhysicsQGSC, created by J.P. Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSC_CHIPS.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSC_CHIPS::HadronPhysicsQGSC_CHIPS(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name)  , QuasiElastic(quasiElastic)
{}

void HadronPhysicsQGSC_CHIPS::CreateModels()
{
  theNeutrons = new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSCNeutron=new G4QGSC_CHIPSNeutronBuilder(QuasiElastic));
  //theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);

  // M.K. ???
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  

  theQGSCNeutron->SetMinEnergy(0.0*GeV);
  //theBertiniNeutron->SetMinEnergy(0.0*GeV);
  //theBertiniNeutron->SetMaxEnergy(9.0*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSCPro=new G4QGSC_CHIPSProtonBuilder(QuasiElastic));
  //thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);

  theQGSCPro->SetMinEnergy(0.0*GeV);
  //theBertiniPro->SetMaxEnergy(9.0*GeV);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSCPiK=new G4QGSC_CHIPSPiKBuilder(QuasiElastic));
  //thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  
  theQGSCPiK->SetMinEnergy(0.0*GeV);
  //theBertiniPiK->SetMaxEnergy(9.0*GeV);
   
  theMiscLHEP=new G4MiscLHEPBuilder;            // To be replaced by QGSC
  //theMiscQGSC=new G4MiscQGSCBuilder;
}

HadronPhysicsQGSC_CHIPS::~HadronPhysicsQGSC_CHIPS() 
{
   delete theQGSCNeutron;
   //delete theBertiniNeutron;
   delete theLEPNeutron;
   delete theNeutrons;

   delete theQGSCPro;
   //delete theBertiniPro;
   delete thePro;

   delete theQGSCPiK;
   //delete theBertiniPiK;
   delete thePiK;

   delete theMiscLHEP;
   //delete theMiscQGSC;
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
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
  //theMiscQGSC->Build();
}

