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
// $Id: HadronPhysicsQGSC_QGSC.cc,v 1.4 2009/04/14 07:23:08 mkossov Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSC_QGSC
//
// Author: 2007  G.Folger
//           created from HadronPhysicsQGSC, created by J.P. Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSC_QGSC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSC_QGSC::HadronPhysicsQGSC_QGSC(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name)  , QuasiElastic(quasiElastic)
{}

void HadronPhysicsQGSC_QGSC::CreateModels()
{
  theNeutrons = new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSCNeutron=new G4QGSC_QGSCNeutronBuilder(QuasiElastic));
  //theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);

  // M.K. ???
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);  

  theQGSCNeutron->SetMinEnergy(0.0*GeV);
  //theBertiniNeutron->SetMinEnergy(0.0*GeV);
  //theBertiniNeutron->SetMaxEnergy(9.0*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSCPro=new G4QGSC_QGSCProtonBuilder(QuasiElastic));
  //thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);

  theQGSCPro->SetMinEnergy(0.0*GeV);
  //theBertiniPro->SetMaxEnergy(9.0*GeV);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSCPiK=new G4QGSC_QGSCPiKBuilder(QuasiElastic));
  //thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  
  theQGSCPiK->SetMinEnergy(0.0*GeV);
  //theBertiniPiK->SetMaxEnergy(9.0*GeV);
   
  //theMiscLHEP=new G4MiscLHEPBuilder;            // To be replaced by QGSC
  theMiscQGSC=new G4MiscQGSCBuilder(0);           // No verbose (@@ to be developed)
}

HadronPhysicsQGSC_QGSC::~HadronPhysicsQGSC_QGSC() 
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

   //delete theMiscLHEP;
   delete theMiscQGSC;
}

void HadronPhysicsQGSC_QGSC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSC_QGSC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  //theMiscLHEP->Build();
  theMiscQGSC->Build();
}

