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
// $Id: HadronPhysicsFTFC.cc,v 1.2 2007/06/01 15:20:06 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 28.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko:remove stopping
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsFTFC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsFTFC::HadronPhysicsFTFC(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{}

void HadronPhysicsFTFC::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theFTFCNeutron=new G4FTFCNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFCNeutron);
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMaxInelasticEnergy(5*GeV);

  thePro=new G4ProtonBuilder;
  theFTFCPro=new G4FTFCProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFCPro);
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMaxEnergy(5*GeV);
  
  thePiK=new G4PiKBuilder;
  theFTFCPiK=new G4FTFCPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFCPiK);
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(5*GeV);
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsFTFC::~HadronPhysicsFTFC() 
{
  delete theNeutrons;
  delete theLEPNeutron;
  delete theFTFCNeutron;
    
  delete thePiK;
  delete theLEPPiK;
  delete theFTFCPiK;
    
  delete thePro;
  delete theLEPPro;
  delete theFTFCPro;    
    
  delete theMiscLHEP;
}

void HadronPhysicsFTFC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsFTFC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

