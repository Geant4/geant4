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
// $Id: HadronPhysicsFTFC.cc,v 1.4 2006/06/29 17:55:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

HadronPhysicsFTFC::HadronPhysicsFTFC(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsFTFC::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theFTFCNeutron=new G4FTFCNeutronBuilder);
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theFTFCPro=new G4FTFCProtonBuilder);
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMaxEnergy(25*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theFTFCPiK=new G4FTFCPiKBuilder);
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(25*GeV);
  
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

