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
// $Id: HadronPhysicsQGSC_LEAD_HP.cc,v 1.1 2006/10/31 11:35:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsQGSC_LEAD_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSC_LEAD_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSC_LEAD_HP::HadronPhysicsQGSC_LEAD_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsQGSC_LEAD_HP::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSCNeutron=new G4QGSCNeutronBuilder);
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);
  theLEPNeutron->SetMinEnergy(19.9*MeV);
  theLEPNeutron->SetMinInelasticEnergy(4.99*GeV);
  theNeutrons->RegisterMe(theLEADNeutron=new G4LEADNeutronBuilder);
  theLEADNeutron->SetMinEnergy(19.9*MeV);
  theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSCPro=new G4QGSCProtonBuilder);
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMaxEnergy(25*GeV);
  theLEPPro->SetMinEnergy(4.99*GeV);
  thePro->RegisterMe(theLEADProton=new G4LEADProtonBuilder);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSCPiK=new G4QGSCPiKBuilder);
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(25*GeV);
  theLEPPiK->SetMinEnergy(4.99*GeV);
  thePiK->RegisterMe(theLEADPiK=new G4LEADPiKBuilder);

  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsQGSC_LEAD_HP::~HadronPhysicsQGSC_LEAD_HP()
{
  delete theNeutrons;
  delete theLEPNeutron;
  delete theQGSCNeutron;
  delete theLEADNeutron;
  delete theHPNeutron;
  
  delete thePiK;
  delete theLEPPiK;
  delete theQGSCPiK;
  delete theLEADPiK;
  
  delete thePro;
  delete theLEPPro;
  delete theQGSCPro;    
  delete theLEADProton;
  
  delete theMiscLHEP;
}

void HadronPhysicsQGSC_LEAD_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSC_LEAD_HP::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

