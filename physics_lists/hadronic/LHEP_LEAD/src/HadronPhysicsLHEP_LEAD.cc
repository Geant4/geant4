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
// $Id: HadronPhysicsLHEP_LEAD.cc,v 1.5 2006/06/29 17:57:13 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsLHEP_LEAD
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsLHEP_LEAD.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_LEAD::HadronPhysicsLHEP_LEAD(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsLHEP_LEAD::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLEADNeutron=new G4LEADNeutronBuilder);
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);
  theLHEPNeutron->SetMinInelasticEnergy(4.99*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theLEADPro=new G4LEADProtonBuilder);
  thePro->RegisterMe(theLHEPPro=new G4LHEPProtonBuilder);
  theLHEPPro->SetMinEnergy(4.99*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLEADPiK=new G4LEADPiKBuilder);
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);
  theLHEPPiK->SetMinEnergy(4.99*GeV);

  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsLHEP_LEAD::~HadronPhysicsLHEP_LEAD()
{
  delete theNeutrons;
  delete theLHEPNeutron;
  delete theLEADNeutron;
  
  delete thePiK;
  delete theLHEPPiK;
  delete theLEADPiK;
  
  delete thePro;
  delete theLHEPPro;
  delete theLEADPro;    
  
  delete theMiscLHEP;
}

void HadronPhysicsLHEP_LEAD::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_LEAD::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

