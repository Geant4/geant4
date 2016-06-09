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
// $Id: HadronPhysicsLHEP_BERT_HP.cc,v 1.1 2006/10/31 11:35:10 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:  HadronPhysicsLHEP_BERT_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsLHEP_BERT_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_BERT_HP::HadronPhysicsLHEP_BERT_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsLHEP_BERT_HP::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);
  theLHEPNeutron->SetMinEnergy(19.9*MeV);
  theLHEPNeutron->SetMinInelasticEnergy(9.5*GeV);
  theBertiniNeutron->SetMaxEnergy(9.9*GeV);
  theBertiniNeutron->SetMinEnergy(19.9*MeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theLHEPPro=new G4LHEPProtonBuilder);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theLHEPPro->SetMinEnergy(9.5*GeV);
  theBertiniPro->SetMaxEnergy(9.9*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsLHEP_BERT_HP::~HadronPhysicsLHEP_BERT_HP()
{
    delete theNeutrons;
    delete theLHEPNeutron;
    delete theBertiniNeutron;
    delete theHPNeutron;
    
    delete thePiK;
    delete theLHEPPiK;
    
    delete thePro;
    delete theLHEPPro;
    delete theBertiniPro;
    
    delete theMiscLHEP;
}

void HadronPhysicsLHEP_BERT_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_BERT_HP::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

