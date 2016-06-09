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
// $Id: HadronPhysicsLHEP_HP.cc,v 1.1 2006/10/31 11:35:10 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsLHEP_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsLHEP_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessManager.hh"

HadronPhysicsLHEP_HP::HadronPhysicsLHEP_HP(const G4String& name)  
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsLHEP_HP::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);
  theLHEPNeutron->SetMinEnergy(19.9*MeV);
  theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);

  theProton=new G4ProtonBuilder;
  theProton->RegisterMe(theLHEPProton=new G4LHEPProtonBuilder);

  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsLHEP_HP::
~HadronPhysicsLHEP_HP()
{
    delete theNeutrons;
    delete theLHEPNeutron;
    delete theHPNeutron;
    
    delete theLHEPPiK;
    delete thePiK;
    delete theLHEPProton;
    delete theProton;
    
    delete theMiscLHEP;
}

void HadronPhysicsLHEP_HP::
ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

void HadronPhysicsLHEP_HP::
ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  theProton->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

