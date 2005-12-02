//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HadronPhysicsLHEP_HP.cc,v 1.3 2005-12-02 17:16:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsLHEP_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
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
  theStoppingHadron=new G4StoppingHadronBuilder;
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
    delete theStoppingHadron;
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
  theStoppingHadron->Build();
}
// 2002 by J.P. Wellisch
