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
// $Id: G4HadronBuilderQGSP.cc,v 1.1 2005-11-11 22:57:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronBuilderQGSP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard 
//
//----------------------------------------------------------------------------
//

#include "G4HadronBuilderQGSP.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

G4HadronBuilderQGSP::G4HadronBuilderQGSP(const G4String& name)
  :  G4VPhysicsConstructor(name), wasActivated(false)
{}

G4HadronBuilderQGSP::~G4HadronBuilderQGSP() 
{
  if(wasActivated) {
    delete theNeutrons;
    delete theLEPNeutron;
    delete theQGSPNeutron;
    
    delete thePiK;
    delete theLEPPiK;
    delete theQGSPPiK;
    
    delete thePro;
    delete theLEPPro;
    delete theQGSPPro;    
    
    delete theMiscLHEP;
    delete theStoppingHadron;
  }
}

void G4HadronBuilderQGSP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

void G4HadronBuilderQGSP::ConstructProcess()
{
  wasActivated = true;

  // Instantiate cmponents
  theNeutrons = new G4NeutronBuilder();
  theLEPNeutron = new G4LEPNeutronBuilder();
  theQGSPNeutron = new G4QGSPNeutronBuilder();

  thePiK = new G4PiKBuilder();
  theLEPPiK = new G4LEPPiKBuilder();
  theQGSPPiK = new G4QGSPPiKBuilder();

  thePro = new G4ProtonBuilder();
  theLEPPro = new G4LEPProtonBuilder();
  theQGSPPro = new G4QGSPProtonBuilder(); 

  theMiscLHEP = new G4MiscLHEPBuilder();
  theStoppingHadron = new G4StoppingHadronBuilder();

  // Initialization
  theNeutrons->RegisterMe(theQGSPNeutron);
  theNeutrons->RegisterMe(theLEPNeutron);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);  

  thePro->RegisterMe(theQGSPPro);
  thePro->RegisterMe(theLEPPro);
  theLEPPro->SetMaxEnergy(25*GeV);
  
  thePiK->RegisterMe(theQGSPPiK);
  thePiK->RegisterMe(theLEPPiK);
  theLEPPiK->SetMaxEnergy(25*GeV);

  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
  theStoppingHadron->Build();
}

