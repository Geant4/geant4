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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsLHEP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 21.11.2005 G.Folger: don't  keep processes as data members, but new these
// 16.10.2012 A.Ribon: renamed stopping builder
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsLHEP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsLHEP);

HadronPhysicsLHEP::HadronPhysicsLHEP(G4int)
    :  G4VPhysicsConstructor("hInelastic LHEP") 
    , theNeutrons(0)
    , theLHEPNeutron(0)
    , thePiK(0)
    , theLHEPPiK(0)
    , thePro(0)
    , theLHEPPro(0)
    , theMiscLHEP(0)
    , theStoppingHadron(0)
{}

HadronPhysicsLHEP::HadronPhysicsLHEP(const G4String& name)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLHEPNeutron(0)
    , thePiK(0)
    , theLHEPPiK(0)
    , thePro(0)
    , theLHEPPro(0)
    , theMiscLHEP(0)
    , theStoppingHadron(0)
{}

void HadronPhysicsLHEP::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theLHEPPro=new G4LHEPProtonBuilder);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);

  theMiscLHEP=new G4MiscLHEPBuilder;
  theStoppingHadron=new G4LHEPStoppingHadronBuilder;  
}

HadronPhysicsLHEP::~HadronPhysicsLHEP()
{
   delete theLHEPNeutron;
   delete theNeutrons;
   delete theLHEPPro;
   delete thePro;
   delete theLHEPPiK;
   delete thePiK;
}

void HadronPhysicsLHEP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
  theStoppingHadron->Build();
}
// 2002 by J.P. Wellisch
