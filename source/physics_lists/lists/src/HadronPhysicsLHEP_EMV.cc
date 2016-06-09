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
// ClassName:   HadronPhysicsLHEP_EMV
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.12.2005 G.Folger: create from HadronPhysicsLHEP_GN
// 08.06.2006 V.Ivanchenko: remove stopping
//
//----------------------------------------------------------------------------

#include "HadronPhysicsLHEP_EMV.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsLHEP_EMV);

HadronPhysicsLHEP_EMV::HadronPhysicsLHEP_EMV(G4int) 
    :  G4VPhysicsConstructor("hInelastic LHEP_EMV")
    , theNeutrons(0)
    , theLHEPNeutron(0)
    , thePiK(0)
    , theLHEPPiK(0)
    , thePro(0)
    , theLHEPPro(0)
    , theMiscLHEP(0)
{}

HadronPhysicsLHEP_EMV::HadronPhysicsLHEP_EMV(const G4String& name)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLHEPNeutron(0)
    , thePiK(0)
    , theLHEPPiK(0)
    , thePro(0)
    , theLHEPPro(0)
    , theMiscLHEP(0)
{}

void HadronPhysicsLHEP_EMV::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theLHEPPro=new G4LHEPProtonBuilder);

  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);
}

HadronPhysicsLHEP_EMV::~HadronPhysicsLHEP_EMV()
{
   delete theLHEPNeutron;
   delete theNeutrons;
   delete theLHEPPro;
   delete thePro;
   delete theLHEPPiK;
   delete thePiK;
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void HadronPhysicsLHEP_EMV::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_EMV::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

