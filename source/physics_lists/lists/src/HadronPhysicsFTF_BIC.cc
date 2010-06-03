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
// $Id: HadronPhysicsFTF_BIC.cc,v 1.3 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsFTF_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsFTF_BIC::HadronPhysicsFTF_BIC(G4int)
                    :  G4VPhysicsConstructor("hInelastic FTF_BIC")
		     , QuasiElastic(false)
{}

HadronPhysicsFTF_BIC::HadronPhysicsFTF_BIC(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{}

void HadronPhysicsFTF_BIC::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;

  theNeutrons->RegisterMe(theFTFBinaryNeutron=new G4FTFBinaryNeutronBuilder(QuasiElastic));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.*eV);
  theLEPNeutron->SetMaxInelasticEnergy(0.*eV);  

  theNeutrons->RegisterMe(theBinaryNeutron=new G4BinaryNeutronBuilder);
  theBinaryNeutron->SetMaxEnergy(5.0*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theFTFBinaryPro=new G4FTFBinaryProtonBuilder(QuasiElastic));

  thePro->RegisterMe(theBinaryPro=new G4BinaryProtonBuilder);
  theBinaryPro->SetMaxEnergy(5.0*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theFTFBinaryPiK=new G4FTFBinaryPiKBuilder(QuasiElastic));
  thePiK->RegisterMe(theBICPiK = new G4BinaryPiKBuilder);
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMinPionEnergy(10*GeV);   // don't use LEP for pion
  theBICPiK->SetMaxEnergy(5*GeV);        //  use Binary up to 5GeV for pion
  theLEPPiK->SetMaxEnergy(5*GeV); 

  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsFTF_BIC::~HadronPhysicsFTF_BIC() 
{
   delete theMiscLHEP;
   delete theFTFBinaryNeutron;
   delete theLEPNeutron;
   delete theBinaryNeutron;
   delete theNeutrons;
   delete theFTFBinaryPro;
   delete theBinaryPro;
   delete thePro;
   delete theFTFBinaryPiK;
   delete theBICPiK;
   delete theLEPPiK;
   delete thePiK;
}

void HadronPhysicsFTF_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsFTF_BIC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

