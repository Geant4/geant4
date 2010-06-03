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
// $Id: HadronPhysicsQGS_BIC.cc,v 1.2 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGS_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGS_BIC::HadronPhysicsQGS_BIC(G4int)
                    :  G4VPhysicsConstructor("hInelastic QGS_BIC")
		     , QuasiElastic(true)
{}

HadronPhysicsQGS_BIC::HadronPhysicsQGS_BIC(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{}

void HadronPhysicsQGS_BIC::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;

  theNeutrons->RegisterMe(theQGSBinaryNeutron=new G4QGSBinaryNeutronBuilder(QuasiElastic));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(9.5*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);  

  theNeutrons->RegisterMe(theBinaryNeutron=new G4BinaryNeutronBuilder);
  theBinaryNeutron->SetMaxEnergy(9.9*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSBinaryPro=new G4QGSBinaryProtonBuilder(QuasiElastic));
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMinEnergy(9.5*GeV);
  theLEPPro->SetMaxEnergy(25*GeV);

  thePro->RegisterMe(theBinaryPro=new G4BinaryProtonBuilder);
  theBinaryPro->SetMaxEnergy(9.9*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSBinaryPiK=new G4QGSBinaryPiKBuilder(QuasiElastic));
  thePiK->RegisterMe(theBICPiK = new G4BinaryPiKBuilder);
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMinPionEnergy(1.2*GeV);
  theLEPPiK->SetMaxEnergy(25*GeV);

  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsQGS_BIC::~HadronPhysicsQGS_BIC() 
{
   delete theMiscLHEP;
   delete theQGSBinaryNeutron;
   delete theLEPNeutron;
   delete theBinaryNeutron;
   delete theQGSBinaryPro;
   delete theLEPPro;
   delete thePro;
   delete theBinaryPro;
   delete theQGSBinaryPiK;
   delete theBICPiK;
   delete theLEPPiK;
   delete thePiK;
}

void HadronPhysicsQGS_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGS_BIC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

