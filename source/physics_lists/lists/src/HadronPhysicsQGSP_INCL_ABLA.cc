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
// $Id: HadronPhysicsQGSP_INCL_ABLA.cc,v 1.2 2010-06-03 10:42:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsQGSP_INCLABLA
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 23.11.2005 G.Folger: migration to non static particles
// 08.06.2006 V.Ivanchenko: remove stopping
// 20.06.2006 G.Folger: InclAbla applies to Kaons, i.e. use SetMinEnergy instead of SetMinPionEnergy
// 25.04.2007 G.Folger: Add code for quasielastic
// 10.12.2007 G.Folger: Add projectilediffrative option for proton/neutron, off by default
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsQGSP_INCL_ABLA.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSP_INCL_ABLA::HadronPhysicsQGSP_INCL_ABLA(G4int)
                    :  G4VPhysicsConstructor("hInelastic QGSP_INCL_ABLA")
		     , QuasiElastic(true)
{
   ProjectileDiffraction=false;
}

HadronPhysicsQGSP_INCL_ABLA::HadronPhysicsQGSP_INCL_ABLA(const G4String& name, G4bool quasiElastic)
                    :  G4VPhysicsConstructor(name) , QuasiElastic(quasiElastic)
{
   ProjectileDiffraction=false;
}

void HadronPhysicsQGSP_INCL_ABLA::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSPNeutron=new G4QGSPNeutronBuilder(QuasiElastic, ProjectileDiffraction));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(9.5*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);  

  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(2.9*GeV);
  theBertiniNeutron->SetMaxEnergy(9.9*GeV);

  theNeutrons->RegisterMe(theInclAblaNeutron=new G4InclAblaNeutronBuilder);
  theInclAblaNeutron->SetMinEnergy(0.0*GeV);
  theInclAblaNeutron->SetMaxEnergy(3.0*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSPPro=new G4QGSPProtonBuilder(QuasiElastic, ProjectileDiffraction));
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMinEnergy(9.5*GeV);
  theLEPPro->SetMaxEnergy(25*GeV);

  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMinEnergy(2.9*GeV);
  theBertiniPro->SetMaxEnergy(9.9*GeV);

  thePro->RegisterMe(theInclAblaPro=new G4InclAblaProtonBuilder);
  theInclAblaPro->SetMinEnergy(0.0*GeV);
  theInclAblaPro->SetMaxEnergy(3.0*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSPPiK=new G4QGSPPiKBuilder(QuasiElastic));
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(25*GeV);
  theLEPPiK->SetMinEnergy(9.5*GeV);

  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMinEnergy(2.9*GeV);
  theBertiniPiK->SetMaxEnergy(9.9*GeV);

  thePiK->RegisterMe(theInclAblaPiK=new G4InclAblaPiKBuilder);
  theInclAblaPiK->SetMinEnergy(0.0*GeV);
  theInclAblaPiK->SetMaxEnergy(3.0*GeV);
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsQGSP_INCL_ABLA::~HadronPhysicsQGSP_INCL_ABLA()
{
   delete theMiscLHEP;
   delete theQGSPNeutron;
   delete theLEPNeutron;
   delete theBertiniNeutron;
   delete theInclAblaNeutron;
   delete theQGSPPro;
   delete theLEPPro;
   delete thePro;
   delete theBertiniPro;
   delete theInclAblaPro;
   delete theQGSPPiK;
   delete theLEPPiK;
   delete theInclAblaPiK;
   delete theBertiniPiK;
   delete thePiK;
}

void HadronPhysicsQGSP_INCL_ABLA::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_INCL_ABLA::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
}

