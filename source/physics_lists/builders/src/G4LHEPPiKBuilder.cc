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
#include "G4LHEPPiKBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

G4LHEPPiKBuilder::
G4LHEPPiKBuilder()  
{
  theM = 0;
  theMinPion = theM;
}

G4LHEPPiKBuilder::
~G4LHEPPiKBuilder() 
{
  delete theLEPionPlusModel;
  delete theHEPionPlusModel;
}

void G4LHEPPiKBuilder::
Build(G4HadronElasticProcess *)
{
}

void G4LHEPPiKBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theHEPionPlusModel = new G4HEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(theMinPion);
  theLEPionPlusModel->SetMaxEnergy(55*GeV);
  theHEPionPlusModel->SetMinEnergy(25*GeV);
  aP->RegisterMe(theLEPionPlusModel);
  aP->RegisterMe(theHEPionPlusModel);
}

void G4LHEPPiKBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theHEPionMinusModel = new G4HEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(theMinPion);
  theLEPionMinusModel->SetMaxEnergy(55*GeV);
  theHEPionMinusModel->SetMinEnergy(25*GeV);
  aP->RegisterMe(theLEPionMinusModel);
  aP->RegisterMe(theHEPionMinusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonPlusInelasticProcess * aP)
{
  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theHEKaonPlusModel = new G4HEKaonPlusInelastic();
  theLEKaonPlusModel->SetMinEnergy(theM);
  aP->RegisterMe(theLEKaonPlusModel);
  aP->RegisterMe(theHEKaonPlusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonMinusInelasticProcess * aP)
{
  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theHEKaonMinusModel = new G4HEKaonMinusInelastic();
  theLEKaonMinusModel->SetMinEnergy(theM);
  aP->RegisterMe(theLEKaonMinusModel);
  aP->RegisterMe(theHEKaonMinusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonZeroLInelasticProcess * aP)
{
  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  theLEKaonZeroLModel->SetMinEnergy(theM);
  aP->RegisterMe(theLEKaonZeroLModel);
  aP->RegisterMe(theHEKaonZeroLModel);
}
 
void G4LHEPPiKBuilder::
Build(G4KaonZeroSInelasticProcess * aP)
{
  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  theLEKaonZeroSModel->SetMinEnergy(theM);
  aP->RegisterMe(theLEKaonZeroSModel);
  aP->RegisterMe(theHEKaonZeroSModel);
}

// 2002 by J.P. Wellisch
