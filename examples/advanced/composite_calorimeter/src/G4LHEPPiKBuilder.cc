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
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

G4LHEPPiKBuilder::
G4LHEPPiKBuilder()  
{
  theElasticModel = new G4LElastic();
  theM = 0;
}

G4LHEPPiKBuilder::
~G4LHEPPiKBuilder() {}

void G4LHEPPiKBuilder::
Build(G4HadronElasticProcess & aP)
{
  aP.RegisterMe(theElasticModel);
}

void G4LHEPPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theHEPionPlusModel = new G4HEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(theM);
  theLEPionPlusModel->SetMaxEnergy(55*GeV);
  theHEPionPlusModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLEPionPlusModel);
  aP.RegisterMe(theHEPionPlusModel);
}

void G4LHEPPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theHEPionMinusModel = new G4HEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(theM);
  theLEPionMinusModel->SetMaxEnergy(55*GeV);
  theHEPionMinusModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLEPionMinusModel);
  aP.RegisterMe(theHEPionMinusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theHEKaonPlusModel = new G4HEKaonPlusInelastic();
  theLEKaonPlusModel->SetMinEnergy(theM);
  aP.RegisterMe(theLEKaonPlusModel);
  aP.RegisterMe(theHEKaonPlusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theHEKaonMinusModel = new G4HEKaonMinusInelastic();
  theLEKaonMinusModel->SetMinEnergy(theM);
  aP.RegisterMe(theLEKaonMinusModel);
  aP.RegisterMe(theHEKaonMinusModel);
}

void G4LHEPPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  theLEKaonZeroLModel->SetMinEnergy(theM);
  aP.RegisterMe(theLEKaonZeroLModel);
  aP.RegisterMe(theHEKaonZeroLModel);
}
 
void G4LHEPPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  theLEKaonZeroSModel->SetMinEnergy(theM);
  aP.RegisterMe(theLEKaonZeroSModel);
  aP.RegisterMe(theHEKaonZeroSModel);
}

// 2002 by J.P. Wellisch
