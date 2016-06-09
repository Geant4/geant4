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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4LEPPiKBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 16.11.2005 G.Folger: don't  keep processes as data members, but new these
// 13.06.2006 G.Folger: (re)move elastic scatterring 
//
//----------------------------------------------------------------------------
//
#include "G4LEPPiKBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

G4LEPPiKBuilder::
G4LEPPiKBuilder() :
   theLEPionPlusModel(0),  theLEPionMinusModel(0), 
   theLEKaonPlusModel(0), theLEKaonMinusModel(0), 
   theLEKaonZeroLModel(0), theLEKaonZeroSModel(0)  
{
  theMin = 0;
  theMax = 55*GeV;
  theMinPion = theMin;
}

G4LEPPiKBuilder::
~G4LEPPiKBuilder() 
{
  delete theLEPionPlusModel;
}

void G4LEPPiKBuilder::
Build(G4HadronElasticProcess *)
{
    G4cout << "Info - G4LEPPiKBuilder::Build() not adding elastic" << G4endl;
}

void G4LEPPiKBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(theMinPion);
  theLEPionPlusModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEPionPlusModel);
}

void G4LEPPiKBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(theMinPion);
  theLEPionMinusModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEPionMinusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonPlusInelasticProcess * aP)
{
  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theLEKaonPlusModel->SetMinEnergy(theMin);
  theLEKaonPlusModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEKaonPlusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonMinusInelasticProcess * aP)
{
  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theLEKaonMinusModel->SetMaxEnergy(theMax);
  theLEKaonMinusModel->SetMinEnergy(theMin);
  aP->RegisterMe(theLEKaonMinusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonZeroLInelasticProcess * aP)
{
  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theLEKaonZeroLModel->SetMaxEnergy(theMax);
  theLEKaonZeroLModel->SetMinEnergy(theMin);
  aP->RegisterMe(theLEKaonZeroLModel);
}
 
void G4LEPPiKBuilder::
Build(G4KaonZeroSInelasticProcess * aP)
{
  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theLEKaonZeroSModel->SetMaxEnergy(theMax);
  theLEKaonZeroSModel->SetMinEnergy(theMin);
  aP->RegisterMe(theLEKaonZeroSModel);
}

// 2002 by J.P. Wellisch
