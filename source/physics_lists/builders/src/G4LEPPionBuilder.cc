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
// ClassName:   G4LEPPionBuilder
//
// Author: 2010 G.Folger
//  devired from G4LEPPiKBuilder
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include "G4LEPPionBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

G4LEPPionBuilder::
G4LEPPionBuilder() : 
   theLEPionPlusModel(0),  theLEPionMinusModel(0) 
{
  theMin = 0;
  theMax = 55*GeV;
}

G4LEPPionBuilder::
~G4LEPPionBuilder() 
{
  delete theLEPionPlusModel;
}

void G4LEPPionBuilder::
Build(G4HadronElasticProcess *)
{
    G4cout << "Info - G4LEPPionBuilder::Build() not adding elastic" << G4endl;
}

void G4LEPPionBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(theMin);
  theLEPionPlusModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEPionPlusModel);
}

void G4LEPPionBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(theMin);
  theLEPionMinusModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEPionMinusModel);
}
