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
// $Id: G4BinaryPionBuilder.cc 62122 2012-10-01 09:33:38Z gcosmo $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4BinaryPionBuilder
//
// Author: 2011 Gunter Folger
//
//
//----------------------------------------------------------------------------
//
#include "G4BinaryPionBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"

G4BinaryPionBuilder::
G4BinaryPionBuilder()
{
  thePiData = new G4CrossSectionPairGG(new G4PiNuclearCrossSection(), 91*GeV);
  theMin = 0*GeV;
  theMax = 1.3*GeV;
  theModel = new G4BinaryCascade;
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax); 
}

G4BinaryPionBuilder::
~G4BinaryPionBuilder()
{
}

void G4BinaryPionBuilder::
Build(G4HadronElasticProcess * ) {}

void G4BinaryPionBuilder::
Build(G4PionPlusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(thePiData);
  aP->RegisterMe(theModel);
}

void G4BinaryPionBuilder::
Build(G4PionMinusInelasticProcess * aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP->AddDataSet(thePiData);
  aP->RegisterMe(theModel);
}
