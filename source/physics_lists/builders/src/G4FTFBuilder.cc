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
// $Id: G4FTFBuilder.cc,v 1.2 2009/11/25 19:33:18 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//---------------------------------------------------------------------------
//
// ClassName:  G4FTFBuilder
//
// Author: 28 June 2009 V.Ivanchenko
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "G4FTFBuilder.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4TheoFSGenerator.hh"
#include "G4QStringChipsParticleLevelInterface.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4LundStringFragmentation.hh"
#include "G4BinaryCascade.hh"


G4FTFBuilder::G4FTFBuilder(const G4String& aName) : G4VHadronModelBuilder(aName)
{}

G4FTFBuilder::~G4FTFBuilder() 
{
  delete theStringDecay;
  //delete theStringModel; // is deleted by the model
}                                     

G4HadronicInteraction* G4FTFBuilder::BuildModel()
{
  G4TheoFSGenerator* theFTFModel = new G4TheoFSGenerator(GetName());
  theStringModel  = new G4FTFModel();
  theStringDecay  = new G4ExcitedStringDecay(new G4LundStringFragmentation());
  theStringModel->SetFragmentationModel(theStringDecay);
  theFTFModel->SetHighEnergyGenerator(theStringModel);

  if(GetName() == "FTFP") {
    theFTFModel->SetTransport(new G4GeneratorPrecompoundInterface());
  } else if(GetName() == "FTFC") {
    theFTFModel->SetTransport(new G4QStringChipsParticleLevelInterface());
  } else if(GetName() == "FTFB") {
    theFTFModel->SetTransport(new G4BinaryCascade());
  } else {
    theFTFModel->SetTransport(new G4GeneratorPrecompoundInterface());
  }

  return theFTFModel;
}
