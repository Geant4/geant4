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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
#include "G4PhantomHeadBuilder.hh"
#include "G4VBodyFactory.hh"
#include "G4MIRDBodyFactory.hh"
#include "G4ORNLFemaleBodyFactory.hh"
#include "G4ORNLMaleBodyFactory.hh"
#include "G4PVPlacement.hh"

G4PhantomHeadBuilder::G4PhantomHeadBuilder(): fModel("MIRD")
{  
  // sex can be "female" or "male"
  fBody = nullptr;
  fMotherVolume = nullptr;
  fHeadVolume = nullptr;
}

G4PhantomHeadBuilder::~G4PhantomHeadBuilder()
{
  delete fBody;
} 

void G4PhantomHeadBuilder::BuildHead(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fMotherVolume == nullptr)
    G4Exception("G4PhantomHeadBuilder::BuildHead()", "human_phantom0011", FatalException, "The moder Volume volume is missing !!!!!");
  else
  {
  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fHeadVolume = fBody -> CreateOrgan("Head", fMotherVolume, colourName, solidVis, sensitivity);
}
}

void G4PhantomHeadBuilder::BuildSkull(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
  if (fHeadVolume == nullptr)
    G4Exception("G4PhantomHeadBuilder::BuildSkull()", "human_phantom0012", FatalException, "The head volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fHeadVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  fBody -> CreateOrgan( "Skull", fHeadVolume, colourName, solidVis, sensitivity);
}

void G4PhantomHeadBuilder::BuildBrain(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{ 
 if (fHeadVolume == nullptr)
   G4Exception("G4PhantomHeadBuilder::BuildBrain()", "human_phantom0013", FatalException, "The head volume is missing !!!!!");
    else fBody -> CreateOrgan("Brain", fHeadVolume, colourName, solidVis, sensitivity);
}

G4VPhysicalVolume* G4PhantomHeadBuilder::GetPhantom()
{
  return fMotherVolume;
}

void G4PhantomHeadBuilder::SetMotherVolume(G4VPhysicalVolume* mother)
{
  fMotherVolume = mother;
}


void G4PhantomHeadBuilder::SetModel(G4String modelFlag)
{
  fModel = modelFlag;

  if(fModel=="MIRD") fBody = new G4MIRDBodyFactory();
  else if(fModel=="ORNLFemale") 
{
#ifdef G4LIB_USE_GDML
fBody = new G4ORNLFemaleBodyFactory();
#else
G4cout << fModel << " Working with GDML only! set G4LIB_USE_GDML 1" << G4endl;
#endif
} 
else if(fModel=="ORNLMale") 
{
#ifdef G4LIB_USE_GDML
fBody = new G4ORNLMaleBodyFactory();
#else
G4cout << fModel << " Working with GDML only! set G4LIB_USE_GDML 1" << G4endl;
#endif
}

}

