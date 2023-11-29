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
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
//
//
#include"G4MaleBuilder.hh"
//#include "G4MIRDBodyFactory.hh"
//#include "G4ORNLMaleBodyFactory.hh"
#include "G4VBodyFactory.hh"

G4MaleBuilder::~G4MaleBuilder()
{
  delete fBody;
} 

void G4MaleBuilder::BuildMaleGenitalia(const G4String& colourName, G4bool solidVis, G4bool sensitivity)
{
if (fMotherVolume == nullptr)
    G4Exception("G4PhantomBuilder::BuildMaleGenitalia()", "human_phantom0048", FatalException, "The world volume is missing !!!!!");
  
  G4cout <<"MotherVolume: " << fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 

  fMaleGenitaliaVolume = fBody -> CreateOrgan("MaleGenitalia", fMotherVolume, colourName, solidVis, sensitivity); 
}

void G4MaleBuilder::BuildLeftTeste(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fMaleGenitaliaVolume == nullptr)
    G4Exception("G4FemaleBuilder::BuildLeftTeste()", "human_phantom0049", FatalException, "The maleGenitaliaVolume volume is missing !!!!!");

  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  
  fBody -> CreateOrgan("LeftTeste", fMaleGenitaliaVolume, colourName,solidVis, sensitivity); 
}

void G4MaleBuilder::BuildRightTeste(const G4String& colourName, G4bool solidVis, G4bool sensitivity )
{ 
  if (fMaleGenitaliaVolume == nullptr)
    G4Exception("G4FemaleBuilder::BuildRightTeste()", "human_phantom0050", FatalException, "The maleGenitaliaVolume volume is missing !!!!!");

  G4cout <<"MotherVolume: " <<  fMotherVolume -> GetName()<< G4endl;
  G4cout << "sensitivity : "<< sensitivity << G4endl; 
  
  fBody -> CreateOrgan("RightTeste", fMaleGenitaliaVolume, colourName,solidVis, sensitivity); 
}


