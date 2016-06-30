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
// $Id: G4ElementData.cc 96794 2016-05-09 10:09:30Z gcosmo $
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections, shell cross sections
//              isotope cross sections
//
// Author:      V.Ivanchenko 10.03.2011
//
// Modifications:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ElementData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElementData::G4ElementData()
{
  name = "";
  for(G4int i=0; i<maxNumElements; ++i) {
    elmData[i] = nullptr;
    elm2Data[i] = nullptr;
    compLength[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElementData::~G4ElementData()
{
  for(G4int i=0; i<maxNumElements; ++i) {
    delete elmData[i];
    delete elm2Data[i];
    size_t n = compLength[i]; 
    //G4cout << "Z= " << i << " L= " << n << G4endl;
    for(size_t j=0; j<n; ++j) {
      //G4cout << "j= " << j << "  " << (compData[i])[j] << G4endl;
      delete (compData[i])[j];
    }
  }
}

void G4ElementData::InitialiseForElement(G4int Z, G4PhysicsVector* v)
{
  if(Z < 1 || Z >= maxNumElements) {
    G4cout << "G4ElementData::InitialiseForElement ERROR for " << name 
	   << "  Z = " << Z << " is out of range!" << G4endl;
    G4Exception("G4ElementData::InitialiseForElement()", "mat601", 
                 FatalException, "Wrong data handling");
    return;
  } 
  if(elmData[Z]) { delete elmData[Z]; }
  elmData[Z] = v;
}

void G4ElementData::InitialiseForElement(G4int Z, G4Physics2DVector* v)
{
  if(Z < 1 || Z >= maxNumElements) {
    G4cout << "G4ElementData::InitialiseForElement ERROR for " << name 
	   << "  Z = " << Z << " is out of range!" << G4endl;
    G4Exception("G4ElementData::InitialiseForElement()", "mat601", 
                 FatalException, "Wrong data handling");
    return;
  } 
  if(elm2Data[Z]) { delete elm2Data[Z]; }
  elm2Data[Z] = v;
}

void G4ElementData::InitialiseForComponent(G4int Z, G4int nComponents)
{
  if(Z < 1 || Z >= maxNumElements) {
    G4cout << "G4ElementData::InitialiseForComponent ERROR for " << name 
	   << "  Z = " << Z << " is out of range!" << G4endl;
    G4Exception("G4ElementData::InitialiseForComponent()", "mat602", 
                 FatalException, "Wrong data handling");	   
    return;
  }

  // reserve a new structure
  size_t n = compLength[Z];
  if(0 < n) { 
    for(size_t i=0; i<n; ++i) { delete (compData[Z])[i]; }
    (compData[Z]).clear();
    (compID[Z]).clear();
  }
  (compData[Z]).reserve(nComponents);
  (compID[Z]).reserve(nComponents);
}

void 
G4ElementData::AddComponent(G4int Z, G4int id, G4PhysicsVector* v)
{
  if(Z < 1 || Z >= maxNumElements) {
    G4cout << "G4ElementData::AddComponent ERROR for " << name 
	   << "  Z = " << Z << " is out of range!" << G4endl;
    G4Exception("G4ElementData::AddComponent()", "mat603", 
                 FatalException, "Wrong data handling");	   
    return;
  }
  (compData[Z]).push_back(v);
  (compID[Z]).push_back(id);
  ++compLength[Z];
}
