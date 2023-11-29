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

#include "G4ElementData.hh"

G4ElementData::G4ElementData()
{
  for (G4int i = 0; i < maxNumElm; ++i) {
    elmData[i] = nullptr;
    elm2Data[i] = nullptr;
    compLength[i] = 0;
    compData[i] = nullptr;
    compID[i] = nullptr;
  }
}


G4ElementData::~G4ElementData()
{
  for (G4int i = 0; i < maxNumElm; ++i) {
    delete elmData[i];
    delete elm2Data[i];
    if (nullptr != compID[i]) {
      for (size_t j = 0; j < compID[i]->size(); ++j) {
        delete (*(compData[i]))[j];
      }
      delete compID[i];
      delete compData[i];
    }
  }
}

void G4ElementData::InitialiseForElement(G4int Z, G4PhysicsVector* v)
{
  if (Z < 1 || Z >= maxNumElm) {
    G4cout << "G4ElementData::InitialiseForElement ERROR for " << name << "  Z = " << Z
           << " is out of range!" << G4endl;
    G4Exception(
      "G4ElementData::InitialiseForElement()", "mat601", FatalException, "Wrong data handling");
    return;
  }
  if (nullptr != elmData[Z]) {
    delete elmData[Z];
  }
  elmData[Z] = v;
}

void G4ElementData::InitialiseForElement(G4int Z, G4Physics2DVector* v)
{
  if (Z < 1 || Z >= maxNumElm) {
    G4cout << "G4ElementData::InitialiseForElement ERROR for " << name << "  Z = " << Z
           << " is out of range!" << G4endl;
    G4Exception(
      "G4ElementData::InitialiseForElement()", "mat601", FatalException, "Wrong data handling");
    return;
  }
  if (nullptr != elm2Data[Z]) {
    delete elm2Data[Z];
  }
  elm2Data[Z] = v;
}

void G4ElementData::InitialiseForComponent(G4int Z, G4int nComponents)
{
  if (Z < 1 || Z >= maxNumElm || nComponents < 0) {
    G4cout << "G4ElementData::InitialiseForComponent ERROR for " << name << "  Z= " << Z
           << "  Ncomp= " << nComponents << " is out of range!" << G4endl;
    G4Exception(
      "G4ElementData::InitialiseForComponent()", "mat602", FatalException, "Wrong data handling");
    return;
  }

  // reserve a new structure
  if (nullptr == compID[Z]) {
    compID[Z] = new std::vector<G4int>();
    compData[Z] = new std::vector<G4PhysicsVector*>();
  }
  compID[Z]->resize(nComponents, -1);
  compData[Z]->resize(nComponents, nullptr);
  compLength[Z] = 0;
}

void G4ElementData::AddComponent(G4int Z, G4int id, G4PhysicsVector* v)
{
  if (Z < 1 || Z >= maxNumElm || (G4int)compID[Z]->size() == compLength[Z]) {
    G4cout << "G4ElementData::AddComponent ERROR for " << name << "  Z = " << Z
           << " is out of range!" << G4endl;
    G4Exception("G4ElementData::AddComponent()", "mat603", FatalException, "Wrong data handling");
    return;
  }
  (*(compData[Z]))[compLength[Z]] = v;
  (*(compID[Z]))[compLength[Z]] = id;
  compLength[Z] = compLength[Z] + 1;
}
