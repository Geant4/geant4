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
#include "G4ElementDataRegistry.hh"
#include "G4PhysicsFreeVector.hh"

G4ElementData::G4ElementData(G4int length)
  : maxNumElm(length)
{
  elmData.resize((std::size_t)maxNumElm, nullptr);
  fRegistry = G4ElementDataRegistry::Instance();
  fRegistry->RegisterMe(this);
}

G4ElementData::~G4ElementData()
{
  for (auto const & p : elmData) { delete p; }
  elmData.clear();
 
  for (auto const & p : elm2Data) {
    delete p;
  }
  elm2Data.clear();

  for (auto const & p : compData) {
    if (nullptr != p ) {
      for (auto const & q : *p) { delete q.second; }
      delete p;
    }
  }
  compData.clear();

  for (auto const & p : comp2D) {
    if (nullptr != p ) {
      for (auto const & q : *p) { delete q.second; }
      delete p;
    }
  }
  comp2D.clear();
  fRegistry->RemoveMe(this);
}

void G4ElementData::Reserve1D(std::size_t n)
{
  if (!elmData.empty()) { return; }
  elmData.reserve(n);
  compData.reserve(n);
}

void G4ElementData::Reserve2D(std::size_t n)
{
  if (!elm2Data.empty()) { return; }
  elm2Data.reserve(n);
  comp2D.reserve(n);
}
  
void G4ElementData::InitialiseForElement(G4int Z, G4PhysicsVector* v)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "InitialiseForElement");
    return;
  }
  delete elmData[Z];
  elmData[Z] = v;
}

void G4ElementData::InitialiseForElement(G4int Z, G4Physics2DVector* v)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "InitialiseFor2DElement");
    return;
  }
  if (elm2Data.empty()) {
    elm2Data.resize((std::size_t)maxNumElm, nullptr);
  }
  delete elm2Data[Z];
  elm2Data[Z] = v;
}

void G4ElementData::InitialiseForComponent(G4int Z, G4int nComponents)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "InitialiseForComponent");
    return;
  }
  if (compData.empty()) {
    compData.resize((std::size_t)maxNumElm, nullptr);
  }
  delete compData[Z];
  compData[Z] = new std::vector<std::pair<G4int, G4PhysicsVector*> >;
  if (0 < nComponents) { compData[Z]->reserve((std::size_t)nComponents); }
}

void G4ElementData::InitialiseFor2DComponent(G4int Z, G4int nComponents)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "InitialiseFor2DComponent");
    return;
  }
  if (comp2D.empty()) {
    comp2D.resize((std::size_t)maxNumElm, nullptr);
  }
  delete comp2D[Z];
  comp2D[Z] = new std::vector<std::pair<G4int, G4Physics2DVector*> >;
  if (0 < nComponents) { comp2D[Z]->reserve((std::size_t)nComponents); }
}

void G4ElementData::AddComponent(G4int Z, G4int id, G4PhysicsVector* v)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "AddComponent");
    return;
  }
  if (compData.empty()) {
    compData.resize((std::size_t)maxNumElm, nullptr);
  }
  if (nullptr == compData[Z]) {
    compData[Z] = new std::vector<std::pair<G4int, G4PhysicsVector*> >;
  } 
  compData[Z]->emplace_back(id, v);
}

void G4ElementData::Add2DComponent(G4int Z, G4int id, G4Physics2DVector* v)
{
  if (Z < 0 || Z >= maxNumElm) {
    DataError(Z, "Add2DComponent");
    return;
  }
  if (comp2D.empty()) {
    compData.resize((std::size_t)maxNumElm, nullptr);
  }
  if (nullptr == comp2D[Z]) {
    comp2D[Z] = new std::vector<std::pair<G4int, G4Physics2DVector*> >;
  } 
  comp2D[Z]->emplace_back(id, v);
}

G4PhysicsVector* G4ElementData::New1DVector(G4int iz, G4int ne)
{
  auto ptr = new G4PhysicsFreeVector(ne);
  InitialiseForElement(iz, ptr);
  return ptr;
}

G4Physics2DVector* G4ElementData::New2DVector(G4int iz, G4int ny, G4int ne)
{
  auto ptr = new G4Physics2DVector(ny, ne);
  InitialiseForElement(iz, ptr);
  return ptr;
}

void G4ElementData::DataError(G4int Z, const G4String& type)
{
  G4cout << "G4ElementData::" << type << " ERROR for G4ElementData <"
         << name << ">  Z = " << Z << " is out of range!" << G4endl;
  G4Exception("G4ElementData", "mat603", FatalException, "Wrong data handling");
}
