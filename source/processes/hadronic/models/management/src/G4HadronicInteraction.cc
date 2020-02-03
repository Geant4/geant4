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
// Hadronic Interaction  base class
// original by H.P. Wellisch
// modified by J.L. Chuma, TRIUMF, 21-Mar-1997
// Last modified: 04-Apr-1997
// reimplemented 1.11.2003 JPW.
// 23-Jan-2009 V.Ivanchenko move constructor and destructor to the body

#include <iostream>

#include "G4HadronicInteraction.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"

G4HadronicInteraction::G4HadronicInteraction(const G4String& modelName) :
  verboseLevel(0), theMinEnergy(0.0), 
  isBlocked(false), recoilEnergyThreshold(0.0), theModelName(modelName),
  epCheckLevels(DBL_MAX, DBL_MAX)
{ 
  theMaxEnergy = G4HadronicParameters::Instance()->GetMaxEnergy();
  registry = G4HadronicInteractionRegistry::Instance();
  registry->RegisterMe(this);
}

G4HadronicInteraction::~G4HadronicInteraction()
{
  registry->RemoveMe(this);
}

void G4HadronicInteraction::BuildPhysicsTable(const G4ParticleDefinition&)
{}

void G4HadronicInteraction::InitialiseModel()
{}

G4double 
G4HadronicInteraction::SampleInvariantT(const G4ParticleDefinition*, 
					G4double, G4int, G4int)
{
  return 0.0;
}

G4bool G4HadronicInteraction::IsApplicable(const G4HadProjectile&, 
					   G4Nucleus&)
{ 
  return true;
}

G4double G4HadronicInteraction::GetMinEnergy(
   const G4Material *aMaterial, const G4Element *anElement ) const
{
  if(!IsBlocked()) { return theMinEnergy; } 
  if( IsBlocked(aMaterial) || IsBlocked(anElement) ) { return DBL_MAX; }
  if(!theMinEnergyListElements.empty()) {
    for(auto const& elmlist : theMinEnergyListElements) {
	if( anElement == elmlist.second )
	  { return elmlist.first; }
    }
  }
  if(!theMinEnergyList.empty()) {
    for(auto const & matlist : theMinEnergyList) {
      if( aMaterial == matlist.second )
	{ return matlist.first; }
    }
  }
  return theMinEnergy;
}
 
void G4HadronicInteraction::SetMinEnergy(G4double anEnergy,
					 const G4Element *anElement )
{
  Block(); 
  if(!theMinEnergyListElements.empty()) {
    for(auto & elmlist : theMinEnergyListElements) {
      if( anElement == elmlist.second ) {
	elmlist.first = anEnergy;
	return;
      }
    }
  }
  theMinEnergyListElements.push_back(std::pair<G4double, const G4Element *>(anEnergy, anElement));
}
 
void G4HadronicInteraction::SetMinEnergy(G4double anEnergy,
					 const G4Material *aMaterial )
{
  Block(); 
  if(!theMinEnergyList.empty()) {
    for(auto & matlist : theMinEnergyList) {
      if( aMaterial == matlist.second ) {
	matlist.first = anEnergy;
	return;
      }
    }
  }
  theMinEnergyList.push_back(std::pair<G4double, const G4Material *>(anEnergy, aMaterial));
}
 
G4double G4HadronicInteraction::GetMaxEnergy(const G4Material *aMaterial, 
					     const G4Element *anElement ) const
{
  if(!IsBlocked()) { return theMaxEnergy; } 
  if( IsBlocked(aMaterial) || IsBlocked(anElement) ) { return 0.0; }
  if(!theMaxEnergyListElements.empty()) {
    for(auto const& elmlist : theMaxEnergyListElements) {
      if( anElement == elmlist.second )
	{ return elmlist.first; }
    }
  }
  if(!theMaxEnergyList.empty()) {
    for(auto const& matlist : theMaxEnergyList) {
      if( aMaterial == matlist.second )
	{ return matlist.first; }
    }
  }
  return theMaxEnergy;
}
 
void G4HadronicInteraction::SetMaxEnergy(G4double anEnergy,
					 const G4Element *anElement ) 
{
  Block(); 
  if(!theMaxEnergyListElements.empty()) {
    for(auto & elmlist : theMaxEnergyListElements) {
      if( anElement == elmlist.second ) {
        elmlist.first = anEnergy;
        return;
      }
    }
  }
  theMaxEnergyListElements.push_back(std::pair<G4double, const G4Element *>(anEnergy, anElement));
}

void G4HadronicInteraction::SetMaxEnergy(G4double anEnergy, const G4Material *aMaterial )
{
  Block(); 
  if(!theMaxEnergyList.empty()) {
    for(auto & matlist: theMaxEnergyList) {
      if( aMaterial == matlist.second ) {
	matlist.first = anEnergy;
	return;
      }
    }
  }
  theMaxEnergyList.push_back(std::pair<G4double, const G4Material *>(anEnergy, aMaterial));
}

void G4HadronicInteraction::DeActivateFor( const G4Material *aMaterial )
{
  Block(); 
  theBlockedList.push_back(aMaterial);
}

void G4HadronicInteraction::DeActivateFor( const G4Element *anElement )
{
  Block(); 
  theBlockedListElements.push_back(anElement);
}


G4bool G4HadronicInteraction::IsBlocked(const G4Material* aMaterial) const
{
  for (auto const& mat : theBlockedList) {
    if (aMaterial == mat) return true;
  }
  return false;
}


G4bool G4HadronicInteraction::IsBlocked(const G4Element* anElement) const
{
  for (auto const& elm : theBlockedListElements) {
    if (anElement == elm) return true;
  }
  return false;
}

const std::pair<G4double, G4double> G4HadronicInteraction::GetFatalEnergyCheckLevels() const
{
  // default level of Check
  return std::pair<G4double, G4double>(2.*perCent, 1. * GeV);
}

std::pair<G4double, G4double>
G4HadronicInteraction::GetEnergyMomentumCheckLevels() const
{
  return epCheckLevels;
}

void G4HadronicInteraction::ModelDescription(std::ostream& outFile) const
{
  outFile << "The description for this model has not been written yet.\n";
}

