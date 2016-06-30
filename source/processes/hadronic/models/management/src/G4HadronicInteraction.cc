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
// $Id: G4HadronicInteraction.cc 96490 2016-04-19 06:57:04Z gcosmo $
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
#include "G4HadronicException.hh"

G4HadronicInteraction::G4HadronicInteraction(const G4String& modelName) :
  verboseLevel(0), theMinEnergy(0.0), theMaxEnergy(25.0*GeV), 
  isBlocked(false), recoilEnergyThreshold(0.0), theModelName(modelName),
  epCheckLevels(DBL_MAX, DBL_MAX)
{ 
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
  if( IsBlocked(aMaterial) ) { return 0.0; }
  if( IsBlocked(anElement) ) { return 0.0; }
  size_t length = theMinEnergyListElements.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
	if( anElement == theMinEnergyListElements[i].second )
	  { return theMinEnergyListElements[i].first; }
    }
  }
  length = theMinEnergyList.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( aMaterial == theMinEnergyList[i].second )
	{ return theMinEnergyList[i].first; }
    }
  }
  if(IsBlocked()) { return 0.0; }
  if( verboseLevel > 1 ) {
    G4cout << "*** Warning from HadronicInteraction::GetMinEnergy" << G4endl
           << "    material " << aMaterial->GetName()
           << " not found in min energy List" << G4endl;
  } 
  return theMinEnergy;
}
 
void G4HadronicInteraction::SetMinEnergy(G4double anEnergy,
					 const G4Element *anElement )
{
  if( IsBlocked(anElement) ) {
    G4cout << "*** Warning from HadronicInteraction::SetMinEnergy" << G4endl
           << "    The model is not active for the Element  "
           << anElement->GetName() << "." << G4endl;
  }
  size_t length = theMinEnergyListElements.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( anElement == theMinEnergyListElements[i].second )
	{
	  theMinEnergyListElements[i].first = anEnergy;
	  return;
	}
    }
  }
  theMinEnergyListElements.push_back(std::pair<G4double, const G4Element *>(anEnergy, anElement));
}
 
void G4HadronicInteraction::SetMinEnergy(G4double anEnergy,
					 const G4Material *aMaterial )
{
  if( IsBlocked(aMaterial) ) {
    G4cout << "*** Warning from HadronicInteraction::SetMinEnergy" << G4endl
           << "    The model is not active for the Material "
           << aMaterial->GetName() << "." << G4endl;
  }
  size_t length = theMinEnergyList.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( aMaterial == theMinEnergyList[i].second )
	{
	  theMinEnergyList[i].first = anEnergy;
	  return;
	}
    }
  }
  theMinEnergyList.push_back(std::pair<G4double, const G4Material *>(anEnergy, aMaterial));
}
 
G4double G4HadronicInteraction::GetMaxEnergy(const G4Material *aMaterial, 
					     const G4Element *anElement ) const
{
  if( IsBlocked(aMaterial) ) { return 0.0; }
  if( IsBlocked(anElement) ) { return 0.0; }
  size_t length = theMaxEnergyListElements.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
	if( anElement == theMaxEnergyListElements[i].second )
	  { return theMaxEnergyListElements[i].first; }
    }
  }
  length = theMaxEnergyList.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( aMaterial == theMaxEnergyList[i].second )
	{ return theMaxEnergyList[i].first; }
    }
  }
  if(IsBlocked()) { return 0.0; }
  if( verboseLevel > 1 ) {
    G4cout << "*** Warning from HadronicInteraction::GetMaxEnergy" << G4endl
           << "    material " << aMaterial->GetName()
           << " not found in min energy List" << G4endl;
  } 
  return theMaxEnergy;
}
 
void G4HadronicInteraction::SetMaxEnergy(G4double anEnergy,
					 const G4Element *anElement ) 
{
  if( IsBlocked(anElement) ) {
    G4cout << "*** Warning from HadronicInteraction::SetMaxEnergy" << G4endl
           << "Warning: The model is not active for the Element  "
           << anElement->GetName() << "." << G4endl;
  }
  size_t length = theMaxEnergyListElements.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( anElement == theMaxEnergyListElements[i].second )
      {
        theMaxEnergyListElements[i].first = anEnergy;
        return;
      }
    }
  }
  theMaxEnergyListElements.push_back(std::pair<G4double, const G4Element *>(anEnergy, anElement));
}

void G4HadronicInteraction::SetMaxEnergy(G4double anEnergy,
					 const G4Material *aMaterial )
{
  if( IsBlocked(aMaterial) ) {
    G4cout << "*** Warning from HadronicInteraction::SetMaxEnergy" << G4endl
           << "Warning: The model is not active for the Material "
           << aMaterial->GetName() << "." << G4endl;
  }
  size_t length = theMaxEnergyList.size();
  if(0 < length) {
    for(size_t i=0; i<length; ++i ) {
      if( aMaterial == theMaxEnergyList[i].second )
	{
	  theMaxEnergyList[i].first = anEnergy;
	  return;
	}
    }
  }
  theMaxEnergyList.push_back(std::pair<G4double, const G4Material *>(anEnergy, aMaterial));
}

void G4HadronicInteraction::DeActivateFor( const G4Material *aMaterial )
{
  theBlockedList.push_back(aMaterial);
}

void G4HadronicInteraction::DeActivateFor( const G4Element *anElement )
{
  theBlockedListElements.push_back(anElement);
}


G4bool G4HadronicInteraction::IsBlocked(const G4Material* aMaterial) const
{
  for (size_t i=0; i<theBlockedList.size(); ++i) {
    if (aMaterial == theBlockedList[i]) return true;
  }
  return false;
}


G4bool G4HadronicInteraction::IsBlocked(const G4Element* anElement) const
{
  for (size_t i=0; i<theBlockedListElements.size(); ++i) {
    if (anElement == theBlockedListElements[i]) return true;
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

/*
G4HadronicInteraction::G4HadronicInteraction(const G4HadronicInteraction &right )
{ 
  *this = right; 
}
    
const G4HadronicInteraction& 
G4HadronicInteraction::operator=(const G4HadronicInteraction &right )
{ 
  G4String text = "unintended use of G4HadronicInteraction::operator=";
  throw G4HadronicException(__FILE__, __LINE__, text); 
  return right;
}
 */
/* end of file */
 
