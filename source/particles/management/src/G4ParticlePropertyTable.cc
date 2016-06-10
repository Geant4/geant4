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
// $Id: G4ParticlePropertyTable.cc 72955 2013-08-14 14:23:14Z gcosmo $
//
// class G4ParticlePropertyTable
//
// Implementation
//
// History:
// first implementation by H Kurashige 9 June 2003
// Add   magnetic moment    by H Kurashige   Mar 2007

#include "G4ios.hh"
#include "globals.hh"
#include "G4StateManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticlePropertyTable.hh"

// Static class variable: ptr to single instance of class
G4ThreadLocal G4ParticlePropertyTable* G4ParticlePropertyTable::fgParticlePropertyTable =0;

////////////////////
G4ParticlePropertyTable* G4ParticlePropertyTable::GetParticlePropertyTable()
{
  if (!fgParticlePropertyTable)
  {
    fgParticlePropertyTable = new  G4ParticlePropertyTable;
  }
  return fgParticlePropertyTable;
}

/////////////////////////////////////////////////////////////
G4ParticlePropertyTable::~G4ParticlePropertyTable()
{
  for (size_t idx=0; idx<arrayDataObject.size(); idx++){
    delete arrayDataObject[idx];
  }
  arrayDataObject.clear();
}

/////////////////////////////////////////////////////////////
G4ParticlePropertyTable::G4ParticlePropertyTable():
  verboseLevel(1)
{
  fParticleTable = G4ParticleTable::GetParticleTable();   
}

////////////////////////
G4ParticlePropertyTable::G4ParticlePropertyTable(const G4ParticlePropertyTable &right)
{
  fParticleTable = G4ParticleTable::GetParticleTable();   
  *this = right;
}
      
////////////////////////
G4ParticlePropertyTable & G4ParticlePropertyTable::operator=(const G4ParticlePropertyTable &right)
{
  if (this != &right) {
    fParticleTable = right.fParticleTable;
    verboseLevel   = right.verboseLevel; 
  }
  return *this;
}
  
////////////////////////
G4int G4ParticlePropertyTable::operator==(const G4ParticlePropertyTable &) const
{
  return true;
}

////////////////////////
G4int G4ParticlePropertyTable::operator!=(const G4ParticlePropertyTable &) const
{
  return false;
}

/////////////////////////////////////////////////////////////
void G4ParticlePropertyTable::Clear()
{
  for (size_t idx=0; idx<arrayDataObject.size(); idx++){
    delete arrayDataObject[idx];
  }
  arrayDataObject.clear();
}

/////////////////////////////////////////////////////////////////////
G4ParticlePropertyData* G4ParticlePropertyTable::GetParticleProperty(const G4String& aParticleName)
{
  G4ParticleDefinition* aParticle = fParticleTable->FindParticle(aParticleName);
  if (aParticle ==0 ) return 0;

  return GetParticleProperty(aParticle);
}

//////////////////////////////////////
G4ParticlePropertyData*  G4ParticlePropertyTable::GetParticleProperty(const G4ParticleDefinition* aParticle)
{
  if (aParticle ==0 ) return 0;
  G4ParticlePropertyData* pData = new G4ParticlePropertyData(aParticle->GetParticleName());
  pData->thePDGMass        = aParticle->GetPDGMass();
  pData->thePDGWidth       = aParticle->GetPDGWidth();
  pData->thePDGCharge      = aParticle->GetPDGCharge();
  pData->thePDGiSpin       = aParticle->GetPDGiSpin();
  pData->thePDGiParity     = aParticle->GetPDGiParity();
  pData->thePDGiConjugation  = aParticle->GetPDGiConjugation();
  pData->thePDGiGParity    = aParticle->GetPDGiGParity();
  pData->thePDGiIsospin    = aParticle->GetPDGiIsospin();
  pData->thePDGiIsospin3   = aParticle->GetPDGiIsospin3();
  pData->thePDGMagneticMoment   = aParticle->GetPDGMagneticMoment();
  pData->theLeptonNumber   = aParticle->GetLeptonNumber();
  pData->theBaryonNumber   = aParticle->GetBaryonNumber();
  pData->thePDGEncoding    = aParticle->GetPDGEncoding();
  pData->theAntiPDGEncoding  = aParticle->GetAntiPDGEncoding();
  pData->thePDGLifeTime    = aParticle->GetPDGLifeTime(); 
  for (size_t flv=0; flv<G4ParticlePropertyData::NumberOfQuarkFlavor; ++flv) {
    pData->theQuarkContent[flv]     = aParticle->theQuarkContent[flv];
    pData->theAntiQuarkContent[flv] = aParticle->theAntiQuarkContent[flv];
  }

  arrayDataObject.push_back(pData);
 
  return pData;
}

////////////////////////// 
G4bool G4ParticlePropertyTable::SetParticleProperty(const G4ParticlePropertyData& pData)
{
  G4StateManager* pStateMan = G4StateManager::GetStateManager();
  if (pStateMan->GetCurrentState() != G4State_PreInit){
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << "G4ParticlePropertyTable::GetParticleProperty() ";
      G4cout << " for " << pData.theParticleName << G4endl;
      G4cout << " Particle properties can be modified only in Pre_Init state";
      G4cout << G4endl;
    }
#endif
    return false;
  } 

  G4ParticleDefinition* aParticle = fParticleTable->FindParticle(pData.theParticleName);
  if (aParticle ==0 ) {
#ifdef G4VERBOSE
    if (verboseLevel>1){
      G4cout << "G4ParticlePropertyTable::GetParticleProperty() ";
      G4cout << " for " << pData.theParticleName << G4endl;
      G4cout << " Particle does not exist" << G4endl;
    }
#endif
    return false;
  }

  if (pData.fPDGMassModified) { 
    aParticle->thePDGMass = pData.thePDGMass;
  }
  if (pData.fPDGWidthModified) {
    aParticle->thePDGMass = pData.thePDGMass;
  }
  if (pData.fPDGChargeModified) {
    aParticle->thePDGCharge  = pData.thePDGCharge;
  }
  if (pData.fPDGiSpinModified) {
    aParticle->thePDGiSpin = pData.thePDGiSpin;
    aParticle->thePDGSpin  = 0.5*pData.thePDGiSpin;
  }
  if (pData.fPDGiParityModified) {
    aParticle->thePDGiParity = pData.thePDGiParity;
  }
  if (pData.fPDGiConjugationModified) {
    aParticle->thePDGiConjugation = pData.thePDGiConjugation;
  }
  if (pData.fPDGiGParityModified) {
    aParticle->thePDGiGParity = pData.thePDGiGParity;
  }
  if (pData.fPDGiIsospinModified) {
    aParticle->thePDGiIsospin = pData.thePDGiIsospin;
    aParticle->thePDGIsospin  = 0.5*pData.thePDGiIsospin;
  }
  if (pData.fPDGiIsospin3Modified) {
    aParticle->thePDGiIsospin3 = pData.thePDGiIsospin3;
    aParticle->thePDGIsospin3  = 0.5*pData.thePDGiIsospin3;
  }
  if (pData.fPDGMagneticMomentModified) {
    aParticle->thePDGMagneticMoment = pData.thePDGMagneticMoment;
   }
  if (pData.fLeptonNumberModified) {
    aParticle->theLeptonNumber   = pData.theLeptonNumber;
  }
  if (pData.fBaryonNumberModified) {
    aParticle->theBaryonNumber   = pData.theBaryonNumber;
  }
  if (pData.fPDGEncodingModified) {
    aParticle->thePDGEncoding   =  pData.thePDGEncoding;
  }
  if (pData.fAntiPDGEncodingModified) {
    aParticle->theAntiPDGEncoding   =  pData.theAntiPDGEncoding;
  }
  if (pData.fPDGLifeTimeModified) {
    aParticle->thePDGLifeTime = pData.thePDGLifeTime; 
  }
  for (size_t flv=0; flv<G4ParticlePropertyData::NumberOfQuarkFlavor; ++flv) {
    if (pData.fQuarkContentModified){
      aParticle->theQuarkContent[flv] = pData.theQuarkContent[flv];
    }
    if (pData.fAntiQuarkContentModified){
      aParticle->theAntiQuarkContent[flv] = pData.theAntiQuarkContent[flv];
    }
  }
  
  return true;
}
