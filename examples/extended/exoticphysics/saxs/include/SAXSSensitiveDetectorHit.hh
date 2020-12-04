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
/// \file SAXSSensitiveDetectorHit.hh
/// \brief Definition of the SAXSSensitiveDetectorHit class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSSensitiveDetectorHit_h
#define SAXSSensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Sensitive Detector hit
///
/// It records:
/// - the spatial coordinate of the hit
/// - the time at which the hit occurred
/// - the type, ID, weight, momentum and energy of the impinging particle 
/// - the ID of the impinging particle parent

class SAXSSensitiveDetectorHit : public G4VHit
{
public:
  SAXSSensitiveDetectorHit();    
  SAXSSensitiveDetectorHit(const SAXSSensitiveDetectorHit &right);
  virtual ~SAXSSensitiveDetectorHit(){;}

  const SAXSSensitiveDetectorHit& operator=(const SAXSSensitiveDetectorHit &right);
  
  int operator==(const SAXSSensitiveDetectorHit &right) const;
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  
  virtual void Draw();
  virtual void Print();
  
private:
  G4int fTrackID;
  G4int fTrackIDP;
  G4double fTime;
  G4ThreeVector fPos;
  G4ThreeVector fMom;
  G4double fEnergy;
  G4int fType;
  G4double fWeight;
  
public:
  inline void SetTrackID(G4int z) {fTrackID = z;}
  inline G4int GetTrackID() const {return fTrackID;}
  
  inline void SetTrackIDP(G4int z) {fTrackIDP = z;}
  inline G4int GetTrackIDP() const {return fTrackIDP;}
  
  inline void SetTime(G4double t) {fTime = t;}
  inline G4double GetTime() const {return fTime;}
  
  inline void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  inline G4ThreeVector GetPos() const {return fPos;}
  
  inline void SetMom(G4ThreeVector xyz) {fMom = xyz;}
  inline G4ThreeVector GetMom() const {return fMom;}
  
  inline void SetEnergy(G4double energy) {fEnergy = energy;}
  inline G4double GetEnergy() const {return fEnergy;}
  
  inline void SetType(G4int type) {fType = type;}
  inline G4int GetType() const {return fType;}
  
  inline void SetWeight(G4double weight) {fWeight = weight;}
  inline G4double GetWeight() const {return fWeight;}
};

using SensitiveDetectorHitsCollection = G4THitsCollection<SAXSSensitiveDetectorHit>;

extern G4ThreadLocal G4Allocator<SAXSSensitiveDetectorHit> *hitAllocator;

inline void* SAXSSensitiveDetectorHit::operator new(size_t)
{
  if (!hitAllocator)
  {
      hitAllocator = new G4Allocator<SAXSSensitiveDetectorHit>;
  }
  return hitAllocator->MallocSingle();
}

inline void SAXSSensitiveDetectorHit::operator delete(void *aHit)
{
    if (!hitAllocator)
    {
        hitAllocator = new G4Allocator<SAXSSensitiveDetectorHit>;
    }
    hitAllocator->FreeSingle((SAXSSensitiveDetectorHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

