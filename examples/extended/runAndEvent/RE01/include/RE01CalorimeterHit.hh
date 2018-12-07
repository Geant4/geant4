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
/// \file runAndEvent/RE01/include/RE01CalorimeterHit.hh
/// \brief Definition of the RE01CalorimeterHit class
//
//

#ifndef RE01CalorimeterHit_h
#define RE01CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "RE01TrackInformation.hh"

#include "G4Types.hh"

class G4AttDef;
class G4AttValue;

class RE01CalorimeterHit : public G4VHit
{
public:
  RE01CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi);
  virtual ~RE01CalorimeterHit();
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  virtual void Draw();
  virtual void Print();
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;

  inline G4int GetZ() { return fZCellID; }
  inline G4int GetPhi() { return fPhiCellID; }
  inline void SetEdep(G4double de)
  { fEdep = de; fEdepByATrack = de; }
  inline void AddEdep(G4double de)
  { fEdep += de; fEdepByATrack += de; }
  inline G4double GetEdep()
  { return fEdep; }
  inline G4double GetEdepByATrack()
  { return fEdepByATrack; }
  inline void ClearEdepByATrack()
  { fEdepByATrack = 0.; fTrackInfo = RE01TrackInformation(); }
  inline void SetPos(G4ThreeVector xyz)
  { fPos = xyz; }
  inline void SetRot(G4RotationMatrix rmat)
  { fRot = rmat; }
  inline void SetTrackInformation(const G4Track* aTrack)
  {
    RE01TrackInformation* anInfo = 
      (RE01TrackInformation*)(aTrack->GetUserInformation());
    fTrackInfo = *anInfo;
  }
  inline RE01TrackInformation* GetTrackInformation()
  { return &fTrackInfo; }

private:
  G4int fZCellID;
  G4int fPhiCellID;
  G4double fEdep;
  G4ThreeVector fPos;
  G4RotationMatrix fRot;
  const G4LogicalVolume* fPLogV;
  G4double fEdepByATrack;
  RE01TrackInformation fTrackInfo;

};

typedef G4THitsCollection<RE01CalorimeterHit> RE01CalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<RE01CalorimeterHit> * RE01CalorimeterHitAllocator;

inline void* RE01CalorimeterHit::operator new(size_t)
{
  if(!RE01CalorimeterHitAllocator)
    RE01CalorimeterHitAllocator = new G4Allocator<RE01CalorimeterHit>;
  return (void *) RE01CalorimeterHitAllocator->MallocSingle();
}

inline void RE01CalorimeterHit::operator delete(void *aHit)
{
  RE01CalorimeterHitAllocator->FreeSingle((RE01CalorimeterHit*) aHit);
}

#endif


