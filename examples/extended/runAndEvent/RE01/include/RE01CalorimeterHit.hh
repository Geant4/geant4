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
// $Id: RE01CalorimeterHit.hh,v 1.3 2006-06-29 17:42:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G4AttDef;
class G4AttValue;

class RE01CalorimeterHit : public G4VHit
{
  public:

      RE01CalorimeterHit();
      RE01CalorimeterHit(G4LogicalVolume* logVol,G4int z,G4int phi);
      ~RE01CalorimeterHit();
      RE01CalorimeterHit(const RE01CalorimeterHit &right);
      const RE01CalorimeterHit& operator=(const RE01CalorimeterHit &right);
      G4int operator==(const RE01CalorimeterHit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      void Print();

  private:
      G4int ZCellID;
      G4int PhiCellID;
      G4double edep;
      G4ThreeVector pos;
      G4RotationMatrix rot;
      const G4LogicalVolume* pLogV;
      G4double edepByATrack;
      RE01TrackInformation trackInfo;

  public:
      inline void SetCellID(G4int z,G4int phi)
      {
        ZCellID = z;
        PhiCellID = phi;
      }
      inline G4int GetZ() { return ZCellID; }
      inline G4int GetPhi() { return PhiCellID; }
      inline void SetEdep(G4double de)
      { edep = de; edepByATrack = de; }
      inline void AddEdep(G4double de)
      { edep += de; edepByATrack += de; }
      inline G4double GetEdep()
      { return edep; }
      inline G4double GetEdepByATrack()
      { return edepByATrack; }
      inline void ClearEdepByATrack()
      { edepByATrack = 0.; trackInfo = RE01TrackInformation(); }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetRot(G4RotationMatrix rmat)
      { rot = rmat; }
      inline G4RotationMatrix GetRot()
      { return rot; }
      inline const G4LogicalVolume * GetLogV()
      { return pLogV; }
      inline void SetTrackInformation(const G4Track* aTrack)
      {
        RE01TrackInformation* anInfo = (RE01TrackInformation*)(aTrack->GetUserInformation());
        trackInfo = *anInfo;
      }
      inline RE01TrackInformation* GetTrackInformation()
      { return &trackInfo; }
};

typedef G4THitsCollection<RE01CalorimeterHit> RE01CalorimeterHitsCollection;

extern G4Allocator<RE01CalorimeterHit> RE01CalorimeterHitAllocator;

inline void* RE01CalorimeterHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) RE01CalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void RE01CalorimeterHit::operator delete(void *aHit)
{
  RE01CalorimeterHitAllocator.FreeSingle((RE01CalorimeterHit*) aHit);
}

#endif


