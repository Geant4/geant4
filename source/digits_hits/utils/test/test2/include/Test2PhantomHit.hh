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

#ifndef Test2PhantomHit_h
#define Test2PhantomHit_h 1

#include "globals.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ParticleDefinition.hh"


class Test2PhantomHit : public G4VHit {
public:
  Test2PhantomHit();
  Test2PhantomHit(G4int & x, G4int & y, G4int & z,G4int segment[3]);
  ~Test2PhantomHit();
  Test2PhantomHit(const Test2PhantomHit & right);
  const Test2PhantomHit & operator=(const Test2PhantomHit & right);
  G4int operator==(const Test2PhantomHit & right) const;

  inline void * operator new(size_t);
  inline void operator delete(void * aHit);

  void Print();

private:
  G4int fXCellID, fYCellID, fZCellID;
  G4int fSegment[3];
  G4double fEdep;
  G4double fDose;
  G4double fTrackLength;
  G4String fParticleName;

  G4double flatSurfCurr;
  G4double passageCellCurr;
  G4double flatSurfaceFlux;
  G4double cellFlux;
  G4double passageCellFlux;
  G4double noOfSecondary;
  G4double cellCharge;

public:
  inline void SetCellID(G4int x, G4int y, G4int z) {
    fXCellID = x;
    fYCellID = y;
    fZCellID = z;
  }
  inline G4int GetX() { return fXCellID; }
  inline G4int GetY() { return fYCellID; }
  inline G4int GetZ() { return fZCellID; }
  inline G4int GetID();
  inline void SetEdep(G4double de) { fEdep = de; }
  inline void AddEdep(G4double de) { fEdep += de; }
  inline G4double GetEdep() { return fEdep; }
  inline void SetDose(G4double dose) { fDose = dose; }
  inline void AddDose(G4double dose) { fDose += dose; }
  inline G4double GetDose() { return fDose; }
  inline void SetTrackLength(G4double tl) { fTrackLength = tl; }
  inline void AddTrackLength(G4double tl) { fTrackLength += tl; }
  inline G4double GetTrackLength() { return fTrackLength; }
  inline void SetParticleName(G4String pn) { fParticleName = pn; }
  inline G4String GetParticleName() { return fParticleName; }

  inline void SetFlatSurfaceCurrent(G4double v){ flatSurfCurr=v; }
  inline void SetPassageCellCurrent(G4double v){ passageCellCurr=v;}
  inline void SetFlatSurfaceFlux(G4double v){  flatSurfaceFlux=v; }
  inline void SetCellFlux(G4double v) {  cellFlux=v; }
  inline void SetPassageCellFlux(G4double v){  passageCellFlux=v; }
  inline void SetNofSecondary(G4double v){  noOfSecondary=v; }
  inline void SetCellCharge(G4double v){  cellCharge=v; }

  inline G4double GetFlatSurfaceCurrent(){ return flatSurfCurr; }
  inline G4double GetPassageCellCurrent(){ return passageCellCurr;}
  inline G4double GetFlatSurfaceFlux(){ return flatSurfaceFlux; }
  inline G4double GetCellFlux() { return cellFlux; }
  inline G4double GetPassageCellFlux(){ return passageCellFlux; }
  inline G4double GetNofSecondary(){ return noOfSecondary; }
  inline G4double GetCellCharge(){ return cellCharge; }
};

typedef G4THitsCollection<Test2PhantomHit> Test2PhantomHitsCollection;

extern G4Allocator<Test2PhantomHit> Test2PhantomHitAllocator;

inline void* Test2PhantomHit::operator new(size_t) {
  void *aHit;
  aHit = (void *) Test2PhantomHitAllocator.MallocSingle();
  return aHit;
}

inline void Test2PhantomHit::operator delete(void *aHit) {
  Test2PhantomHitAllocator.FreeSingle((Test2PhantomHit*) aHit);
}

inline G4int Test2PhantomHit::GetID(){
  return fXCellID*fSegment[1]*fSegment[2]+fYCellID*fSegment[2]+fZCellID;
}
#endif


