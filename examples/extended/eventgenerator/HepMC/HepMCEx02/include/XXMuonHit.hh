//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

// ====================================================================
//
//   XXMuonHit.hh
//   $Id: XXMuonHit.hh,v 1.1 2002-04-29 20:44:39 asaim Exp $
//
// ====================================================================
#ifndef XX_MUON_HIT_H
#define XX_MUON_HIT_H
 
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class XXMuonHit : public G4VHit {
private:
  G4int moduleID;
  G4String pname;
  G4ThreeVector momentum;
  G4ThreeVector position;
  G4double tof;

public:
  XXMuonHit();
  XXMuonHit(G4int imod, G4String aname, const G4ThreeVector& pxyz,
	    const G4ThreeVector& xyz, G4double atof);
  ~XXMuonHit();

  XXMuonHit(const XXMuonHit& right);
  const XXMuonHit& operator=(const XXMuonHit& right);
  int operator==(const XXMuonHit& right) const;
  
  void* operator new(size_t);
  void operator delete(void* aHit);

  // set/get functions...
  void SetModuleID(G4int i);
  G4int GetModuleID() const;

  void SetParticle(G4String aname);
  G4String GetParticle() const;

  void SetMomentum(const G4ThreeVector& pxyz);
  G4ThreeVector GetMomentum() const;

  void SetPosition(const G4ThreeVector& xyz);
  G4ThreeVector GetPosition() const;

  void SetTOF(G4double atof);
  G4double GetTOF() const;

  // methods...
  virtual void Draw();
  virtual void Print();  
};

// ====================================================================
// inline functions
// ====================================================================

inline void XXMuonHit::SetModuleID(G4int i) { moduleID=i; }
inline G4int XXMuonHit::GetModuleID() const { return moduleID; }

inline void XXMuonHit::SetParticle(G4String aname) { pname=aname; }
inline G4String XXMuonHit::GetParticle() const { return pname; }

inline void XXMuonHit::SetMomentum(const G4ThreeVector& pxyz) 
{ momentum=pxyz; }
inline G4ThreeVector XXMuonHit::GetMomentum() const { return momentum; }

inline void XXMuonHit::SetPosition(const G4ThreeVector& xyz) { position=xyz; }
inline G4ThreeVector XXMuonHit::GetPosition() const { return position; }

inline void XXMuonHit::SetTOF(G4double atof) { tof=atof; }
inline G4double XXMuonHit::GetTOF() const { return tof; }

typedef G4THitsCollection<XXMuonHit> XXMuonHitsCollection;
extern G4Allocator<XXMuonHit> XXMuonHitAllocator;

inline void* XXMuonHit::operator new(size_t)
{
  void* aHit;
  aHit= (void*)XXMuonHitAllocator.MallocSingle();
  return aHit;
}

inline void XXMuonHit::operator delete(void* aHit)
{
  XXMuonHitAllocator.FreeSingle((XXMuonHit*) aHit);
}

#endif
