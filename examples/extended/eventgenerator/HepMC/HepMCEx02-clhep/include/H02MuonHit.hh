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
//   H02MuonHit.hh
//   $Id: H02MuonHit.hh,v 1.1 2002-11-19 10:35:49 murakami Exp $
//
// ====================================================================
#ifndef H02_MUON_HIT_H
#define H02_MUON_HIT_H
 
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class H02MuonHit : public G4VHit {
private:
  G4int moduleID;
  G4String pname;
  G4ThreeVector momentum;
  G4ThreeVector position;
  G4double tof;

public:
  H02MuonHit();
  H02MuonHit(G4int imod, G4String aname, const G4ThreeVector& pxyz,
	    const G4ThreeVector& xyz, G4double atof);
  ~H02MuonHit();

  H02MuonHit(const H02MuonHit& right);
  const H02MuonHit& operator=(const H02MuonHit& right);
  int operator==(const H02MuonHit& right) const;
  
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

inline void H02MuonHit::SetModuleID(G4int i) { moduleID=i; }
inline G4int H02MuonHit::GetModuleID() const { return moduleID; }

inline void H02MuonHit::SetParticle(G4String aname) { pname=aname; }
inline G4String H02MuonHit::GetParticle() const { return pname; }

inline void H02MuonHit::SetMomentum(const G4ThreeVector& pxyz) 
{ momentum=pxyz; }
inline G4ThreeVector H02MuonHit::GetMomentum() const { return momentum; }

inline void H02MuonHit::SetPosition(const G4ThreeVector& xyz) { position=xyz; }
inline G4ThreeVector H02MuonHit::GetPosition() const { return position; }

inline void H02MuonHit::SetTOF(G4double atof) { tof=atof; }
inline G4double H02MuonHit::GetTOF() const { return tof; }

typedef G4THitsCollection<H02MuonHit> H02MuonHitsCollection;
extern G4Allocator<H02MuonHit> H02MuonHitAllocator;

inline void* H02MuonHit::operator new(size_t)
{
  void* aHit;
  aHit= (void*)H02MuonHitAllocator.MallocSingle();
  return aHit;
}

inline void H02MuonHit::operator delete(void* aHit)
{
  H02MuonHitAllocator.FreeSingle((H02MuonHit*) aHit);
}

#endif
