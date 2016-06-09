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
/// \file analysis/A01/include/A01HodoscopeHit.hh
/// \brief Definition of the A01HodoscopeHit class
//
// $Id$
// --------------------------------------------------------------
//
#ifndef A01HodoscopeHit_h
#define A01HodoscopeHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class A01HodoscopeHit : public G4VHit
{
  public:
      A01HodoscopeHit(G4int i,G4double t);
      virtual ~A01HodoscopeHit();
      A01HodoscopeHit(const A01HodoscopeHit &right);
      const A01HodoscopeHit& operator=(const A01HodoscopeHit &right);
      int operator==(const A01HodoscopeHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void*aHit);

      void Draw();
      virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
      virtual std::vector<G4AttValue>* CreateAttValues() const;
      void Print();

  private:
      G4int fId;
      G4double fTime;
      G4ThreeVector fPos;
      G4RotationMatrix fRot;
      const G4LogicalVolume* fPLogV;

  public:
      inline G4int GetID() const { return fId; }
      inline G4double GetTime() const { return fTime; }
      inline void SetTime(G4double val) { fTime = val; }
      inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
      inline G4ThreeVector GetPos() const { return fPos; }
      inline void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
      inline G4RotationMatrix GetRot() const { return fRot; }
      inline void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
      inline const G4LogicalVolume* GetLogV() const { return fPLogV; }
};

typedef G4THitsCollection<A01HodoscopeHit> A01HodoscopeHitsCollection;

extern G4Allocator<A01HodoscopeHit> A01HodoscopeHitAllocator;

inline void* A01HodoscopeHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*)A01HodoscopeHitAllocator.MallocSingle();
  return aHit;
}

inline void A01HodoscopeHit::operator delete(void*aHit)
{
  A01HodoscopeHitAllocator.FreeSingle((A01HodoscopeHit*) aHit);
}

#endif


