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
/// \file Par01/include/Par01CalorimeterHit.hh
/// \brief Definition of the Par01CalorimeterHit class
//
//
//

#ifndef Par01CalorimeterHit_h
#define Par01CalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class Par01CalorimeterHit : public G4VHit
{
public:
  
  Par01CalorimeterHit();
  Par01CalorimeterHit(G4LogicalVolume* logVol);
  ~Par01CalorimeterHit();
  Par01CalorimeterHit(const Par01CalorimeterHit &right);
  const Par01CalorimeterHit& operator=(const Par01CalorimeterHit &right);
  G4bool operator==(const Par01CalorimeterHit &right) const;
  
  
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
  void operator delete(void *,void*){}
#endif
  
  // methods from base class
  virtual void Draw();
  virtual void Print();
  
private:
  G4double fEdep;
  G4ThreeVector fPosition;
  G4RotationMatrix fRot;
  const G4LogicalVolume* fLogV;
  
public:
  inline void SetEdep(G4double de)
  { fEdep = de; };
  inline void AddEdep(G4double de)
  { fEdep += de; };
  inline G4double GetEdep()
  { return fEdep; };
  inline void SetPos(G4ThreeVector xyz)
  { fPosition = xyz; };
  inline G4ThreeVector GetPos()
  { return fPosition; };
  inline void SetRot(G4RotationMatrix rmat)
  { fRot = rmat; };
  inline G4RotationMatrix GetRot()
  { return fRot; };
  inline const G4LogicalVolume * GetLogV()
  { return fLogV; };
  
};

typedef G4THitsCollection<Par01CalorimeterHit> Par01CalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<Par01CalorimeterHit>* Par01CalorimeterHitAllocator;

inline void* Par01CalorimeterHit::operator new(size_t)
{
  if(!Par01CalorimeterHitAllocator) Par01CalorimeterHitAllocator =
                                      new G4Allocator<Par01CalorimeterHit>;
  return (void *) Par01CalorimeterHitAllocator->MallocSingle();
}

inline void Par01CalorimeterHit::operator delete(void *aHit)
{
  Par01CalorimeterHitAllocator->FreeSingle((Par01CalorimeterHit*) aHit);
}

#endif
