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
/// \file eventgenerator/HepMC/HepMCEx01/include/ExN04CalorimeterHit.hh
/// \brief Definition of the ExN04CalorimeterHit class
//
//

#ifndef ExN04CalorimeterHit_h
#define ExN04CalorimeterHit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"

class ExN04CalorimeterHit : public G4VHit {
public:
  ExN04CalorimeterHit();
  ExN04CalorimeterHit(G4LogicalVolume* logVol, G4int z, G4int phi);
  ExN04CalorimeterHit(const ExN04CalorimeterHit& right);

  virtual ~ExN04CalorimeterHit();

  const ExN04CalorimeterHit& operator=(const ExN04CalorimeterHit& right);
  G4bool operator==(const ExN04CalorimeterHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  virtual void Draw();
  virtual void Print();

  void SetCellID(G4int z,G4int phi);
  G4int GetZ();
  G4int GetPhi();
  void SetEdep(G4double de);
  void AddEdep(G4double de);
  G4double GetEdep();
  void SetPos(G4ThreeVector xyz);
  G4ThreeVector GetPos();
  void SetRot(G4RotationMatrix rmat);
  G4RotationMatrix GetRot();
  const G4LogicalVolume* GetLogV();

private:
  G4int fZCellID;
  G4int fPhiCellID;
  G4double fedep;
  G4ThreeVector fpos;
  G4RotationMatrix frot;
  const G4LogicalVolume* fpLogV;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void ExN04CalorimeterHit::SetCellID(G4int z,G4int phi)
 {
  fZCellID = z;
  fPhiCellID = phi;
}

inline G4int ExN04CalorimeterHit::GetZ()
{
  return fZCellID;
}

inline G4int ExN04CalorimeterHit::GetPhi()
{
  return fPhiCellID;
}

inline void ExN04CalorimeterHit::SetEdep(G4double de)
{
  fedep = de;
}

inline void ExN04CalorimeterHit::AddEdep(G4double de)
{
  fedep += de;
}

inline G4double ExN04CalorimeterHit::GetEdep()
{
  return fedep;
}

inline void ExN04CalorimeterHit::SetPos(G4ThreeVector xyz)
{
  fpos = xyz;
}

inline G4ThreeVector ExN04CalorimeterHit::GetPos()
{
  return fpos;
}

inline void ExN04CalorimeterHit::SetRot(G4RotationMatrix rmat)
{
  frot = rmat;
}

inline G4RotationMatrix ExN04CalorimeterHit::GetRot()
{
  return frot;
}

inline const G4LogicalVolume * ExN04CalorimeterHit::GetLogV()
{
  return fpLogV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
typedef G4THitsCollection<ExN04CalorimeterHit> ExN04CalorimeterHitsCollection;

extern G4Allocator<ExN04CalorimeterHit> ExN04CalorimeterHitAllocator;

inline void* ExN04CalorimeterHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) ExN04CalorimeterHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04CalorimeterHit::operator delete(void* aHit)
{
  ExN04CalorimeterHitAllocator.FreeSingle((ExN04CalorimeterHit*) aHit);
}

#endif
