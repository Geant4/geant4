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
/// \file eventgenerator/HepMC/HepMCEx01/include/ExN04MuonHit.hh
/// \brief Definition of the ExN04MuonHit class
//
//

#ifndef ExN04MuonHit_h
#define ExN04MuonHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class ExN04MuonHit : public G4VHit {
public:
  ExN04MuonHit();
  ~ExN04MuonHit();
  ExN04MuonHit(const ExN04MuonHit& right);
  const ExN04MuonHit& operator=(const ExN04MuonHit &right);
  G4bool operator==(const ExN04MuonHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  virtual void Draw();
  virtual void Print();

  inline void SetEdep(G4double de) { fedep = de; }
  inline void AddEdep(G4double de) { fedep += de; }
  inline G4double GetEdep() { return fedep; }
  inline void SetPos(G4ThreeVector xyz) { fpos = xyz; }
  inline G4ThreeVector GetPos() { return fpos; }

private:
  G4double fedep;
  G4ThreeVector fpos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
typedef G4THitsCollection<ExN04MuonHit> ExN04MuonHitsCollection;

extern G4Allocator<ExN04MuonHit> ExN04MuonHitAllocator;

inline void* ExN04MuonHit::operator new(size_t)
{
  void* aHit;
  aHit = (void *) ExN04MuonHitAllocator.MallocSingle();
  return aHit;
}

inline void ExN04MuonHit::operator delete(void *aHit)
{
  ExN04MuonHitAllocator.FreeSingle((ExN04MuonHit*) aHit);
}

#endif
