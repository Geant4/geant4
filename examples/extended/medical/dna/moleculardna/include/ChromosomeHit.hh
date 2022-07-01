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
//
/// \file molecularChromosomeHit.hh
/// \brief Hit class for an event that interacts with DNA
#ifndef MOLECULAR_CHROMOSOME_HIT_HH
#define MOLECULAR_CHROMOSOME_HIT_HH

#include "G4VHit.hh"

#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ChromosomeHit : public G4VHit
{
 public:
  explicit ChromosomeHit(G4String  key);

  ChromosomeHit(const ChromosomeHit&);

  ~ChromosomeHit() override;

  // operators
  inline void* operator new(size_t);

  inline void operator delete(void*);

  const ChromosomeHit& operator=(const ChromosomeHit&);

  G4int operator==(const ChromosomeHit&) const;

  // Adders
  void AddChromosomeEdep(G4double edep) { fEdepChromosome += edep; };

  void AddDNAEdep(G4double edep) { fEdepDNA += edep; };

  // Getters
  G4double GetChromosomeEdep() const { return fEdepChromosome; };

  G4double GetDNAEdep() const { return fEdepDNA; };

  const G4String& GetName() const { return fName; };

 private:
  G4String fName;
  G4double fEdepChromosome = 0.;
  G4double fEdepDNA = 0.;
};

using MolecularChromosomeHitsCollection = G4THitsCollection<ChromosomeHit>;

extern G4ThreadLocal G4Allocator<ChromosomeHit>*
  MolecularChromosomeHitAllocator;

inline void* ChromosomeHit::operator new(size_t)
{
  if(MolecularChromosomeHitAllocator == nullptr)
  {
    MolecularChromosomeHitAllocator = new G4Allocator<ChromosomeHit>;
  }
  return (void*) MolecularChromosomeHitAllocator->MallocSingle();
}

inline void ChromosomeHit::operator delete(void* hit)
{
  MolecularChromosomeHitAllocator->FreeSingle((ChromosomeHit*) hit);
}

#endif  // MOLECULAR_CHROMOSOME_HIT_HH
