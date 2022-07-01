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
/// \file molecularDNAHit.hh
/// \brief Hit class for an event that interacts with DNA
#ifndef MOLECULAR_DNA_HIT_HH
#define MOLECULAR_DNA_HIT_HH

#include <utility>

#include "G4VHit.hh"

#include "MoleculeList.hh"

#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4MolecularConfiguration;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DNAHit : public G4VHit
{
 public:
  // Main constructors
  DNAHit();

  DNAHit(const molecule&, const G4int&, const G4int&, const G4int&,
         const int64_t&,
         G4ThreeVector,                                              // dousatsu
         G4ThreeVector, const G4double&, const G4double&, G4String,  // dousatsu
         const G4MolecularConfiguration*);                           // dousatsu

  // Constructor delegation
  // Physical Hit
  DNAHit(const molecule& mol, const G4int& placement_idx, const G4int& chain,
         const G4int& strand, const int64_t& bp, const G4ThreeVector& pos,
         const G4ThreeVector& localpos, G4double energy, G4double d,
         G4String chromo)
    : DNAHit(mol, placement_idx, chain, strand, bp, pos, localpos, energy, d,
             std::move(chromo), nullptr){};

  // Chemical hit
  DNAHit(const molecule& mol, const G4int& placement_idx, const G4int& chain,
         const G4int& strand, const int64_t& bp, const G4ThreeVector& pos,
         const G4ThreeVector& localpos, G4String chromo,
         const G4MolecularConfiguration* radical)
    : DNAHit(mol, placement_idx, chain, strand, bp, pos, localpos, 0, 0,
             std::move(chromo), radical){};

  DNAHit(const DNAHit&);

  ~DNAHit() override;

  // Add Method, not as an operator to prevent people making mistakes
  void AddHit(const DNAHit&);

  // operators
  inline void* operator new(size_t);

  inline void operator delete(void*);

  const DNAHit& operator=(const DNAHit&);

  G4int operator==(const DNAHit&) const;

  // setters
  // Avoid setters in general and use the constructor, as it also sets
  // the computed quantites
  void SetMolecule(const molecule& mol) { fMoleculeEnum = mol; };

  void SetPlacementIdx(const G4int& place_idx) { fPlacementIdx = place_idx; };

  void SetChainIdx(const G4int& chainidx) { fChainIdx = chainidx; };

  void SetStrandIdx(const G4int& strandidx) { fStrandIdx = strandidx; };

  void SetBasePairIdx(const int64_t& bpidx)
  {
    fBasePairIdx = bpidx;
  };  // dousatsu

  void SetPosition(const G4ThreeVector& p) { fPosition = p; };

  void SetLocalPosition(const G4ThreeVector& p) { fLocalPosition = p; };

  void SetEnergy(const G4double& energy) { fEnergy = energy; };

  void SetDistance(const G4double& dist) { fDistance = dist; };

  void SetChromosome(const G4String& chrom) { fChromosome = chrom; };

  void SetRadical(const G4MolecularConfiguration* r) { fRadical = r; };

  // getters
  molecule GetMolecule() const { return fMoleculeEnum; };

  G4int GetPlacementIdx() const { return fPlacementIdx; };

  G4int GetChainIdx() const { return fChainIdx; };

  G4int GetStrandIdx() const { return fStrandIdx; };

  int64_t GetBasePairIdx() const
  {
    return fBasePairIdx;
  };  // dousatsu
  G4ThreeVector GetPosition() const { return fPosition; };

  G4ThreeVector GetLocalPosition() const
  {
    return fLocalPosition;
  };

  G4double GetEnergy() const { return fEnergy; };

  G4double GetDistance() const { return fDistance; };

  G4String GetChromosome() const { return fChromosome; };

  const G4MolecularConfiguration* GetRadical() const
  {
    return fRadical;
  };

  // Get Computed Quantites
  const G4MolecularConfiguration* GetStrand1Rad() const
  {
    return fStrand1Rad;
  };

  const G4MolecularConfiguration* GetBase1Rad() const
  {
    return fBase1Rad;
  };

  const G4MolecularConfiguration* GetStrand2Rad() const
  {
    return fStrand2Rad;
  };

  const G4MolecularConfiguration* GetBase2Rad() const
  {
    return fBase2Rad;
  };

  G4double GetStrand1Energy() const { return fStrand1Energy; };

  G4double GetStrand2Energy() const { return fStrand2Energy; };

  G4double GetBP1Energy() const { return fBP1Energy; };

  G4double GetBP2Energy() const { return fBP2Energy; };

 private:
  molecule fMoleculeEnum = UNSPECIFIED;
  G4int fPlacementIdx = -1;   // ORG
  G4int fChainIdx = -1;       // ORG
  G4int fStrandIdx = -1;      // ORG
  int64_t fBasePairIdx = -1;  // dousatsu
  G4ThreeVector fPosition = G4ThreeVector();
  G4ThreeVector fLocalPosition = G4ThreeVector();
  G4double fEnergy = 0.;
  G4double fDistance = 0.;
  G4String fChromosome = "";
  const G4MolecularConfiguration* fRadical = nullptr;
  const G4MolecularConfiguration* fStrand1Rad = nullptr;
  const G4MolecularConfiguration* fBase1Rad = nullptr;
  const G4MolecularConfiguration* fStrand2Rad = nullptr;
  const G4MolecularConfiguration* fBase2Rad = nullptr;

  G4double fStrand1Energy = 0.;
  G4double fStrand2Energy = 0.;
  G4double fBP1Energy = 0.;
  G4double fBP2Energy = 0.;
};

// typedefs
using MolecularDNAHitsCollection = G4THitsCollection<DNAHit>;

// memory management
extern G4ThreadLocal G4Allocator<DNAHit>* MolecularDNAHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* DNAHit::operator new(size_t)
{
  if(MolecularDNAHitAllocator == nullptr)
  {
    MolecularDNAHitAllocator = new G4Allocator<DNAHit>;
  }
  return (void*) MolecularDNAHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DNAHit::operator delete(void* hit)
{
  MolecularDNAHitAllocator->FreeSingle((DNAHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DNA_HIT_HH
