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
/// \file GB03HodoHit.hh
/// \brief Definition of the GB03HodoHit class

/// The  hit class.
///
/// It defines data members to store the the kinetic energy,
/// Hit postion, particle PDG id:
/// - fEkin, fEncoding, fPos
///

#ifndef GB03HodoHit_h
#define GB03HodoHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class GB03HodoHit : public G4VHit
{
  public:
    GB03HodoHit() = default;
    GB03HodoHit(const GB03HodoHit&);
    ~GB03HodoHit() override = default;

    // operators
    const GB03HodoHit& operator=(const GB03HodoHit&);
    G4bool operator==(const GB03HodoHit&) const;

    inline void* operator new(size_t);
    inline void operator delete(void*);

    // methods from base class
    void Draw() override;
    void Print() override;

    // methods to handle data
    void Set(G4int id, G4double ek, G4double pw, G4ThreeVector xyz)
    {
      fEkin = ek;
      fWeight = pw;
      fEncoding = id;
      fPos = xyz;
    }
    void SetEkin(G4double ek) { fEkin = ek; }
    void SetWeight(G4double pw) { fWeight = pw; }
    void SetId(G4int id) { fEncoding = id; }
    void SetPos(G4ThreeVector xyz) { fPos = xyz; }

    // get methods
    G4double GetEkin() const { return fEkin; }
    G4double GetWeight() const { return fWeight; }
    G4int GetId() const { return fEncoding; }
    G4ThreeVector GetPos() const { return fPos; }

  private:
    G4double fEkin{0.0};  ///< Kinetic energy
    G4double fWeight{0.0};  ///< Track weight
    G4int fEncoding{0};  ///< PDG integer identifier
    G4ThreeVector fPos;  ///< Hit position
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using GB03HodoHitsCollection = G4THitsCollection<GB03HodoHit>;

extern G4ThreadLocal G4Allocator<GB03HodoHit>* GB03HodoHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* GB03HodoHit::operator new(size_t)
{
  if (GB03HodoHitAllocator == nullptr) {
    GB03HodoHitAllocator = new G4Allocator<GB03HodoHit>;
  }
  void* hit;
  hit = (void*)GB03HodoHitAllocator->MallocSingle();
  return hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void GB03HodoHit::operator delete(void* hit)
{
  GB03HodoHitAllocator->FreeSingle((GB03HodoHit*)hit);
}

#endif
