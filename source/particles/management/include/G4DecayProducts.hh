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
// G4DecayProducts
//
// Class description:
//
// Container for decay products of dynamic particles

// Author: H.Kurashige, 12 July 1996
// --------------------------------------------------------------------
#ifndef G4DecayProducts_hh
#define G4DecayProducts_hh 1

#include "G4DynamicParticle.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <vector>

class G4DecayProducts
{
  public:
    // Constructors
    G4DecayProducts();
    G4DecayProducts(const G4DynamicParticle& aParticle);

    // Copy constructor and assignment operator
    // Deep copy: for G4DynamicParticle pointer
    G4DecayProducts(const G4DecayProducts& right);
    G4DecayProducts& operator=(const G4DecayProducts& right);

    // Destructor
    ~G4DecayProducts();

    // Equality operators
    inline G4bool operator==(const G4DecayProducts& right) const;
    inline G4bool operator!=(const G4DecayProducts& right) const;

    // Set-get methods for the parent particle; theParentPaticle is used
    // to get information of parent particle when decay products are filled.
    // New G4DynamicParticle object is created in set methods
    inline const G4DynamicParticle* GetParentParticle() const;
    void SetParentParticle(const G4DynamicParticle& aParticle);

    // Boost all products
    void Boost(G4double totalEnergy, const G4ThreeVector& momentumDirection);
    void Boost(G4double betax, G4double betay, G4double betaz);

    // Push/pop methods for decay products pointer
    G4DynamicParticle* PopProducts();
    G4int PushProducts(G4DynamicParticle* aParticle);

    G4DynamicParticle* operator[](G4int anIndex) const;

    inline G4int entries() const { return numberOfProducts; }

    // Check energy/momentum of products
    G4bool IsChecked() const;

    void DumpInfo() const;

    using G4DecayProductVector = std::vector<G4DynamicParticle*>;

  private:
    G4int numberOfProducts = 0;
    G4DynamicParticle* theParentParticle = nullptr;
    G4DecayProductVector* theProductVector = nullptr;
};

// ------------------------
// Inline methods
// ------------------------

inline G4bool G4DecayProducts::operator==(const G4DecayProducts& right) const
{
  return (this == (G4DecayProducts*)&right);
}

inline G4bool G4DecayProducts::operator!=(const G4DecayProducts& right) const
{
  return (this != (G4DecayProducts*)&right);
}

inline const G4DynamicParticle* G4DecayProducts::GetParentParticle() const
{
  return theParentParticle;
}

#endif
