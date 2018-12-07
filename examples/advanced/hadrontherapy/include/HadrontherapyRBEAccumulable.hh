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
// HadrontherapyRBEAccumulable.hh;
//
#ifndef HADRONTHERAPYRBEACCUMULABLE_HH
#define HADRONTHERAPYRBEACCUMULABLE_HH

#include <G4VAccumulable.hh>

#include <valarray>

/**
 * @brief Accumulable of RBE-related data (that must be thread-local).
 *
 * It keeps the sum of alpha and beta numerators/denominator, as well as energy deposits.
 * The class is closely tied with the singleton HadrontherapyRBE that is used both
 * to calculate alphas and betas, and also to store results.
 *
 * This is implemented as a customized G4VAccumulable with non-scalar data.
 *
 * @note There are two levels of merging (accumulation):
 *   1) From more threads in one run (G4VAccumulable merging is applied)
 *   2) (Optional) inter-run merging of data (implemented in HadrontherapyRBE).
 *
 * @note std::valarray is used (instead of C arrays or std::vectors)
 *    to accumulate data for its logical simplicity.
 */
class HadrontherapyRBEAccumulable : public G4VAccumulable
{
public:
    HadrontherapyRBEAccumulable();
    HadrontherapyRBEAccumulable(const HadrontherapyRBEAccumulable& other) = default;

    // G4VAccumulable virtual methods
    void Merge(const G4VAccumulable &rhs) override;
    void Reset() override;

    // Store information from a single step
    void Accumulate(G4double E, G4double energyDeposit, G4double dX, G4int Z, G4int i, G4int j, G4int k);

    // Type alias for numerical arrays
    using array_type = std::valarray<G4double>;

    // Access to stored data (to be called on the merged data)
    const array_type GetEnergyDeposit() const;
    const array_type GetAlphaNumerator() const { return fAlphaNumerator; }
    const array_type GetBetaNumerator() const { return fBetaNumerator; }
    const array_type GetDenominator() const { return fDenominator; }

    // Verbosity, shared with HadrontherapyRBE
    G4int GetVerboseLevel() const;

private:
    /** @brief Helper function to get the 1D index in the 3D array */
    inline G4int GetIndex(G4int i, G4int j, G4int k) const {return (i * fVoxelsAlongY + j) * fVoxelsAlongZ + k; }

    // Apply configuration from the HadrontherapyRBE class and prepare matrices
    void Initialize();
    G4bool fInitialized { false };

    array_type fAlphaNumerator;
    array_type fBetaNumerator;
    array_type fDenominator;
    array_type fEnergyDeposit;

    // How many voxels do we have?
    // ...along each axis
    size_t fVoxelsAlongX;
    size_t fVoxelsAlongY;
    size_t fVoxelsAlongZ;
    // ...and in total
    size_t fVoxels;
};

#endif // HADRONTHERAPYRBEACCUMULABLE_HH
