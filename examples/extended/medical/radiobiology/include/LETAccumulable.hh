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
/// \file radiobiology/include/LETAccumulable.hh
/// \brief Definition of the RadioBio::LETAccumulable class

#ifndef RadiobiologyLETACCUMULABLE_HH
#define RadiobiologyLETACCUMULABLE_HH

#include "IonLet.hh"
#include "VRadiobiologicalAccumulable.hh"
#include <G4VAccumulable.hh>

#include <valarray>
#include <vector>

namespace RadioBio
{

// Forward declariation of other radiobiology classes
class Hit;
class VoxelizedSensitiveDetector;

/**
 * @brief Accumulable of LET-related data (that must be thread-local).
 *
 * It keeps the sum of alpha and beta numerators/denominator, as well as energy deposits.
 * The class is closely tied with the singLETon LET that is used both
 * to calculate alphas and betas, and also to store results.
 *
 * This is implemented as a customized G4VAccumulable with non-scalar data.
 *
 * @note There are two levels of merging (accumulation):
 *   1) From more threads in one run (G4VAccumulable merging is applied)
 *   2) (Optional) inter-run merging of data (implemented in LET).
 *
 * @note std::valarray is used (instead of C arrays or std::vectors)
 *    to accumulate data for its logical simplicity.
 */
class LETAccumulable : public VRadiobiologicalAccumulable
{
  public:
    LETAccumulable();
    LETAccumulable(const LETAccumulable& other) = default;

    // G4VAccumulable virtual methods to override
    void Merge(const G4VAccumulable& rhs) override;
    void Reset() override;

    // Store information from a single step
    void Accumulate(Hit* hit) override;

    // Type alias for numerical arrays
    using array_type = std::valarray<G4double>;

    // Access to stored data (to be called on the merged data)
    const array_type GetTotalLETT() const { return fTotalLETT; }
    const array_type GetTotalLETD() const { return fTotalLETD; }
    const array_type GetDTotalLETT() const { return fDTotalLETT; }
    const array_type GetDTotalLETD() const { return fDTotalLETD; }

    const std::vector<IonLet> GetIonLetStore() const { return fIonLetStore; }

    // Verbosity, shared with LET
    G4int GetVerboseLevel() const;

  private:
    // Apply configuration from the LET class and prepare matrices
    void Initialize();
    G4bool fInitialized{false};

    array_type fTotalLETT = {};
    array_type fTotalLETD = {};
    array_type fDTotalLETT = {};
    array_type fDTotalLETD = {};

    std::vector<IonLet> fIonLetStore = {};
};

}  // namespace RadioBio

#endif  // LETACCUMULABLE_HH
