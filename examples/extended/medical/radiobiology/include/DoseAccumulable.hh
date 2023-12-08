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
/// \file radiobiology/include/DoseAccumulable.hh
/// \brief Definition of the RadioBio::DoseAccumulable class

#ifndef RadiobiologyDoseACCUMULABLE_HH
#define RadiobiologyDoseACCUMULABLE_HH

#include "VRadiobiologicalAccumulable.hh"
#include <G4VAccumulable.hh>

#include <valarray>
#include <vector>

namespace RadioBio
{

// Forward declariation of other radiobiology classes
class VoxelizedSensitiveDetector;

/**
 * @brief Accumulable of Dose-related data (that must be thread-local).
 *
 * It keeps the energy deposit data divided in voxels for dose calculation
 *
 * This is implemented as a customized G4VAccumulable with non-scalar data.
 *
 * @note There are two levels of merging (accumulation):
 *   1) From more threads in one run (G4VAccumulable merging is applied)
 *   2) (Optional) inter-run merging of data (implemented in Dose).
 *
 * @note std::valarray is used (instead of C arrays or std::vectors)
 *    to accumulate data for its logical simplicity.
 */
class DoseAccumulable : public VRadiobiologicalAccumulable
{
  public:
    DoseAccumulable();
    DoseAccumulable(const DoseAccumulable& other) = default;

    // G4VAccumulable virtual methods to override
    void Merge(const G4VAccumulable& rhs) override;
    void Reset() override;

    // Store information from a single step
    void Accumulate(Hit* hit) override;

    // Type alias for numerical arrays
    using array_type = std::valarray<G4double>;

    // Access to stored data (to be called on the merged data)
    const array_type GetEnDeposit() const { return fEnDep; }

    // Verbosity, shared with Dose
    G4int GetVerboseLevel() const;

  private:
    // Apply configuration from the Dose class and prepare matrices
    void Initialize();
    G4bool fInitialized = false;

    array_type fEnDep;
};

}  // namespace RadioBio

#endif  // DoseACCUMULABLE_HH
