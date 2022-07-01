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
#ifdef USE_INFERENCE
#ifndef PAR04INFEERENCESETUP_HH
#define PAR04INFEERENCESETUP_HH

#include <G4String.hh>                   // for G4String
#include <G4SystemOfUnits.hh>            // for mm
#include <G4Types.hh>                    // for G4int, G4double, G4bool, G4f...
#include <memory>                        // for unique_ptr
#include <vector>                        // for vector
#include "CLHEP/Units/SystemOfUnits.h"   // for mm
#include "G4ThreeVector.hh"              // for G4ThreeVector
class Par04DetectorConstruction;
class Par04InferenceInterface;
class Par04InferenceMessenger;

/**
 * @brief Inference setup.
 *
 * Constructs the input vector of size b+c to run the inference, b represents the size of
 * the latent space (or the encoded space in a Variational Autoencoder based model),
 * c represents the size of the conditional vector. The b values of the input vector
 * are randomly sampled from b-dimensional Gaussian distribution. The c values
 * represent respectively the condition values of the particle energy, angle and
 * detector geometry. These condition values are user-specific application.
 * The energy rescaling is used to retrieve the original energy scale in MeV.
 * Computes the cell position in the detector of each inferred energy value.
 *
 **/

class Par04InferenceSetup
{
 public:
  Par04InferenceSetup();
  ~Par04InferenceSetup();

  /// Geometry setup
  /// Check if inference should be performed for the particle
  /// @param[in] aEnergy Particle's energy
  G4bool IfTrigger(G4double aEnergy);
  /// Set mesh size.
  /// @param aSize (x,y,x) size for Carthesian coordinates, or (R, phi, z) for
  /// cylindrical coordinates.
  inline void SetMeshSize(const G4ThreeVector& aSize) { fMeshSize = aSize; };
  /// Get mesh size.
  /// @return G4ThreeVector (x,y,x) size for Carthesian coordinates, or (R, phi,
  /// z) for cylindrical coordinates.
  inline G4ThreeVector GetMeshSize() const { return fMeshSize; };
  /// Set number of mesh cells.
  /// @param aSize (x,y,x) size for Carthesian coordinates, or (R, phi, z) for
  /// cylindrical coordinates.
  inline void SetMeshNumber(const G4ThreeVector& aSize) { fMeshNumber = aSize; };
  /// Get number of mesh cells.
  /// @return G4ThreeVector (x,y,x) size for Carthesian coordinates, or (R, phi,
  /// z) for cylindrical coordinates.
  inline G4ThreeVector GetMeshNumber() const { return fMeshNumber; };
  /// Set size of the condition vector
  inline void SetSizeConditionVector(G4int aNumber) { fSizeConditionVector = aNumber; };
  /// Get size of the condition vector
  inline G4int GetSizeConditionVector() const { return fSizeConditionVector; };
  /// Set size of the latent space vector
  inline void SetSizeLatentVector(G4int aNumber) { fSizeLatentVector = aNumber; };
  /// Get size of the latent space vector
  inline G4int GetSizeLatentVector() const { return fSizeLatentVector; };
  /// Set path and name of the model
  inline void SetModelPathName(G4String aName) { fModelPathName = aName; };
  /// Get path and name of the model
  inline G4String GetModelPathName() const { return fModelPathName; };
  /// Set profiling flag
  inline void SetProfileFlag(G4int aNumber) { fProfileFlag = aNumber; };
  /// Get profiling flag
  inline G4int GetProfileFlag() const { return fProfileFlag; };
  /// Set optimization flag
  inline void SetOptimizationFlag(G4int aNumber) { fOptimizationFlag = aNumber; };
  /// Get optimization flag
  inline G4int GetOptimizationFlag() const { return fOptimizationFlag; };
  /// Get name of the inference library
  inline G4String GetInferenceLibrary() const { return fInferenceLibrary; };
  /// Set name of the inference library and create a pointer to chosen inference interface
  void SetInferenceLibrary(G4String aName);
  /// Check settings of the inference library
  void CheckInferenceLibrary();
  /// Set number of Mesh cells in cylindrical coordinates (r, phi, z)
  inline void SetMeshNbOfCells(G4ThreeVector aNb) { fMeshNumber = aNb; };
  /// Set number of Mesh cells in cylindrical coordinates
  /// @param[in] aIndex index of cylindrical axis (0,1,2) = (r, phi, z)
  inline void SetMeshNbOfCells(G4int aIndex, G4double aNb) { fMeshNumber[aIndex] = aNb; };
  /// Get number of Mesh cells in cylindrical coordinates (r, phi, z)
  inline G4ThreeVector GetMeshNbOfCells() const { return fMeshNumber; };
  /// Set size of Mesh cells in cylindrical coordinates (r, phi, z)
  inline void SetMeshSizeOfCells(G4ThreeVector aNb) { fMeshSize = aNb; };
  /// Set size of Mesh cells in cylindrical coordinates
  /// @param[in] aIndex index of cylindrical axis (0,1,2) = (r, phi, z)
  inline void SetMeshSizeOfCells(G4int aIndex, G4double aNb) { fMeshSize[aIndex] = aNb; };
  /// Get size of Mesh cells in cylindrical coordinates (r, phi, z)
  inline G4ThreeVector GetMeshSizeOfCells() const { return fMeshSize; };

  /// Execute inference
  /// @param[out] aDepositsEnergies of inferred energies deposited in the
  /// detector
  /// @param[in] aParticleEnergy Energy of initial particle
  void GetEnergies(std::vector<G4double>& aEnergies, G4double aParticleEnergy,
                   G4float aInitialAngle);

  /// Calculate positions
  /// @param[out] aDepositsPositions Vector of positions corresponding to
  /// energies deposited in the detector
  /// @param[in] aParticlePosition Initial particle position which is centre of
  /// transverse plane of the mesh
  ///            and beginning of the mesh in the longitudinal direction
  /// @param[in] aParticleDirection Initial particle direction for the mesh
  /// rotation
  void GetPositions(std::vector<G4ThreeVector>& aDepositsPositions, G4ThreeVector aParticlePosition,
                    G4ThreeVector aParticleDirection);

 private:
  /// Cell's size: (x,y,x) for Carthesian, and (R, phi, z) for cylindrical
  /// coordinates Can be changed with UI command `/example/mesh/size <x y z>/<r
  /// phi z> <unit>`. For cylindrical coordinates phi is ignored and calculated
  /// from fMeshNumber.
  G4ThreeVector fMeshSize = G4ThreeVector(2.325 * CLHEP::mm, 1, 3.4 * CLHEP::mm);
  /// Number of cells: (x,y,x) for Carthesian, and (R, phi, z) for cylindrical
  /// coordinates. Can be changed with UI command `/example/mesh/number <Nx Ny
  /// Nz>/<Nr Nphi Nz>`
  G4ThreeVector fMeshNumber = G4ThreeVector(18, 50, 45);
  /// Inference interface
  std::unique_ptr<Par04InferenceInterface> fInferenceInterface;
  /// Inference messenger
  Par04InferenceMessenger* fInferenceMessenger;
  /// Maximum particle energy value (in MeV) in the training range
  float fMaxEnergy = 1024000.0;
  /// Maximum particle angle (in degrees) in the training range
  float fMaxAngle = 90.0;
  /// Name of the inference library
  G4String fInferenceLibrary = "ONNX";
  /// Size of the latent space vector
  G4int fSizeLatentVector = 10;
  /// Size of the condition vector
  G4int fSizeConditionVector = 4;
  /// Name of the inference library
  G4String fModelPathName = "MLModels/Generator.onnx";
  /// ONNX specific
  /// Profiling flag
  G4bool fProfileFlag = false;
  /// Optimization flag
  G4bool fOptimizationFlag = false;
  /// Intra-operation number of threads
  G4int fIntraOpNumThreads = 1;
};

#endif /* PAR04INFEERENCESETUP_HH */
#endif
