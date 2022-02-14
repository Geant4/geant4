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
#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4String.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "CLHEP/Units/SystemOfUnits.h"

class G4ParticleGun;
class G4Event;
class G4Box;
class PrimaryGeneratorMessenger;

#ifdef WITHROOT
#include <RtypesCore.h>
class TFile;
class TTreeReader;
template <typename T> class TTreeReaderValue;
#endif

/**
 * @brief Primary generator
 *
 * Primary generator action.
 *
 * By default particle gun is used. It can be controlled with standard UI
 * commands (/gun/) and with additional ones introduced by the messenger.
 * "/HGCalTestbeam/generator/momentumSpread <VALUE>" to change constant
 * particle energy to Gaussian distribution with sigma expressed in units of the
 * initial energy (e.g. value 0.05 means sigma of 0.05 * E).
 * By default it equals to 0 and constant energy value is used.
 * "/HGCalTestbeam/generator/beamSpread <none/Gaussian/flat>" to define type of
 * beam position spread. By default none is used.
 * "/HGCalTestbeam/generator/beamSpreadX <SIZE>" to define size of beam spread
 * along x axis. It is sigma of a Gaussian distribution, or half-width of a
 * flat distribution.
 * "/HGCalTestbeam/generator/beamSpreadY <SIZE>" to define size of beam spread
 * along y axis. It is sigma of a Gaussian distribution, or half-width of a
 * flat distribution.
 * "/HGCalTestbeam/generator/fBeamZ0 <POSITION>" to define beam position along z
 * axis. By default edge of the world volume is used.
 *
 * If installation was done with ROOT package (CMake was able to locate it),
 * an additional option of input read from the ROOT file is enabled.
 * It can be activated with "/HGCalTestbeam/generator/fReadInputFile true".
 * "/HGCalTestbeam/generator/fPathInputFile <FILE>" sets the path to the input
 * file.
 * "/HGCalTestbeam/generator/startFromEvent <N>" allows to start simulation from
 * Nth event.
 * Please note that in current implementation input from file needs to be
 * executed in a non-multithreaded mode (or with 1 thread).
 *
 */

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event *);

  const G4ParticleGun *GetParticleGun() const { return fParticleGun; }

#ifdef WITHROOT
  /// Open input file with list of particles
  void OpenInput();

  /// Set flag indicating that particles should be read from file
  inline void SetIfUseInputFiles(G4bool aUseInputFiles) {
    fReadInputFile = aUseInputFiles;
  };
  /// Set the flag indicating that particles should be read from file
  inline G4bool GetIfUseInputFiles() { return fReadInputFile; }
  /// Set the path to the input file
  inline void SetInputFiles(G4String aInputFiles) {
    fPathInputFile = aInputFiles;
  };
  /// Get the path to the input file
  inline G4String GetInputFiles() const { return fPathInputFile; }
  /// Set ID of the first event to be read for the simulation
  inline void SetStartFromEvent(G4int aStartFromEvent) {
    fStartFromEvent = aStartFromEvent;
  };
  /// Get ID of the first event to be read for the simulation
  inline G4int GetStartFromEvent() const { return fStartFromEvent; }
#endif
  /// Set sigma of the Gaussian distribution for the momentum spread
  /// @param[in] aMomentumSpread sigma of Gaussian distribution expressed in
  /// units of initial energy (e.g. 0.05 means sigma = 0.05 * E)
  inline void SetMomentumSpread(G4double aMomentumSpread) {
    fMomentumGaussianSpread = aMomentumSpread;
  };
  /// Get sigma of the Gaussian distribution for the momentum spread
  inline G4double GetMomentumSpread() const { return fMomentumGaussianSpread; }
  /// Set type of beam position spread
  /// @param[in] aType Type of beam position spread: "none", "Gaussian" or
  /// "flat". By default "none" is used.
  inline void SetBeamSpreadType(G4String aType) {
    if (aType == "none")
      fBeamType = eNone;
    if (aType == "Gaussian")
      fBeamType = eGaussian;
    if (aType == "flat")
      fBeamType = eFlat;
  }
  /// Get type of beam position spread
  inline G4String GetBeamSpreadType() const {
    switch (fBeamType) {
    case eNone:
      return "none";
      break;
    case eGaussian:
      return "Gaussian";
      break;
    case eFlat:
      return "flat";
      break;
    }
    return "";
  }
  /// Set size of beam position spread along X axis
  /// @param[in] aBeamSpreadX Size of beam position spread. Sigma for Gaussian
  /// distribution or half-width of flat distribution.
  inline void SetBeamSpreadX(G4double aBeamSpreadX) {
    fSigmaBeamX = aBeamSpreadX;
  }
  /// Get size of beam position spread along X axis
  inline G4double GetBeamSpreadX() const { return fSigmaBeamX; }
  /// Set size of beam position spread along Y axis
  /// @param[in] aBeamSpreadY Size of beam position spread. Sigma for Gaussian
  /// distribution or half-width of flat distribution.
  inline void SetBeamSpreadY(G4double aBeamSpreadY) {
    fSigmaBeamY = aBeamSpreadY;
  }
  /// Get size of beam position spread along Y axis
  inline G4double GetBeamSpreadY() const { return fSigmaBeamY; }
  /// Set initial beam position along Z axis
  /// By default edge of world volume is used
  inline void SetBeamZ0(G4double aBeamZ0) { fBeamZ0 = aBeamZ0; }
  /// Get initial beam position along Z axis
  inline G4double GetBeamZ0() const { return fBeamZ0; }

private:
  /// Pointer to the particle gun
  G4ParticleGun *fParticleGun;
  /// Pointer to the world volume for initial beam position
  G4Box *fEnvelopeBox;
  /// Pointer to the messenger with custom UI commands
  PrimaryGeneratorMessenger *fMessenger;
  /// enum describing the beam position spread in transverse plane
  enum eBeamType { eNone, eGaussian, eFlat };
  /// Type of beam position spread in transverse dimension
  eBeamType fBeamType = eBeamType::eNone;
  /// Size of beam position spread along X axis
  /// Sigma for Gaussian, and half-width for flat distribution
  G4double fSigmaBeamX = 0;
  /// Size of beam position spread along Y axis
  /// Sigma for Gaussian, and half-width for flat distribution
  G4double fSigmaBeamY = 0;
  /// Initial beam position along Z axis
  G4double fBeamZ0 = -999 * CLHEP::m;
  /// Sigma of Gaussian momentum spread
  G4double fMomentumGaussianSpread = 0;

#ifdef WITHROOT
  /// Flag indicating if primaries should be read from file instead of using
  /// the particle gun
  G4bool fReadInputFile = false;
  /// Path to the input file
  G4String fPathInputFile = "";
  /// ID of the first event in the file to be used in this simulation
  G4int fStartFromEvent = 0;
  /// Counter of event
  G4int fEventCounter = -1;
  /// Pointer to the input file
  TFile *fInputFile = nullptr;
  /// Pointer to the tree containing particles
  TTreeReader *fHgcalReader = nullptr;
  /// Reader of event ID
  TTreeReaderValue<Float_t> *fHgcalEventId;
  /// Reader of particle PDG
  TTreeReaderValue<Float_t> *fHgcalPdgId;
  /// Reader of particle X position (in mm)
  TTreeReaderValue<Float_t> *fHgcalPosX;
  /// Reader of particle Y position (in mm)
  TTreeReaderValue<Float_t> *fHgcalPosY;
  /// Reader of particle Z position (in mm)
  // TTreeReaderValue<Float_t> *fHgcalPosZ;
  /// Reader of particle X momentum (in MeV)
  TTreeReaderValue<Float_t> *fHgcalMomX;
  /// Reader of particle Y momentum (in MeV)
  TTreeReaderValue<Float_t> *fHgcalMomY;
  /// Reader of particle Z momentum (in MeV)
  TTreeReaderValue<Float_t> *fHgcalMomZ;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
