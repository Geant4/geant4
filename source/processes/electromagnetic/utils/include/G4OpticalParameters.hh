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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4OpticalParameters
//
// Author:        Daren Sawkey based on G4EmParameters
//
// Creation date: 14.07.2020
//
// Modifications:
//
//
// Class Description:
//
// A utility static class, responsable for keeping parameters
// for all optical physics processes and models.
//
// It is initialized by the master thread but can be updated
// at any moment. Parameters may be used in run time or at
// initialisation
//
// -------------------------------------------------------------------
//

#ifndef G4OpticalParameters_h
#define G4OpticalParameters_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include <vector>

class G4OpticalParametersMessenger;
class G4StateManager;

enum G4OpticalProcessIndex
{
  kCerenkov,       ///< Cerenkov process index
  kScintillation,  ///< Scintillation process index
  kAbsorption,     ///< Absorption process index
  kRayleigh,       ///< Rayleigh scattering process index
  kMieHG,          ///< Mie scattering process index
  kBoundary,       ///< Boundary process index
  kWLS,            ///< Wave Length Shifting process index
  kWLS2,           ///< Second Wave Length Shifting process index
  kNoProcess       ///< Number of processes, no selected process
};

/// Return the name for a given optical process index
G4String G4OpticalProcessName(G4int);

inline G4String G4OpticalProcessName(G4int processNumber)
{
  switch(processNumber)
  {
    case kCerenkov:
      return "Cerenkov";
    case kScintillation:
      return "Scintillation";
    case kAbsorption:
      return "OpAbsorption";
    case kRayleigh:
      return "OpRayleigh";
    case kMieHG:
      return "OpMieHG";
    case kBoundary:
      return "OpBoundary";
    case kWLS:
      return "OpWLS";
    case kWLS2:
      return "OpWLS2";
    default:
      return "NoProcess";
  }
}

class G4OpticalParameters
{
 public:
  static G4OpticalParameters* Instance();

  ~G4OpticalParameters();

  void SetDefaults();

  // printing
  void StreamInfo(std::ostream& os) const;
  void Dump() const;
  friend std::ostream& operator<<(std::ostream& os, const G4OpticalParameters&);

  void  SetVerboseLevel(G4int);
  G4int GetVerboseLevel() const;

  void   SetProcessActivation(const G4String&, G4bool);
  G4bool GetProcessActivation(const G4String&) const;

  // Cerenkov
  void  SetCerenkovMaxPhotonsPerStep(G4int);
  G4int GetCerenkovMaxPhotonsPerStep() const;
  void  SetCerenkovVerboseLevel(G4int);
  G4int GetCerenkovVerboseLevel() const;
  void     SetCerenkovMaxBetaChange(G4double);
  G4double GetCerenkovMaxBetaChange() const;
  void   SetCerenkovTrackSecondariesFirst(G4bool);
  G4bool GetCerenkovTrackSecondariesFirst() const;
  void   SetCerenkovStackPhotons(G4bool);
  G4bool GetCerenkovStackPhotons() const;

  // Scintillation
  void   SetScintByParticleType(G4bool);
  G4bool GetScintByParticleType() const;
  void   SetScintTrackInfo(G4bool);
  G4bool GetScintTrackInfo() const;
  void   SetScintTrackSecondariesFirst(G4bool);
  G4bool GetScintTrackSecondariesFirst() const;
  void   SetScintFiniteRiseTime(G4bool);
  G4bool GetScintFiniteRiseTime() const;
  void   SetScintStackPhotons(G4bool);
  G4bool GetScintStackPhotons() const;
  void   SetScintVerboseLevel(G4int);
  G4int  GetScintVerboseLevel() const;
  void   SetScintEnhancedTimeConstants(G4bool);
  G4bool GetScintEnhancedTimeConstants() const;

  // WLS
  void     SetWLSTimeProfile(const G4String&);
  G4String GetWLSTimeProfile() const;
  void  SetWLSVerboseLevel(G4int);
  G4int GetWLSVerboseLevel() const;

  // WLS2
  void     SetWLS2TimeProfile(const G4String&);
  G4String GetWLS2TimeProfile() const;
  void  SetWLS2VerboseLevel(G4int);
  G4int GetWLS2VerboseLevel() const;

  // boundary
  void  SetBoundaryVerboseLevel(G4int);
  G4int GetBoundaryVerboseLevel() const;
  void   SetBoundaryInvokeSD(G4bool);
  G4bool GetBoundaryInvokeSD() const;

  // absorption
  void  SetAbsorptionVerboseLevel(G4int);
  G4int GetAbsorptionVerboseLevel() const;

  // rayleigh
  void  SetRayleighVerboseLevel(G4int);
  G4int GetRayleighVerboseLevel() const;

  // mie
  void  SetMieVerboseLevel(G4int);
  G4int GetMieVerboseLevel() const;

 private:
  G4OpticalParameters();
  void Initialise();
  G4bool IsLocked() const;
  void PrintWarning(G4ExceptionDescription& ed) const;

  static G4OpticalParameters* theInstance;

  G4OpticalParametersMessenger* theMessenger;
  G4StateManager* fStateManager;

  G4int verboseLevel;

  // Whether to activate each process
  std::map<G4String, G4bool> processActivation;

  // cerenkov/////////////////
  G4bool cerenkovStackPhotons;
  G4bool cerenkovTrackSecondariesFirst;
  G4int cerenkovVerboseLevel;
  G4int cerenkovMaxPhotons;
  G4double cerenkovMaxBetaChange;

  // scintillation /////////////////

  /// option to set a finite rise-time; Note: the G4Scintillation
  /// process expects the user to have set the constant material
  /// property SCINTILLATIONRISETIME{1,2,3}
  G4bool scintFiniteRiseTime;

  /// option to  allow for the light yield to be a function of
  /// particle type and deposited energy in case of non-linear
  /// light emission in scintillators
  G4bool scintByParticleType;

  /// option to allow for G4ScintillationTrackInformation
  /// to be attached to a scintillation photon's track
  G4bool scintTrackInfo;

  /// option to allow stacking of secondary Scintillation photons
  G4bool scintStackPhotons;

  G4int scintVerboseLevel;
  G4bool scintTrackSecondariesFirst;

  ///////////////// WLS
  G4String wlsTimeProfileName;
  G4int wlsVerboseLevel;

  ///////////////// WLS2
  G4String wls2TimeProfileName;
  G4int wls2VerboseLevel;

  //////////////// absorption
  G4int absorptionVerboseLevel;

  //////////////// rayleigh
  G4int rayleighVerboseLevel;

  //////////////// mie
  G4int mieVerboseLevel;

  //////////////// boundary
  /// G4OpBoundaryProcess to call InvokeSD method
  G4bool boundaryInvokeSD;
  G4int boundaryVerboseLevel;

#ifdef G4MULTITHREADED
  static G4Mutex opticalParametersMutex;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
