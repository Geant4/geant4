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
// GEANT4 Class file
//
// File name:     G4OpticalParameters
//
// Author:        Daren Sawkey based on G4EmParameters
//
// Creation date: 14.07.2020
//
// Modifications:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4OpticalParameters.hh"
#include "G4OpticalParametersMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"

G4OpticalParameters* G4OpticalParameters::theInstance = nullptr;

#ifdef G4MULTITHREADED
G4Mutex G4OpticalParameters::opticalParametersMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4OpticalParameters* G4OpticalParameters::Instance()
{
  if(nullptr == theInstance)
  {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&opticalParametersMutex);
    if(nullptr == theInstance)
    {
#endif
      static G4OpticalParameters manager;
      theInstance = &manager;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&opticalParametersMutex);
#endif
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4OpticalParameters::~G4OpticalParameters() { delete theMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4OpticalParameters::G4OpticalParameters()
{
  theMessenger = new G4OpticalParametersMessenger(this);
  Initialise();

  fStateManager = G4StateManager::GetStateManager();
}

void G4OpticalParameters::SetDefaults()
{
  if(!IsLocked())
  {
    Initialise();
  }
}

void G4OpticalParameters::Initialise()
{
  verboseLevel = 0;

  cerenkovStackPhotons          = true;
  cerenkovTrackSecondariesFirst = true;
  cerenkovVerboseLevel          = 0;
  cerenkovMaxPhotons            = 100;
  cerenkovMaxBetaChange         = 10.;

  scintByParticleType        = false;
  scintTrackInfo             = false;
  scintStackPhotons          = true;
  scintFiniteRiseTime        = false;
  scintTrackSecondariesFirst = true;
  scintVerboseLevel          = 0;

  wlsTimeProfileName = "delta";
  wlsVerboseLevel    = 0;

  wls2TimeProfileName = "delta";
  wls2VerboseLevel    = 0;

  absorptionVerboseLevel = 0;

  rayleighVerboseLevel = 0;

  mieVerboseLevel = 0;

  boundaryInvokeSD     = false;
  boundaryVerboseLevel = 0;

  processActivation["OpRayleigh"]    = true;
  processActivation["OpBoundary"]    = true;
  processActivation["OpMieHG"]       = true;
  processActivation["OpAbsorption"]  = true;
  processActivation["OpWLS"]         = true;
  processActivation["OpWLS2"]        = true;
  processActivation["Cerenkov"]      = true;
  processActivation["Scintillation"] = true;
}

void G4OpticalParameters::SetVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  verboseLevel = val;
  SetCerenkovVerboseLevel(verboseLevel);
  SetScintVerboseLevel(verboseLevel);
  SetRayleighVerboseLevel(verboseLevel);
  SetAbsorptionVerboseLevel(verboseLevel);
  SetMieVerboseLevel(verboseLevel);
  SetBoundaryVerboseLevel(verboseLevel);
  SetWLSVerboseLevel(verboseLevel);
  SetWLS2VerboseLevel(verboseLevel);
}

G4int G4OpticalParameters::GetVerboseLevel() const { return verboseLevel; }

void G4OpticalParameters::SetProcessActivation(const G4String& process,
                                               G4bool val)
{
  // Configure the physics constructor to use/not use a selected process.
  // This method can only be called in PreInit> phase (before execution of
  // ConstructProcess). The process is not added to particle's process manager
  // and so it cannot be re-activated later in Idle> phase with the command
  // /process/activate.

  if(IsLocked())
  {
    return;
  }
  if(processActivation[process] == val)
    return;

  // processActivation keys defined at initialisation
  if(processActivation.find(process) != processActivation.end())
  {
    processActivation[process] = val;
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "Process name " << process << " out of bounds.";
    G4Exception("G4OpticalParameters::SetProcessActivation()", "Optical013",
                FatalException, ed);
  }
}

G4bool G4OpticalParameters::GetProcessActivation(const G4String& process) const
{
  return processActivation.find(process)->second;
}

void G4OpticalParameters::SetCerenkovStackPhotons(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  cerenkovStackPhotons = val;
}

G4bool G4OpticalParameters::GetCerenkovStackPhotons() const
{
  return cerenkovStackPhotons;
}

void G4OpticalParameters::SetCerenkovVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  cerenkovVerboseLevel = val;
}

G4int G4OpticalParameters::GetCerenkovVerboseLevel() const
{
  return cerenkovVerboseLevel;
}

void G4OpticalParameters::SetCerenkovMaxPhotonsPerStep(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  cerenkovMaxPhotons = val;
}

G4int G4OpticalParameters::GetCerenkovMaxPhotonsPerStep() const
{
  return cerenkovMaxPhotons;
}

void G4OpticalParameters::SetCerenkovMaxBetaChange(G4double val)
{
  if(IsLocked())
  {
    return;
  }
  cerenkovMaxBetaChange = val;
}

G4double G4OpticalParameters::GetCerenkovMaxBetaChange() const
{
  return cerenkovMaxBetaChange;
}

void G4OpticalParameters::SetCerenkovTrackSecondariesFirst(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  cerenkovTrackSecondariesFirst = val;
}

G4bool G4OpticalParameters::GetCerenkovTrackSecondariesFirst() const
{
  return cerenkovTrackSecondariesFirst;
}

void G4OpticalParameters::SetScintByParticleType(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  scintByParticleType = val;
}

G4bool G4OpticalParameters::GetScintByParticleType() const
{
  return scintByParticleType;
}

void G4OpticalParameters::SetScintTrackInfo(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  scintTrackInfo = val;
}

G4bool G4OpticalParameters::GetScintTrackInfo() const { return scintTrackInfo; }

void G4OpticalParameters::SetScintTrackSecondariesFirst(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  scintTrackSecondariesFirst = val;
}

G4bool G4OpticalParameters::GetScintTrackSecondariesFirst() const
{
  return scintTrackSecondariesFirst;
}

void G4OpticalParameters::SetScintFiniteRiseTime(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  scintFiniteRiseTime = val;
}

G4bool G4OpticalParameters::GetScintFiniteRiseTime() const
{
  return scintFiniteRiseTime;
}

void G4OpticalParameters::SetScintStackPhotons(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  scintStackPhotons = val;
}

G4bool G4OpticalParameters::GetScintStackPhotons() const
{
  return scintStackPhotons;
}

void G4OpticalParameters::SetScintVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  scintVerboseLevel = val;
}

G4int G4OpticalParameters::GetScintVerboseLevel() const
{
  return scintVerboseLevel;
}

void G4OpticalParameters::SetWLSTimeProfile(const G4String& val)
{
  if(IsLocked())
  {
    return;
  }
  wlsTimeProfileName = val;
}

G4String G4OpticalParameters::GetWLSTimeProfile() const
{
  return wlsTimeProfileName;
}

void G4OpticalParameters::SetWLSVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  wlsVerboseLevel = val;
}

G4int G4OpticalParameters::GetWLSVerboseLevel() const
{
  return wlsVerboseLevel;
}

void G4OpticalParameters::SetWLS2TimeProfile(const G4String& val)
{
  if(IsLocked())
  {
    return;
  }
  wls2TimeProfileName = val;
}

G4String G4OpticalParameters::GetWLS2TimeProfile() const
{
  return wls2TimeProfileName;
}

void G4OpticalParameters::SetWLS2VerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  wls2VerboseLevel = val;
}

G4int G4OpticalParameters::GetWLS2VerboseLevel() const
{
  return wls2VerboseLevel;
}

void G4OpticalParameters::SetBoundaryVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  boundaryVerboseLevel = val;
}

G4int G4OpticalParameters::GetBoundaryVerboseLevel() const
{
  return boundaryVerboseLevel;
}

void G4OpticalParameters::SetBoundaryInvokeSD(G4bool val)
{
  if(IsLocked())
  {
    return;
  }
  boundaryInvokeSD = val;
}

G4bool G4OpticalParameters::GetBoundaryInvokeSD() const
{
  return boundaryInvokeSD;
}

void G4OpticalParameters::SetAbsorptionVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  absorptionVerboseLevel = val;
}

G4int G4OpticalParameters::GetAbsorptionVerboseLevel() const
{
  return absorptionVerboseLevel;
}

void G4OpticalParameters::SetRayleighVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  rayleighVerboseLevel = val;
}

G4int G4OpticalParameters::GetRayleighVerboseLevel() const
{
  return rayleighVerboseLevel;
}

void G4OpticalParameters::SetMieVerboseLevel(G4int val)
{
  if(IsLocked())
  {
    return;
  }
  mieVerboseLevel = val;
}

G4int G4OpticalParameters::GetMieVerboseLevel() const
{
  return mieVerboseLevel;
}

void G4OpticalParameters::PrintWarning(G4ExceptionDescription& ed) const
{
  G4Exception("G4EmParameters", "Optical0020", JustWarning, ed);
}

void G4OpticalParameters::StreamInfo(std::ostream& os) const
{
  G4long prec = os.precision(5);
  os
    << "======================================================================="
    << "\n";
  os
    << "======                         Optical Physics Parameters      ========"
    << "\n";
  os
    << "======================================================================="
    << "\n";

  os << " Cerenkov process active:               "
     << GetProcessActivation("Cerenkov") << "\n";
  os << " Cerenkov maximum photons per step:     " << cerenkovMaxPhotons
     << "\n";
  os << " Cerenkov maximum beta change per step: " << cerenkovMaxBetaChange
     << " %\n";
  os << " Cerenkov stack photons:                " << cerenkovStackPhotons
     << "\n";
  os << " Cerenkov track secondaries first:      "
     << cerenkovTrackSecondariesFirst << "\n";
  os << " Scintillation process active:          "
     << GetProcessActivation("Scintillation") << "\n";
  os << " Scintillation finite rise time:        " << scintFiniteRiseTime
     << "\n";
  os << " Scintillation by particle type:        " << scintByParticleType
     << "\n";
  os << " Scintillation record track info:       " << scintTrackInfo << "\n";
  os << " Scintillation stack photons:           " << scintStackPhotons << "\n";
  os << " Scintillation track secondaries first: " << scintTrackSecondariesFirst
     << "\n";
  os << " WLS process active:                    "
     << GetProcessActivation("OpWLS") << "\n";
  os << " WLS time profile name:                 " << wlsTimeProfileName
     << "\n";
  os << " WLS2 process active:                   "
     << GetProcessActivation("OpWLS2") << "\n";
  os << " WLS2 time profile name:                " << wls2TimeProfileName
     << "\n";
  os << " Boundary process active:               "
     << GetProcessActivation("OpBoundary") << "\n";
  os << " Boundary invoke sensitive detector:    " << boundaryInvokeSD << "\n";
  os << " Rayleigh process active:               "
     << GetProcessActivation("OpRayleigh") << "\n";
  os << " MieHG process active:                  "
     << GetProcessActivation("OpMieHG") << "\n";
  os << " Absorption process active:             "
     << GetProcessActivation("OpAbsorption") << "\n";
  os
    << "======================================================================="
    << "\n";
  os.precision(prec);
}

void G4OpticalParameters::Dump() const
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&opticalParametersMutex);
#endif
  StreamInfo(G4cout);
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&opticalParametersMutex);
#endif
}

std::ostream& operator<<(std::ostream& os, const G4OpticalParameters& par)
{
  par.StreamInfo(os);
  return os;
}

G4bool G4OpticalParameters::IsLocked() const
{
  return (!G4Threading::IsMasterThread() ||
          (fStateManager->GetCurrentState() != G4State_PreInit &&
           fStateManager->GetCurrentState() != G4State_Init &&
           fStateManager->GetCurrentState() != G4State_Idle));
}
