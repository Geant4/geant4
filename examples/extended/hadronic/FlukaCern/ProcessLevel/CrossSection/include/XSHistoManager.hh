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
///  \file XSHistoManager.hh
///  \brief Create a set of profiles for XS study.
//
//  Adapted from hadronic/Hadr00/include/HistoManager.hh
//  Author: G.Hugo, 06 January 2023
//
// ***************************************************************************
//
//      XSHistoManager
//
///  Create a set of profiles for XS study.
///
///  All profiles are G4H1. 
///  They are created and filled via G4VAnalysisManager.
///
///  The profiles can be dumped to all usual formats, including ROOT 
///  (via G4VAnalysisManager).
///  An interesting added feature here, is that the plots, while being allocated 
///  and filled via G4VAnalysisManager, are also dumped 
///  in a Flair-compatible format (via tools::histo::flair).
///
///  NB: tools::histo::flair code, which allows the dump of any G4H1 
///  into Flair-compatible format, is fully application-agnostic, 
///  and is placed in FlukaCern/utils. 
///  It could also be added as an extension of core G4 Analysis Manager.
//
// ***************************************************************************

#ifndef XS_HISTO_MANAGER_HH
#define XS_HISTO_MANAGER_HH

#include <unordered_map>

#include "globals.hh"

#include "G4SystemOfUnits.hh"

#include "XSHistoManagerMessenger.hh"


class G4ParticleDefinition;
class G4Element;
class G4Material;
class G4VAnalysisManager;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class XSHistoManager {

public:
  XSHistoManager();

  void SetOutputFileName(const G4String& outputFileName);
  void SetParticle(const G4String& particleName);
  void SetElement(const G4String& elementName);
  void SetMaterial(const G4String& materialName);
  void SetNumberOfBins(const G4int numBins) { fNumBins = numBins; }
  void SetMinKinEnergy(const G4double minKineticEnergy) { fMinKineticEnergy = minKineticEnergy; }
  void SetMaxKinEnergy(const G4double maxKineticEnergy) { fMaxKineticEnergy = maxKineticEnergy; }

  void Book();
  void EndOfRun();


private:
  void CheckInput();
  void DumpAllG4H1IntoRootFile() const;
  void DumpAllG4H1IntoFlairFile() const;

  XSHistoManagerMessenger* fMessenger = nullptr;

  G4String fOutputFileName = "all_XS";
  G4String fRootOutputFileName = "all_XS.root";
  G4String fFlairOutputFileName = "all_XS.hist";
  const G4ParticleDefinition* fParticle = nullptr;
  const G4Element* fElement = nullptr;
  const G4Material* fMaterial = nullptr;
  G4int fNumBins = 10000;
  G4double fMinKineticEnergy = 1.*keV;
  G4double fMaxKineticEnergy = 10.*TeV;
  G4String fFunctionName = "none";
  G4String fBinSchemeName = "log";
  G4String fRootEnergyUnit = "MeV";

  G4VAnalysisManager* fAnalysisManager = nullptr;
  std::unordered_map<G4int, G4int> fXSProfileIndex; // key: XS index
                                                // value: histo index from G4VAnalysisManager

  G4int fElasticXSIndex = 0;
  G4int fInelasticXSIndex = 1;
  G4int fCaptureXSIndex = 2;
  G4int fFissionXSIndex = 3;
  G4int fChargeExchangeXSIndex = 4;
  G4int fTotalXSIndex = 5;
  G4int fElasticPerVolumeXSIndex = 6;
  G4int fInelasticPerVolumeXSIndex = 7;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
