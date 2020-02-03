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
//
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  11 Jul 2003 A.Mantero, code cleaning / Plotter-XML addiction
//     Sep 2002 A.Mantero, AIDA3.0 Migration
//  06 Dec 2001 A.Pfeiffer updated for singleton
//  30 Nov 2001 Guy Barrand : migrate to AIDA-2.2.
//  28 Nov 2001 Elena Guardincerri   Created
//
// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"
#include "XrayFluoDataSet.hh"
#include "XrayFluoAnalysisMessenger.hh"

class G4Step;
class XrayFluoAnalysisMessenger;

//....oooOO0OOoo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class XrayFluoAnalysisManager
{
public:
 
  virtual ~XrayFluoAnalysisManager();
  
  void book();
  
  void finish();
  
  //fill histograms with data from XrayFluoSteppingAction
  void analyseStepping(const G4Step* aStep);
  
 //fill histograms with data from XrayFluoEventAction
  void analyseEnergyDep(G4double eDep);
  
 //fill histograms with data from XrayFluoPrimarygeneratorAction
  void analysePrimaryGenerator(G4double energy);
  
  //method to call to create an instance of this class
  static XrayFluoAnalysisManager* getInstance();

  // methods to set the flag for the storage of the space of phases into ntuple
  inline void PhaseSpaceOn(){phaseSpaceFlag = true;}

  inline void PhaseSpaceOff(){phaseSpaceFlag = false;}

  //method to chenge the name of the output file
  void SetOutputFileName(G4String);

  const std::pair<G4double,G4String> GetEmittedParticleEnergyAndType();

  void LoadGunData(G4String, G4bool);

  void SetPhysicFlag(G4bool);

private:
  //private constructor in order to create a singleton
  XrayFluoAnalysisManager();

  G4String outputFileName;

  G4bool phaseSpaceFlag;

  G4bool physicFlag;

  std::vector<G4double>* gunParticleEnergies;
  std::vector<G4String>* gunParticleTypes;

  //Instance for singleton implementation this is the returned 
  static XrayFluoAnalysisManager* instance;
  
  //pointer to the analysis messenger
  XrayFluoAnalysisMessenger* analisysMessenger;

  G4bool dataLoaded;
 
  G4int fParticleEnergyAndTypeIndex;

};
#endif



