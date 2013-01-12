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
// $Id: XrayFluoAnalysisManager.hh
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  28 Nov 2001  Elena Guardincerri   Created
//
// -------------------------------------------------------------------



#ifdef G4ANALYSIS_USE
#ifndef XrayFluoAnalysisManager_h
#define XrayFluoAnalysisManager_h 1

#include "globals.hh"

// Histogramming from AIDA 
#include "Interfaces/IHistogram1D.h"
#include "Interfaces/IHistogram2D.h"

// Histogramming from Anaphe
#include "Interfaces/IHistoManager.h"

// Ntuples from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
#include "NtupleTag/LizardQuantity.h"


class G4Step;
class XrayFluoAnalysisMessenger;
class NTuple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
  
 //fill histograms with data from XrayFluoPrimarygeneratoraction
  void analysePrimaryGenerator(G4double energy);
  
  //method to call to create an instance of this class
  static XrayFluoAnalysisManager* getInstance();
 

private:
  
  //private constructor in order to create a singleton
  XrayFluoAnalysisManager();
 
  static XrayFluoAnalysisManager* instance;
  
  IHistoManager* histoManager;
 
  //nTuple factory
  Lizard::NTupleFactory* factory;
  Lizard::NTuple* ntuple;
  
  //pointer to the analysis messenger
  XrayFluoAnalysisMessenger* analysisMessenger;

  // Quantities for the ntuple
  Lizard::Quantity<float> eDep;
  Lizard::Quantity<float> counts;
};
#endif
#endif



