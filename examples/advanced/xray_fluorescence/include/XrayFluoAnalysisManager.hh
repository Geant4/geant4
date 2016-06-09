//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include "AIDA/AIDA.h" // Headers for AIDA interfaces
#include "XrayFluoAnalysisMessenger.hh"

class G4Step;
class XrayFluoAnalysisMessenger;

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
  
 //fill histograms with data from XrayFluoPrimarygeneratorAction
  void analysePrimaryGenerator(G4double energy);
  
  //method to call to create an instance of this class
  static XrayFluoAnalysisManager* getInstance();

  //method to create hbook or xml file for persistency
  
  void CreatePersistency(G4String fileName,G4String persistencyType,
		       G4bool readOnly = false, G4bool createNew = true);

  inline  void CreatePersistency() {CreatePersistency(outputFileName,persistencyType);}

  //method to chenge the name of the output file
  void SetOutputFileName(G4String);

  //method to chenge the type of the output file
  void SetOutputFileType(G4String);
 
  // method used by the messenger 
  G4bool GetDeletePersistencyFileFlag();

  // methods used by RunManager and EvenManager to visualize partial results
  //  void InitializePlotter();

  //  void PlotCurrentResults();

private:
  //private constructor in order to create a singleton
  XrayFluoAnalysisManager();

  G4String outputFileName;

  //  G4bool visPlotter;

  G4String persistencyType;

  G4bool deletePersistencyFile;

  //Instance for singleton implementation this is the returned 
  static XrayFluoAnalysisManager* instance;
  
  //pointer to the analysis messenger
  XrayFluoAnalysisMessenger* analisysMessenger;

  //XrayFluoEventAction* pEvent;

  // analysis data members 

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::IHistogramFactory* histogramFactory;
  //  AIDA::IPlotterFactory* plotterFactory;
  //  AIDA::IPlotter* plotter;

  AIDA::IHistogram1D*   histo_1;
  AIDA::IHistogram1D*   histo_2;
  AIDA::IHistogram1D*   histo_3;
  AIDA::IHistogram1D*   histo_4;
  AIDA::IHistogram1D*   histo_5;
  AIDA::IHistogram1D*   histo_6;
  AIDA::IHistogram1D*   histo_7;
  AIDA::IHistogram1D*   histo_8;
  AIDA::IHistogram1D*   histo_9;
  //AIDA::IHistogram1D*   histo_10; //
  //AIDA::IHistogram1D*   histo_12; // Created for debuggig purpose
  //AIDA::IHistogram1D*   histo_11; //
};

#endif



