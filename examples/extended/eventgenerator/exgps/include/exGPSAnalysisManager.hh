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
#ifndef exGPSAnalysisManager_h
#define exGPSAnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"

#ifdef G4ANALYSIS_USE_AIDA
namespace AIDA{
class IAnalysisFactory;
class ITree;
class IHistogramFactory;
class ITupleFactory;
class IPlotter;
}
#endif

#ifdef G4ANALYSIS_USE_ROOT
class TFile;
class TH1D;
class TH2D;
class TNtuple;
#endif

class exGPSAnalysisMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSAnalysisManager
{

private:
  exGPSAnalysisManager ();

public:
  virtual ~exGPSAnalysisManager ();
  static exGPSAnalysisManager* getInstance ();
  static void dispose();
  
#ifdef G4ANALYSIS_USE_AIDA
  AIDA::IHistogramFactory* getHistogramFactory();
  AIDA::ITupleFactory* getTupleFactory();
  AIDA::IPlotter* createPlotter();
#endif

public:
  void BeginOfRun();
  void EndOfRun();

  void SetFileName(G4String filename) {fileName = filename;};
  void SetFileType(G4String filetype) {fileType = filetype;};

  void SetPosMax(G4double pmax) {maxpos = pmax;};
  void SetPosMin(G4double pmin) {minpos = pmin;};
  void SetEngMax(G4double emax) {maxeng = emax;};
  void SetEngMin(G4double emin) {mineng = emin;};
  
  void Fill(G4double, G4double, G4double, G4double, G4double, G4double, G4double, G4double);

private:

  static exGPSAnalysisManager* instance;

  G4String fileName;
  G4String fileType;

  G4double minpos, maxpos;
  G4double mineng, maxeng;
  
 #ifdef G4ANALYSIS_USE_AIDA 
  IAnalysisFactory* analysisFactory;
  IHistogramFactory* hFactory;
  ITupleFactory* tFactory;
  ITree* tree;
  IHistogram1D* enerHisto;
  IHistogram2D* posiXY;
  IHistogram2D* posiXZ;
  IHistogram2D* posiYZ;
  IHistogram2D* anglCTP;
  IHistogram2D* anglTP;
  ITuple* tuple;

  IPlotter* plotter;
#endif

#ifdef G4ANALYSIS_USE_ROOT
  TFile* hfileroot; // the file for histograms, tree ...
  TH1D* enerHistoroot;
  TH2D* posiXYroot;
  TH2D* posiXZroot;
  TH2D* posiYZroot;
  TH2D* anglCTProot;
  TH2D* anglTProot;
  TNtuple* tupleroot; 	
#endif

  exGPSAnalysisMessenger* analysisMessenger;

};

#endif // G4ANALYSIS_USE

#endif











