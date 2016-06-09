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

namespace AIDA {

class IAnalysisFactory;
class ITree;
class IPlotter;
class IHistogram1D;
class IHistogram2D;
class ITuple;
}

class exGPSAnalysisMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSAnalysisManager
{

private:
  exGPSAnalysisManager ();
  virtual ~exGPSAnalysisManager ();

public:
  static exGPSAnalysisManager* getInstance ();
  static void dispose();

public:
  void BeginOfRun();
  void EndOfRun();

  void SetFileName(G4String filename) {fileName = filename;};
  void SetFileType(G4String filetype) {fileType = filetype;};

  void SetPosMax(G4double pmax) {maxpos = pmax;};
  void SetPosMin(G4double pmin) {minpos = pmin;};
  void SetEngMax(G4double emax) {maxeng = emax;};
  void SetEngMin(G4double emin) {mineng = emin;};
  
  void Fill(G4String, G4double, G4double, G4double, G4double, G4double, G4double, G4double);

private:

  static exGPSAnalysisManager* instance;

  G4String fileName;
  G4String fileType;

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::IPlotter* plotter;

  G4double minpos, maxpos;
  G4double mineng, maxeng;

  AIDA::IHistogram1D* enerHisto;
  AIDA::IHistogram2D* posiXY;
  AIDA::IHistogram2D* posiXZ;
  AIDA::IHistogram2D* posiYZ;
  AIDA::IHistogram2D* anglCTP;
  AIDA::IHistogram2D* anglTP;
  AIDA::ITuple* tuple;

  exGPSAnalysisMessenger* analysisMessenger;

};

#endif // G4ANALYSIS_USE

#endif











