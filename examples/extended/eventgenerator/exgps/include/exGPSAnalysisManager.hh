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
#ifndef exGPSAnalysisManager_h
#define exGPSAnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"

using namespace AIDA;

class AIDA::IAnalysisFactory;
class AIDA::ITree;
class AIDA::IHistogramFactory;
class AIDA::ITupleFactory;
class AIDA::IPlotter;

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

  IHistogramFactory* getHistogramFactory();
  ITupleFactory* getTupleFactory();
  IPlotter* createPlotter();


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

  IAnalysisFactory* analysisFactory;
  IHistogramFactory* hFactory;
  ITupleFactory* tFactory;
  ITree* tree;

  G4double minpos, maxpos;
  G4double mineng, maxeng;

  IHistogram1D* enerHisto;
  IHistogram2D* posiXY;
  IHistogram2D* posiXZ;
  IHistogram2D* posiYZ;
  IHistogram2D* anglCTP;
  IHistogram2D* anglTP;
  ITuple* tuple;

  IPlotter* plotter;

  exGPSAnalysisMessenger* analysisMessenger;

};

#endif // G4ANALYSIS_USE

#endif











