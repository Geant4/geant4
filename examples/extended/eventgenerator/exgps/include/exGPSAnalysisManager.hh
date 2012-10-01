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
/// \file eventgenerator/exgps/include/exGPSAnalysisManager.hh
/// \brief Definition of the exGPSAnalysisManager class
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
  static exGPSAnalysisManager* GetInstance ();
  static void Dispose();

public:
  void BeginOfRun();
  void EndOfRun();

  void SetFileName(G4String filename) {fIleName = filename;};
  void SetFileType(G4String filetype) {fIleType = filetype;};

  void SetPosMax(G4double pmax) {fMaxpos = pmax;};
  void SetPosMin(G4double pmin) {fMinpos = pmin;};
  void SetEngMax(G4double emax) {fMaxeng = emax;};
  void SetEngMin(G4double emin) {fMineng = emin;};
  
  void Fill(G4String, G4double, G4double, G4double, G4double, G4double,
                                                                                                                                          G4double, G4double);

private:

  static exGPSAnalysisManager* fInstance;

  G4String fIleName;
  G4String fIleType;

  AIDA::IAnalysisFactory* fAnalysisFactory;
  AIDA::ITree* fTree;
  AIDA::IPlotter* fPlotter;

  G4double fMinpos, fMaxpos;
  G4double fMineng, fMaxeng;

  AIDA::IHistogram1D* fEnerHisto;
  AIDA::IHistogram2D* fPosiXY;
  AIDA::IHistogram2D* fPosiXZ;
  AIDA::IHistogram2D* fPosiYZ;
  AIDA::IHistogram2D* fAnglCTP;
  AIDA::IHistogram2D* fAnglTP;
  AIDA::ITuple* fTuple;

  exGPSAnalysisMessenger* fAnalysisMessenger;

};

#endif // G4ANALYSIS_USE

#endif











