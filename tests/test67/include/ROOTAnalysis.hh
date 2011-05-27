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
// $Id: ROOTAnalysis.hh,v 1.1 2009/03/21 18:37:27 vnivanch Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ROOTANALYSIS_h
#define ROOTANALYSIS_h 1

#include "globals.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TObjString.h"
#include "TGraphErrors.h"
#include <map>


class ROOTAnalysis
{
public:
  virtual ~ROOTAnalysis();
  void BookNewHistogram(G4int runID,G4double primaryEnergy);
  void AddEventEnergy(G4int runID,G4double energy);
  void AddSecondary(G4String processName);
  void CloseFile();
  void EndOfRun(G4int runID,G4double energy,G4int nEventsTotal,
		G4int nEventsPeak,G4int nEventsFull);

  //method to call to create an instance of this class
  static ROOTAnalysis* getInstance();
  void SetListName(G4String listName);

private:
   //private constructor in order to create a singleton
  ROOTAnalysis();
  static ROOTAnalysis* instance;
 
  TFile* fFile;
  std::map<G4int,TH1D*> *fHistograms;
  TObjString *fVersionName;
  TObjString *fPhysicsListName;
  TGraphErrors* fGraph;
  TGraphErrors* fGraph2;

  TTree* fTree;

  //Fields of the TTree
  //Counters for fluorescence gammas
  Double_t fPE; //photoelectric
  Double_t fIon; //ionisation
  Double_t fComp; //Compton
  Double_t fEfficiency;
  Double_t fEfficiencyErr;
  Double_t fPrimaryEnergy;
  Double_t fEfficiencyFull;
  Double_t fEfficiencyFullErr;

};


#endif


