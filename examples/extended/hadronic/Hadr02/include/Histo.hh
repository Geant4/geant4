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
/// \file hadronic/Hadr02/include/Histo.hh
/// \brief Definition of the Histo class
//

#ifndef Histo_h
#define Histo_h 1

//---------------------------------------------------------------------------
//
// ClassName:   Histo
//
// Description: Utility class to hold and manipulate histograms/nTuples
//
// Author:      V.Ivanchenko 30/10/03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DataVector.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VAnalysisManager;
class HistoMessenger;

class Histo
{
public:

  Histo();

  ~Histo();

  // Book predefined histogramms 
  void Book();

  // Save histogramms to file
  void Save();

  // In this method 1-D histogramms are predefined
  void Add1D(const G4String&, const G4String&, G4int nb, G4double x1, 
            G4double x2, G4double u=1.);

  // It change bins and boundaries
  void SetHisto1D(G4int, G4int, G4double, G4double, G4double);

  // Histogram activation/deactivation
  void Activate(G4int, G4bool);

  // Histogramms are filled
  void Fill(G4int, G4double, G4double);

  // Histogramms are scaled
  void ScaleH1(G4int, G4double);

  // In this method nTuple is booked
  void AddTuple(const G4String&);

  // In this method nTuple is booked
  void AddTupleI(const G4String&);
  void AddTupleF(const G4String&);
  void AddTupleD(const G4String&);

  // Fill nTuple parameter
  void FillTupleI(G4int, G4int);
  void FillTupleF(G4int, G4float);
  void FillTupleD(G4int, G4double);

  // Save tuple event 
  void AddRow();

  // Set output file
  void SetFileName(const G4String&);
  void SetFileType(const G4String&);

  inline void SetVerbose(G4int val) { fVerbose = val; };

  inline G4bool IsActive() const { return fHistoActive; };

private:

  G4VAnalysisManager* fManager;
  HistoMessenger* fMessenger;
 
  G4String fHistName;
  G4String fHistType;
  G4String fTupleName;
  G4String fTupleTitle;
  G4int    fNHisto;
  G4int    fVerbose;
  G4bool   fDefaultAct;
  G4bool   fHistoActive;
  G4bool   fNtupleActive;

  std::vector<G4int>    fHisto;
  std::vector<G4int>    fTupleI;
  std::vector<G4int>    fTupleF;
  std::vector<G4int>    fTupleD;
  std::vector<G4int>    fBins;
  std::vector<G4bool>   fActive;
  std::vector<G4double> fXmin;
  std::vector<G4double> fXmax;
  std::vector<G4double> fUnit;
  std::vector<G4String> fIds;
  std::vector<G4String> fTitles;
  std::vector<G4String> fNtupleI;
  std::vector<G4String> fNtupleF;
  std::vector<G4String> fNtupleD;

};

#endif
