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
// $Id: exrdmHisto.hh 68007 2013-03-13 11:28:03Z gcosmo $
//
/// \file radioactivedecay/rdecay02/include/exrdmHisto.hh
/// \brief Definition of the exrdmHisto class
//

#ifndef exrdmHisto_h
#define exrdmHisto_h 1

//---------------------------------------------------------------------------
//
// ClassName:   exrdmHisto
//
// Description: Utility class to hold and manipulate histograms/nTuples
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>

#ifdef G4ANALYSIS_USE
namespace AIDA {
 class IAnalysisFactory;
 class ITree;
 class ITuple;
 class IHistogram1D;
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE_ROOT
// Root classes
class TFile;
class TH1D;
class TNtuple;
#endif

class exrdmHistoMessenger;

class exrdmHisto
{

public:
  exrdmHisto();

  ~exrdmHisto();

  void Book();
  // Book predefined histogramms 

  void Save();
  // Save histogramms to file

  void Add1D(const G4String&, const G4String&, G4int nb=100, G4double x1=0., 
                                               G4double x2=1., G4double u=1.);
  // In this method histogramms are predefined

  void SetHisto1D(G4int, G4int, G4double, G4double, G4double);
  // It change fBins and boundaries

  void FillHisto(G4int, G4double, G4double);
  // exrdmHistogramms are filled

  void ScaleHisto(G4int, G4double);

  void AddTuple(const G4String&, const G4String&, const G4String&);
  // In this method fNTuple is booked

  void FillTuple(G4int, const G4String&, G4double);
  // Fill fNTuple parameter with a double

  void FillTuple(G4int, G4int, G4double);
  // Fill fNTuple at a given col with a double
  void FillTuple(G4int, const G4String&, G4String&);
  // Fill fNTuple parameter with a string

  void FillTuple(G4int, const G4String&, G4bool);
  // Fill fNTuple parameter with a bool

  void AddRow(G4int);
  // Save tuple event 

  void SetFileName(const G4String&);
  const G4String& GetFileName() const;  

  void SetFileType(const G4String&);
  const G4String& FileType() const;

private:

  G4String fHistName;
  G4String fHistType;

  G4int    fNHisto;
  G4int    fNTuple;
  G4int    fVerbose;
  G4int    fDefaultAct;
#ifdef G4ANALYSIS_USE
  std::vector<AIDA::IHistogram1D*> fHisto;
  std::vector<AIDA::ITuple*>   fNtup;
  AIDA::IAnalysisFactory* fAida;
  AIDA::ITree*    fTree;
#endif

#ifdef G4ANALYSIS_USE_ROOT
  TFile* fHfileROOT; 
  std::vector<TH1D*> fROOThisto;
  std::vector<TNtuple*>   fROOTntup;
  std::vector< std::vector<float> > fRarray;
  std::vector<G4int> fRcol;
#endif

  exrdmHistoMessenger* fMessenger;

  std::vector<G4int>     fActive;
  std::vector<G4int>     fBins;
  std::vector<G4double>  fXmin;
  std::vector<G4double>  fXmax;
  std::vector<G4double>  fUnit;
  std::vector<G4String>  fIds;
  std::vector<G4String>  fTitles;
  std::vector<G4String>  fTupleName;
  std::vector<G4String>  fTupleId;
  std::vector<G4String>  fTupleList;
  std::vector<G4String>  fTupleListROOT; 
};

#endif
