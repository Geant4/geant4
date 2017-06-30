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
// $Id: HistoManager.hh 104174 2017-05-15 12:12:45Z selles $
// GEANT4 tag $Name: geant4-09-02-ref-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HistoManager_h
#define HistoManager_h 1

#include "g4root.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int MaxHisto = 10;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:

  HistoManager(DetectorConstruction* det, G4double binLength);
  ~HistoManager();
  
  void Book();
  void Add1D(G4int histoId, const G4String& name, G4int nb, 
             G4double x1, G4double x2, 
             const G4String& u1="none", const G4String& u2="none");
  void SetHisto1D(G4int,G4int,G4double,G4double,const G4String& unit1="none");  
  
  void FillHisto(G4int id, G4double e, G4double weight = 1.0);
  void Scale (G4int, G4double);    
  void SetBinLength(G4double binLength);
  void Update(G4double binLength);
  
  G4bool    HistoExist  (G4int id) const { return fExist[id]; }
  G4double  GetHistoUnit(G4int id) const { return fUnit1[id]; }
  G4double  GetBinWidth (G4int id) const { return fWidth[id]; }
  
private:
  
  DetectorConstruction* fDetector;
  G4bool   fVerbose;
  G4double fBinLength;
  
  G4int                                 fNbHisto;
  std::vector<G4int>                    fHistoId;
  
  G4bool                                fNtupleActive;
  G4String fTupleName;
  G4String fTupleTitle;
  
  std::vector<G4bool>                   fExist;
  std::vector<G4String>                 fLabel;
  std::vector<G4String>                 fTitle;
  std::vector<G4int>                    fNbins;
  std::vector<G4double>                 fVmin ;
  std::vector<G4double>                 fVmax ;
  std::vector<G4double>                 fWidth;
  std::vector<G4double>                 fUnit1;
  std::vector<G4double>                 fUnit2;
  std::vector<G4String>                 fIds;
  std::vector<G4bool>                   fActive;
  
  std::vector<G4int>    fTupleI;
  std::vector<G4int>    fTupleF;
  std::vector<G4int>    fTupleD;
  std::vector<G4String> fNtupleI;
  std::vector<G4String> fNtupleF;
  std::vector<G4String> fNtupleD;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

