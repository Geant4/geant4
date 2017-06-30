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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
//#include "G4ThreeVector.hh"

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class RunAction : public G4UserRunAction
{
public:
  
  RunAction(const DetectorConstruction*);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
    
  void  SetRndmFreq(G4int val) {fSaveRndm = val;}
  G4int GetRndmFreq()    const {return fSaveRndm;}

  void AddDoseN(G4double dose) {fDoseN += dose;}
  void SetDoseN(G4double dose) {fDoseN = dose;}
  G4double GetDoseN()    const {return fDoseN;}

  void AddDoseC(G4double dose) {fDoseC += dose;}
  void SetDoseC(G4double dose) {fDoseC = dose;}
  G4double GetDoseC()    const {return fDoseC;}

  G4int GetNumEvent()    const {return fNumEvent;}
  void SetNumEvent(G4int i)    {fNumEvent = i;}

  G4int GetNbOfHitsGas() const {return fNbOfHitsGas;}
  void AddNbOfHitsGas()        {fNbOfHitsGas = fNbOfHitsGas+1;}

  void SetMassNucleus(G4double mN) {fMassNucleus = mN;}
  G4double GetMassNucleus()  const {return fMassNucleus;}

  void SetMassCytoplasm(G4double mC) {fMassCytoplasm = mC;}
  G4double GetMassCytoplasm()  const {return fMassCytoplasm;}

  void AddDoseBox(G4int i, G4double x) {fDose3DDose[i] +=x;}
  G4double GetDoseBox(G4int i)   const {return fDose3DDose[i];}
  
  G4ThreeVector GetVectCell(G4int i) const {return fMapVoxels[i];}

private:

  const DetectorConstruction* fDetector;    

  G4int fSaveRndm;
  G4int fNumEvent;
  G4int fNbOfPixels;
  G4int fNbOfHitsGas;

  G4double fDoseN;
  G4double fDoseC;
  G4double fMassCytoplasm;
  G4double fMassNucleus;
  G4double * fDose3DDose;

  G4ThreeVector * fMapVoxels;

};

#endif
