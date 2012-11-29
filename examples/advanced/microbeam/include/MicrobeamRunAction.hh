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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#ifndef MicrobeamRunAction_h
#define MicrobeamRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"

#include "MicrobeamDetectorConstruction.hh"
#include "MicrobeamHistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class MicrobeamRunAction : public G4UserRunAction
{
public:
  
  MicrobeamRunAction(MicrobeamDetectorConstruction*, MicrobeamHistoManager*);
  ~MicrobeamRunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
  void  SetRndmFreq(G4int   val)  {saveRndm = val;}
  G4int GetRndmFreq()             {return saveRndm;}

  void AddDoseN(G4float dose){ DoseN += dose;}
  void SetDoseN(G4float dose){ DoseN = dose;}
  G4float GetDoseN(){return DoseN;}

  void AddDoseC(G4float dose){ DoseC += dose;}
  void SetDoseC(G4float dose){ DoseC = dose;}
  G4float GetDoseC(){return DoseC;}

  G4int GetNumEvent(){return numEvent;}
  void SetNumEvent(G4int i){numEvent = i;}

  G4int GetNbOfHitsGas(){return nbOfHitsGas;}
  void AddNbOfHitsGas(){nbOfHitsGas = nbOfHitsGas+1;}

  void SetMassNucleus(G4float mN){ massNucleus = mN;}
  G4float GetMassNucleus(){return massNucleus;}

  void SetMassCytoplasm(G4float mC){ massCytoplasm = mC;}
  G4float GetMassCytoplasm(){return massCytoplasm;}

  void AddDoseBox(G4int i, G4float x){ dose3DDose[i] +=x;}
  G4float GetDoseBox(G4int i){ return dose3DDose[i];}
  
  G4ThreeVector GetVectCell(G4int i) {return mapVoxels[i];}

private:

  MicrobeamDetectorConstruction* Detector;    
  MicrobeamHistoManager* Histo;
  MicrobeamPhantomConfiguration myMicrobeamPhantomConfiguration;  

  G4int saveRndm;
  G4int numEvent;
  G4int nbOfPixels;
  G4int nbOfHitsGas;
  G4float SP;
  G4float R;
  G4float RnElec;
  G4float RcElec;
  G4float DoseN;
  G4float DoseC;
  G4float DoseNElec;
  G4float DoseCElec;
  G4float massPhantom;
  G4float massCytoplasm;
  G4float massNucleus;
  G4bool boolSP;
  
  G4int * x3DDose;
  G4int * y3DDose;
  G4int * z3DDose;
  G4float * dose3DDose;
  G4ThreeVector * mapVoxels;

};

#endif
