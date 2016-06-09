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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraPMTSD.hh
//   **********************************************
//
//    Class used to define the Ultra photomultiplier as a sensitive detector.
//    Hits in this sensitive detector are defined in the UltraOpticalHit class
//
#ifndef UltraPMTSD_h
#define UltraPMTSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Material;
class G4HCofThisEvent;
class G4Step;

#include "UltraOpticalHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class UltraPMTSD : public G4VSensitiveDetector
{
public:
  
  UltraPMTSD(G4String);
  ~UltraPMTSD();
  
  void    Initialize(G4HCofThisEvent*);
  G4bool  ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void    EndOfEvent(G4HCofThisEvent*);

  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  UltraOpticalHitsCollection*  OpticalHitsCollection;   


};

#endif
