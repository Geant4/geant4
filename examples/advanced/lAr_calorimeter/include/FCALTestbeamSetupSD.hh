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
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALTestbeamSetupSD.hh,v 1.3 2002-12-12 19:16:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALTestbeamSetupSD_h
#define FCALTestbeamSetupSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class FCALTestbeamSetup;
class G4HCofThisEvent;
class G4Step;
#include "FCALCalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALTestbeamSetupSD : public G4VSensitiveDetector
{
  public:
  
  //      FCALTestbeamSetupSD(G4String, FCALTestbeamSetup* );
      FCALTestbeamSetupSD(G4String); 
     ~FCALTestbeamSetupSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
  //  FCALCalorHitsCollection*  CalCollection;      
  //      FCALTestbeamSetup* Detector;
  //   G4int*                   HitID;
  
  G4int InitBeam;

public:
  
  G4double EBeamS1, EBeamS2, EBeamS3;
  G4double EHoleScint, EBeamHole;
  G4double EBeamDead;
  G4int TailCatcherID;
  G4double ETailVis[8], ETailDep[7];
  

};

#endif

