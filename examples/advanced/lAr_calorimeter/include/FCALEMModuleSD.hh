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
// $Id: FCALEMModuleSD.hh,v 1.5 2002-12-17 15:53:21 pmendez Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALEMModuleSD_h
#define FCALEMModuleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "FCALCalorHit.hh"

class FCALTestbeamSetup;
class FCALEMModule;
class G4HCofThisEvent;
class G4Step;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALEMModuleSD : public G4VSensitiveDetector
{
  public:
  
  //      FCALEMModuleSD(G4String, FCALTestbeamSetup* );
      FCALEMModuleSD(G4String); 
     ~FCALEMModuleSD();

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
  FCALEMModule* EmModule;
  G4int Init_state;
  
  public:

  G4double EvisF1Tile[1131];

};

#endif

