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
// $Id: Pers01CalorimeterSD.hh,v 1.1 2002-12-04 02:44:25 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Pers01CalorimeterSD_h
#define Pers01CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Pers01DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Pers01CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Pers01CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Pers01CalorimeterSD(G4String, Pers01DetectorConstruction* );
     ~Pers01CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      Pers01CalorHitsCollection*  CalCollection;      
      Pers01DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

