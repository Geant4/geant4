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
// $Id: Tst50DetectorConstruction.hh,v 1.2 2002-11-29 11:19:29 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst50DetectorConstruction_h
#define Tst50DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class Tst50DetectorMessenger;
class  Tst50TrackerSD;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst50DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     Tst50DetectorConstruction();
    ~Tst50DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();

     const 
     G4VPhysicalVolume* GetTracker() {return physiTracker;};
  /*     G4double GetTrackerFullLength() {return fTrackerLength;};
     G4double GetTargetFullLength()  {return fTargetLength;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 
  */
     void setTargetMaterial (G4String);
public:
 void      UseUserLimits(G4bool value); 
  void  SetMaxStep(G4double value); 
 
private:
 G4bool           fUseUserLimits;
 G4UserLimits*    theUserLimits; 
 G4double         theMaxStep;
    
  private:

     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
     
     G4Box*             solidTarget;   // pointer to the solid Target
     G4LogicalVolume*   logicTarget;   // pointer to the logical Target
     G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
               
     G4Box*             solidTracker;  // pointer to the solid Tracker
     G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker
     
     G4Material*         TargetMater;  // pointer to the target  material
     G4Material*         ChamberMater; // pointer to the chamber material     

  Tst50TrackerSD* pTargetSD;
     Tst50DetectorMessenger* detectorMessenger;  // pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
