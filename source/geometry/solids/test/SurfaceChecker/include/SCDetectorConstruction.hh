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
//
// $Id: SCDetectorConstruction.hh,v 1.6 2006-12-13 15:43:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SCDetectorConstruction_h
#define SCDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "SCMagneticField.hh"
#include "G4VSolid.hh"

class G4Box ;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class SCDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SCDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
     SCDetectorConstruction();
    ~SCDetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     const 
     G4VPhysicalVolume* GetTracker() {return physiTracker;};
     G4double GetWorldFullLength()   {return fWorldLength;}; 


     G4double GetTrackerR() { return fTrackerR ; } ;
     G4double GetTrackerR1() { return fTrackerR1 ; } ;
     G4double GetTrackerR2() { return fTrackerR2 ; } ;

     G4double GetTrackerpDz()  { return fTrackerpDz;};
     G4double GetTrackerpDx1() { return fTrackerpDx1 ; } ;
     G4double GetTrackerpDx2() { return fTrackerpDx2 ; } ;
     G4double GetTrackerpDx3() { return fTrackerpDx3 ; } ;
     G4double GetTrackerpDx4() { return fTrackerpDx4 ; } ;

     G4double GetTrackerpDy1() { return fTrackerpDy1 ; } ;
     G4double GetTrackerpDy2() { return fTrackerpDy2 ; } ;
     G4double GetTwistAngle() { return fTwistAngle ; } ;

     G4double GetAlpha() { return fAlph ; } ;
     G4double GetPhi() { return fPhi ; } ;
     G4double GetTheta() { return fTheta ; } ;
     G4String GetDetectorType() { return fval ; } ;

     G4double GetPhiSegment()   { return fPhiSegment ; } ;
     G4double GetThetaSegment() { return fThetaSegment ; } ;

     G4double GetSemiAxisX() { return fSemiAxisX ; }
     G4double GetSemiAxisY() { return fSemiAxisY ; }
     G4double GetSemiAxisZ() { return fSemiAxisZ ; }

     G4VSolid * GetSolid()  { return aVolume ; }

     void  SwitchDetector();
     G4VPhysicalVolume* SelectDetector (const G4String& val);

     void setTargetMaterial (G4String);
     void setChamberMaterial(G4String);
     void SetMagField(G4double);
     
  private:

     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
     G4VSolid* aVolume;

  //     G4TwistedTrap*     solidTracker;

     G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker
     
     SCMagneticField* fpMagField;   // pointer to the magnetic field 
     
     SCDetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume


     G4double fTrackerR ;   // a radius
     G4double fTrackerR1 ;   // r1
     G4double fTrackerR2 ;   // r2
    
     G4double fTrackerpDz  ;          // Full length of Tracker (pDz)
     G4double fTrackerpDx1 ;  // twisted Trapezoid
     G4double fTrackerpDx2 ;
     G4double fTrackerpDx3 ;
     G4double fTrackerpDx4 ;
     G4double fTrackerpDy1 ;
     G4double fTrackerpDy2 ;
     G4double fTwistAngle ;
     G4double fPhi ;
     G4double fTheta ;
     G4double fAlph ;
     G4double fSemiAxisX ;
     G4double fSemiAxisY ;
     G4double fSemiAxisZ ;


     G4double fPhiSegment ;    // e.g. for a sphere : phi in [phi,phi+segment]
     G4double fThetaSegment ;
 
     G4String fval ;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
