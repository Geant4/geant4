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
// $Id: SCDetectorConstruction.hh,v 1.1 2005-05-19 13:07:29 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SCDetectorConstruction_h
#define SCDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "SCMagneticField.hh"

class G4Box;
class G4Tubs;
class G4TwistedTubs;
class G4TwistedTrap;
class G4Sphere;
class G4Torus;
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

     void setTargetMaterial (G4String);
     void setChamberMaterial(G4String);
     void SetMagField(G4double);
     
  private:

     G4Box*             solidWorld;    // pointer to the solid envelope 
     G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
     G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

  //     G4TwistedTrap*     solidTracker;
     G4Torus*             solidTracker ;
     G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
     G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker
     
     SCMagneticField* fpMagField;   // pointer to the magnetic field 
     
     SCDetectorMessenger* detectorMessenger;  // pointer to the Messenger
       
     G4double fWorldLength;            // Full length of the world volume

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
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
