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
#ifndef LXeMainVolume_H
#define LXeMainVolume_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "LXeDetectorConstruction.hh"

class LXeWLSFiber : public G4PVPlacement
{
public:
  LXeWLSFiber(G4RotationMatrix *pRot,
		const G4ThreeVector &tlate,
		G4LogicalVolume *pMotherLogical,
		G4bool pMany,
		G4int pCopyNo,
		LXeDetectorConstruction* c);
private:

  void CopyValues();

  static G4LogicalVolume* clad2_log;

  G4bool updated; //does the fiber need to be rebuilt
  
  G4double fiber_rmin;    
  G4double fiber_rmax;    
  G4double fiber_z;
  G4double fiber_sphi;
  G4double fiber_ephi;

  G4double clad1_rmin;
  G4double clad1_rmax;    
  G4double clad1_z;
  G4double clad1_sphi;
  G4double clad1_ephi; 
  
  G4double clad2_rmin;
  G4double clad2_rmax;    
  G4double clad2_z;
  G4double clad2_sphi;
  G4double clad2_ephi;

  LXeDetectorConstruction* constructor;
};

#endif
