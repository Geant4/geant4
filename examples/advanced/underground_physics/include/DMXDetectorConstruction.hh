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
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DMXDetectorConstruction_h
#define DMXDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class G4UserLimits;

class DMXScintSD;
class DMXPmtSD;

class DMXDetectorConstruction : public G4VUserDetectorConstruction 
{
public:

  DMXDetectorConstruction();
  ~DMXDetectorConstruction();

public:

  G4VPhysicalVolume* Construct();

  // methods for UserLimits
  void      UseUserLimits(G4bool value); 
  G4bool    IsUseUserLimits()       { return fUseUserLimits; } 
  G4double  GetMaxTime() const      { return theMaxTimeCuts;  }
  G4double  GetMaxStepSize() const  { return theMaxStepSize;  }

  G4double  GetRoomTime() const      { return theRoomTimeCut;  }

  void  SetMaxTime(G4double value);  
  void  SetMaxStepSize(G4double value);  
  void  SetRoomTime(G4double value);

private:

  void DefineMaterials();

  G4bool           fUseUserLimits;
  G4UserLimits*    theUserLimits; 
  G4double         theMaxTimeCuts;
  G4double         theMaxStepSize;
  
  G4double         theRoomTimeCut;

  G4Material*     world_mat;            // materials used
  G4Material*       lab_mat;        
  G4Material*    jacket_mat;
  G4Material*    vacuum_mat;
  G4Material*    vessel_mat;
  G4Material*  detector_mat;   
  G4Material*  liqPhase_mat;
  G4Material*  CuShield_mat;
  G4Material*     alpha_mat;
  G4Material*    recess_mat;
  G4Material* americium_mat;
  G4Material*      ring_mat;
  G4Material*    mirror_mat;
  G4Material*      grid_mat;
  G4Material*       pmt_mat;
  G4Material*    window_mat;
  G4Material*    phcath_mat;


  G4double worldRadius;                // sizes
  G4double worldHeight;
  G4double sourceZ;
  G4double sourceZ2;


  G4LogicalVolume*   world_log;        // pointers
  G4VPhysicalVolume* world_phys;  
  G4LogicalVolume*   lab_log;
  G4VPhysicalVolume* lab_phys;  
  G4LogicalVolume*   jacket_log;
  G4VPhysicalVolume* jacket_phys;
  G4LogicalVolume*   vacuum_log;
  G4VPhysicalVolume* vacuum_phys;
  G4LogicalVolume*   vessel_log;
  G4VPhysicalVolume* vessel_phys;
  G4LogicalVolume*   detector_log;
  G4VPhysicalVolume* detector_phys;  
  G4LogicalVolume*   CuShield_log; 
  G4VPhysicalVolume* CuShield_phys;  
  G4LogicalVolume*   liqPhase_log; 
  G4VPhysicalVolume* liqPhase_phys;  
  G4LogicalVolume*   alpha_log;   
  G4VPhysicalVolume* alpha_phys; 
  G4LogicalVolume*   recess_log;   
  G4VPhysicalVolume* recess_phys; 
  G4LogicalVolume*   americium_log;   
  G4VPhysicalVolume* americium_phys; 
  G4LogicalVolume*   ring_log;   
  G4VPhysicalVolume* ring_phys_gas[1]; 
  G4VPhysicalVolume* ring_phys_liq[5]; 
  //  G4VPhysicalVolume* ring_phys[8]; 
  G4LogicalVolume*   mirror_log;   
  G4VPhysicalVolume* mirror_phys; 
  G4LogicalVolume*   grid1_log;   
  G4VPhysicalVolume* grid1_phys; 
  G4LogicalVolume*   grid2_log;   
  G4VPhysicalVolume* grid2_phys; 
  G4LogicalVolume*   pmt_log;   
  G4VPhysicalVolume* pmt_phys; 
  /*
  G4LogicalVolume*   window_log;   
  G4VPhysicalVolume* window_phys; 
  */
  G4LogicalVolume*   phcath_log;
  G4VPhysicalVolume* phcath_phys; 


  DMXScintSD*  LXeSD;            //pointer to sensitive detectors
  DMXPmtSD* pmtSD;

};

#endif

