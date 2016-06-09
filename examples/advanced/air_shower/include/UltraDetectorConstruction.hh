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
//   *        UltraDetectorConstruction.hh
//   **********************************************
//
//    Class used in the definition of the Ultra setup consisting of:
//      - the UVscope detector
//      - an optional reflecting surface
//    Optical photons can reach the UVscope either directly or after reflection in the
//    surface, which can be polished or diffusing.
//    The main part of the UVscope definition is the Fresnel lens construction based
//    on the UltraFresnelLens class.
//
#ifndef UltraDetectorConstruction_h
#define UltraDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4SDManager;
class UltraScintSD;
class UltraPMTSD;
class UltraFresnelLens;

class UltraDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    UltraDetectorConstruction();
    ~UltraDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();

  private:
    
  // Methods to build ULTRA
    G4VPhysicalVolume  *ConstructUVscope(G4VPhysicalVolume *);
    G4VPhysicalVolume  *ConstructMirror(G4VPhysicalVolume *);
    G4VPhysicalVolume  *ConstructGround(G4VPhysicalVolume *);
    UltraFresnelLens   *FresnelLens ;

    // Material definitions
    void ConstructTableMaterials();


  private:
    UltraPMTSD*   PMTSD  ;          //pointer to the photomultiplier sensitive detector
    G4SDManager*  SDmanager ;       // Sensitive Detector Manager

};

#endif 
