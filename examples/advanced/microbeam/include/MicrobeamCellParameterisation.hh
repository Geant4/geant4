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
// -------------------------------------------------------------------
// $Id$
// -------------------------------------------------------------------

#ifndef MicrobeamCellParameterisation_H
#define MicrobeamCellParameterisation_H 1

#include "G4VPVParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class Polycone;
class Polyhedra;
class G4Hype;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MicrobeamCellParameterisation : public G4VPVParameterisation
{ 
  public:
  
    MicrobeamCellParameterisation
    (G4int NoBoxes, G4float DimBoxX, G4float DimBoxY, G4float DimBoxZ,
     G4Material* nucleus1, G4Material* cytoplasm1,
     G4Material* nucleus2, G4Material* cytoplasm2,
     G4Material* nucleus3, G4Material* cytoplasm3
     );

    virtual				 
   ~MicrobeamCellParameterisation();
   
    void ComputeTransformation (const G4int copyNo,G4VPhysicalVolume* physVol) const;
    
    void ComputeDimensions(G4Box&, 
				  const G4int, 
				  const G4VPhysicalVolume* ) const;

    void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    G4int GetNoBoxes() {return NoCellBoxes;}

    G4Material* ComputeMaterial(const G4int copyNo,
                                      G4VPhysicalVolume* physVol,
                                const G4VTouchable*);

  private:

    G4int   NoCellBoxes;
    G4float DimCellBoxX;
    G4float DimCellBoxY;
    G4float DimCellBoxZ;

    G4Material * nucleusMaterial1;
    G4Material * cytoplasmMaterial1;
    G4Material * nucleusMaterial2;
    G4Material * cytoplasmMaterial2;
    G4Material * nucleusMaterial3;
    G4Material * cytoplasmMaterial3;
    
    G4VisAttributes * nucleusAttributes1;
    G4VisAttributes * cytoplasmAttributes1;
    G4VisAttributes * nucleusAttributes2;
    G4VisAttributes * cytoplasmAttributes2;
    G4VisAttributes * nucleusAttributes3;
    G4VisAttributes * cytoplasmAttributes3;
    
    // PHANTOM SIZE IN VOXELS
    G4ThreeVector *mapCell ; // VOXEL COORDINATES
    G4float *material      ; // MATERIAL 
    G4float *mass          ; // DENSITY REGION
      
};

#endif


