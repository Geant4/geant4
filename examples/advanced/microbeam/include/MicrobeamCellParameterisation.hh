// -------------------------------------------------------------------
// $Id: MicrobeamCellParameterisation.hh,v 1.1 2006-04-06 15:32:43 sincerti Exp $
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
     G4Material* nucleus, G4Material* cytoplasm);

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

    G4Material * nucleusMaterial;
    G4Material * cytoplasmMaterial ;
    
    G4VisAttributes * nucleusAttributes;
    G4VisAttributes * cytoplasmAttributes;
    
    // PHANTOM SIZE IN VOXELS
    G4ThreeVector *mapCell ; // VOXEL COORDINATES
    G4float *material      ; // MATERIAL   
};

#endif


