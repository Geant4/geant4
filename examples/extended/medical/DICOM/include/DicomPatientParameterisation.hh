#ifndef DicomPatientParameterisation_h
#define DicomPatientParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "g4std/vector"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Box;
class G4Material;
class G4VisAttributes;
class DicomConfiguration;
class DicomPatientParameterisation : public G4VPVParameterisation
{
public:

  DicomPatientParameterisation(G4int NoVoxels, 
			       G4double maxDensity, 
			       G4double minDensity ,
			       G4Material* lunginhale,
			       G4Material* lungexhale,
			       G4Material* adipose_tissue,
			       G4Material* breast,
			       G4Material* phantom,
			       G4Material* muscle,
			       G4Material* liver,
			       G4Material* dense_bone,
			       G4Material* trabecular_bone);

  virtual ~DicomPatientParameterisation();

  virtual void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;

  //  virtual void ComputeDimensions (G4Box & voxels, const G4int copyNo, const G4VPhysicalVolume* physVol) const;

  virtual void ComputeDimensions (G4Box&, 
				  const G4int, 
				  const G4VPhysicalVolume* ) const;

   virtual void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol);

  void GetDensity( G4double maxDensity, G4double minDensity );

private:
 
  //materials ...
  G4Material* lungExhale;
  G4Material* lungInhale;
  G4Material* adiposeTissue;
  G4Material* breastTissue;
  G4Material* phantomTissue;
  G4Material* muscleTissue;
  G4Material* liverTissue;
  G4Material* denseBoneTissue;
  G4Material* trabecularBoneTissue;

  G4VisAttributes* attributeAir;
  G4VisAttributes* attributeLungINhale;
  G4VisAttributes* attributeLungEXhale;
  G4VisAttributes* attributeAdipose;
  G4VisAttributes* attributeBreast;
  G4VisAttributes* attributePhantom;
  G4VisAttributes* attributeMuscle;
  G4VisAttributes* attributeLiver;
  G4VisAttributes* attributeTrabecularBone;
  G4VisAttributes* attributeDenseBone;

  G4int max;
  G4int compression;
  
  FILE* readData; 
  G4int columns,rows;
  G4double pixelSpacingX;
  G4double pixelSpacingY;
  G4double sliceThickness;
  G4double sliceLocation;
 
  //G4double PatientX;
  //G4double PatientY;
  //G4double PatientZ;

  G4std::vector<G4double> Density;
  G4std::vector<G4double> PatientPlacementX;
  G4std::vector<G4double> PatientPlacementY;
  G4std::vector<G4double> PatientPlacementZ;


  G4LogicalVolume* LogicalVolumeParam;

  G4double MiddleLocationValue;
};
#endif



