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

  DicomPatientParameterisation(G4int nVoxels, 
			       G4double maxDensity, 
			       G4double minDensity);

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

  G4int compression;
  
  G4int columns,rows;
  G4double pixelSpacingX;
  G4double pixelSpacingY;
  G4double sliceThickness;
  G4double sliceLocation;

  G4std::vector<G4double> density;
  G4std::vector<G4double> patientPlacementX;
  G4std::vector<G4double> patientPlacementY;
  G4std::vector<G4double> patientPlacementZ;


  //G4LogicalVolume* LogicalVolumeParam;

  G4double middleLocationValue;
};
#endif



