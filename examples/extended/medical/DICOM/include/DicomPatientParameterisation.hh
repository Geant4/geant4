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

  G4double red;
  G4double green;
  G4double blue;
  G4double alpha;

  G4VisAttributes* Attributes_air;
  G4VisAttributes* Attributes_LungINhale;
  G4VisAttributes* Attributes_LungEXhale;
  G4VisAttributes* Attributes_Adipose;
  G4VisAttributes* Attributes_Breast;
  G4VisAttributes* Attributes_Phantom;
  G4VisAttributes* Attributes_Muscle;
  G4VisAttributes* Attributes_Liver;
  G4VisAttributes* Attributes_TrabecularBone;
  G4VisAttributes* Attributes_DenseBone;

  char Densitybuf[300];
  FILE* readConfiguration; //lecturepref;
  G4int max;
  char name[300];
  G4int compression;
  char maxbuf[300];

  FILE* readData; //lectureDon;
  G4int columns,rows;
  G4double pixel_spacing_X,pixel_spacing_Y;
  G4double SliceTickness;
  G4double SliceLocation;
  char rowsbuf[300],columnsbuf[300];
  char pixel_spacing_Xbuf[300],pixel_spacing_Ybuf[300];
  char SliceTicknessbuf[300];
  char SliceLocationbuf[300];
  char compressionbuf[300];
  char fullname[300];

  G4double PatientX;
  G4double PatientY;
  G4double PatientZ;

  G4std::vector<double> Density;
  G4std::vector<double> PatientPlacementX;
  G4std::vector<double> PatientPlacementY;
  G4std::vector<double> PatientPlacementZ;


  G4LogicalVolume* LogicalVolumeParam;

  G4double MiddleLocationValue;

  G4int lenc,lenr;
};
#endif



