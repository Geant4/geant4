/***************************************************************************
                          DicomPatientParameterisation.h  -  description
                             -------------------
    begin                : Sat Feb 1 2003
    copyright            : (C) 2003 by Vincent Hubert-Tremblay
    email                : vihut@phy.ulaval.ca
 ***************************************************************************/
// library G4

#ifndef DicomPatientParameterisation_h
#define DicomPatientParameterisation_h 1


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4VPVParameterisation.hh"

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;

class G4VPhysicalVolume;
class G4Box;

class DicomPatientParameterisation : public G4VPVParameterisation
{
public:
  void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;
  void ComputeDimensions (G4Box & Voxels, const G4int copyNo, const G4VPhysicalVolume* physVol) const;

private:
  G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol);
  G4Material* P_lung_exhale;
  G4Material* P_lung_inhale;
  G4Material* P_adipose;
  G4Material* P_breast;
  G4Material* P_phantom;
  G4Material* P_muscle;
  G4Material* P_liver;
  G4Material* P_dense_bone;
  G4Material* P_trabecular_bone;

  G4double red;
  G4double green;
  G4double blue;
  G4double alpha;

  G4VisAttributes *Attributes_air;
  G4VisAttributes *Attributes_LungINhale;
  G4VisAttributes *Attributes_LungEXhale;
  G4VisAttributes *Attributes_Adipose;
  G4VisAttributes *Attributes_Breast;
  G4VisAttributes *Attributes_Phantom;
  G4VisAttributes *Attributes_Muscle;
  G4VisAttributes *Attributes_Liver;
  G4VisAttributes *Attributes_TrabecularBone;
  G4VisAttributes *Attributes_DenseBone;

public:
  DicomPatientParameterisation(G4int  NoVoxels , double max_density , double min_density ,
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

  void GetDensity( double maxdensity , double mindensity );
  char Densitybuf[300];
  FILE* readConfiguration; //lecturepref;
  int max;
  char name[300];
  int compression;
  char maxbuf[300];

  FILE* readData; //lectureDon;
  int columns,rows;
  double pixel_spacing_X,pixel_spacing_Y;
  double SliceTickness;
  double SliceLocation;
  char rowsbuf[300],columnsbuf[300];
  char pixel_spacing_Xbuf[300],pixel_spacing_Ybuf[300];
  char SliceTicknessbuf[300];
  char SliceLocationbuf[300];
  char compressionbuf[300];
  char fullname[300];

  G4double PatientX;
  G4double PatientY;
  G4double PatientZ;

  vector<double> Density;
  vector<double> PatientPlacementX;
  vector<double> PatientPlacementY;
  vector<double> PatientPlacementZ;


  G4LogicalVolume *LogicalVolumeParam;

  double MiddleLocationValue;

private:
  int lenc,lenr;
};
#endif



