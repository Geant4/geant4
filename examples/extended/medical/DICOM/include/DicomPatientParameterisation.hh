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
// The code was written by :
//	*Louis Archambault louis.archambault@phy.ulaval.ca,
//      *Luc Beaulieu beaulieu@phy.ulaval.ca
//      +Vincent Hubert-Tremblay at tigre.2@sympatico.ca
//
//
// *Centre Hospitalier Universitaire de Quebec (CHUQ),
// Hotel-Dieu de Quebec, departement de Radio-oncologie
// 11 cote du palais. Quebec, QC, Canada, G1R 2J6
// tel (418) 525-4444 #6720
// fax (418) 691 5268
//
// + Université Laval, Québec (QC) Canada
//*******************************************************//
#ifndef DicomPatientParameterisation_h
#define DicomPatientParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include <vector>

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
			       G4Material* adipose,
			       G4Material* breast,
			       G4Material* phantom,
			       G4Material* muscle,
			       G4Material* liver,
			       G4Material* denseBone,
			       G4Material* trabecularBone);

  virtual ~DicomPatientParameterisation();

  virtual void ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const;

  //  virtual void ComputeDimensions (G4Box & voxels, const G4int copyNo, const G4VPhysicalVolume* physVol) const;

  virtual void ComputeDimension (G4Box&, 
				  const G4int, 
				  const G4VPhysicalVolume* ) const;

   virtual void ComputeDimension(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimension(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimension(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}
	
    virtual void ComputeDimension(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimension(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimension(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimension(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimension(G4Hype &,
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

  std::vector<G4double> density;
  std::vector<G4double> patientPlacementX;
  std::vector<G4double> patientPlacementY;
  std::vector<G4double> patientPlacementZ;


  //G4LogicalVolume* LogicalVolumeParam;

  G4double middleLocationValue;
};
#endif



