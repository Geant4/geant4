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
//*******************************************************
//
// DicomGeometry.hh :
//	- Start the building of the geometry
//	- Creation of the world and other "mother"(middle) volume
//	- Initialisation of patient geometry
//	- Initialisation of HeaD geometry
// 	- Functions are in DicomGeometry.cc, PatientConstructor.cc
//*******************************************************

#ifndef DicomGeometry_h
#define DicomGeometry_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class DicomConfiguration;
class DicomPatientConstructor;
class G4Material;
class G4LogicalVolume;
class G4PhysicalVolume;
class G4Box;

class DicomGeometry : public G4VUserDetectorConstruction
{
public:

  DicomGeometry();
  ~DicomGeometry();

  G4VPhysicalVolume* Construct();

private:

  void InitialisationOfMaterials();
  void PatientConstruction();
 
  DicomPatientConstructor* patientConstructor;

  //Materials ...
 
  G4Material* trabecularBone; 
  G4Material* denseBone;
  G4Material* liver; 
  G4Material* muscle; 
  G4Material* phantom; 
  G4Material* breast; 
  G4Material* adiposeTissue; 
  G4Material* lungexhale; 
  G4Material* lunginhale;
  G4Material* air; 

  // World ...

  G4Box* solidWorld;
  G4LogicalVolume* logicWorld;
  G4VPhysicalVolume* physiWorld;
  G4VPhysicalVolume* parameterisedPhysVolume;  
  G4VPhysicalVolume*  physicalLungINhale;
};

#endif

