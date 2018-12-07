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
// This is the *BASIC* version of IORT, a Geant4-based application
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////


#ifndef IORTDetectorConstruction_H
#define IORTDetectorConstruction_H 1

#include "G4Box.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4Tubs.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class IORTDetectorMessenger;

class IORTDetectorConstruction 
{
public:

  IORTDetectorConstruction(G4VPhysicalVolume*);

  ~IORTDetectorConstruction();

  //  G4VPhysicalVolume *detectorPhysicalVolume;  aggiunto

private: 

  void ConstructPhantom();
  void ConstructDetector();
  //  void ConstructDisc();
  void ConstructSensitiveDetector();
  void ParametersCheck();

public: 
// Get detector position relative to WORLD
inline G4ThreeVector GetDetectorToWorldPosition()
  {
    return phantomPosition + detectorPosition;
  }

/////////////////////////////////////////////////////////////////////////////
// Get displacement between phantom and detector by detector position (center of), phantom (center of) and detector sizes
inline G4ThreeVector GetDetectorToPhantomPosition()
{
    return G4ThreeVector(phantomSizeX/2 - detectorSizeX/2 + detectorPosition.getX(),
                         phantomSizeY/2 - detectorSizeY/2 + detectorPosition.getY(),
                         phantomSizeZ/2 - detectorSizeZ/2 + detectorPosition.getZ()
		          );
}

/////////////////////////////////////////////////////////////////////////////
// Calculate (and set) detector position by displacement, phantom and detector sizes
inline void SetDetectorPosition()
  {
	  // Adjust detector position
	  detectorPosition.setX(detectorToPhantomPosition.getX() - phantomSizeX/2 + detectorSizeX/2);
	  detectorPosition.setY(detectorToPhantomPosition.getY() - phantomSizeY/2 + detectorSizeY/2);
	  detectorPosition.setZ(detectorToPhantomPosition.getZ() - phantomSizeZ/2 + detectorSizeZ/2);
     
    //G4cout << "*************** DetectorToPhantomPosition " << detectorToPhantomPosition/cm << "\n";
    //G4cout << "*************** DetectorPosition " << detectorPosition/cm << "\n";
  }
/////////////////////////////////////////////////////////////////////////////
// Check whether detector is inside phantom
inline bool IsInside(G4double detectorX,
		     G4double detectorY,
		     G4double detectorZ,
		     G4double phantomX,
		     G4double phantomY,
		     G4double phantomZ,
		     G4ThreeVector detToPhantomPosition)
{
// Dimensions check... X Y and Z
// Firstly check what dimension we are modifying
	{
	    if (detectorX > phantomX) 
		 {
		    G4cout << "Error: Detector X dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomX - detectorX) < detToPhantomPosition.getX()) 
	         {
		    G4cout << "Error: X dimension doesn't fit with detector to phantom relative position" << G4endl;
		    return false;
	         }
	}

	{
	    if (detectorY > phantomY) 
		 {
		    G4cout << "Error: Detector Y dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomY - detectorY) < detToPhantomPosition.getY()) 
	     {
		   G4cout << "Error: Y dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}			 

	{
	    if (detectorZ > phantomZ) 
		 {
		    G4cout << "Error: Detector Z dimension must be smaller or equal to the correspondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( (phantomZ - detectorZ) < detToPhantomPosition.getZ()) 
	     {
		   G4cout << "Error: Z dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}

	return true;
}
/////////////////////////////////////////////////////////////////////////////

  G4bool  SetPhantomMaterial(G4String material);
  void SetVoxelSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  void SetPhantomPosition(G4ThreeVector);
  void SetDetectorToPhantomPosition(G4ThreeVector DetectorToPhantomPosition);
  void UpdateGeometry();
  void DeleteDisc();
  void ConstructDisc();
  void PrintParameters();
  G4LogicalVolume* GetDetectorLogicalVolume(){ return detectorLogicalVolume;}

  G4bool  SetDiscoMaterialIORT(G4String material);
  void SetOuterRadiusDiscoIORT(G4double outerr);
  void SetinnerRadiusDiscoIORT(G4double innerr);
  void SetheightDiscoIORT(G4double height);
  void SetDiscoXPositionIORT(G4double xpos);
  void SetDiscoYPositionIORT(G4double ypos);
  void SetDiscoZPositionIORT(G4double zpos);

  G4bool  SetDiscoMaterialIORT1(G4String material);
  void SetOuterRadiusDiscoIORT1(G4double outerr);
  void SetinnerRadiusDiscoIORT1(G4double innerr);
  void SetheightDiscoIORT1(G4double height);
  void SetDiscoXPositionIORT1(G4double xpos);

  void SetAngleDiscoIORT0(G4double phi0);

private:

  IORTDetectorMessenger* detectorMessenger; 

  G4VisAttributes* red;

  G4VPhysicalVolume* motherPhys;

  G4Box *phantom , *detector;
  G4LogicalVolume *phantomLogicalVolume, *detectorLogicalVolume; 
  G4VPhysicalVolume *phantomPhysicalVolume,   *detectorPhysicalVolume;
  
  G4double phantomSizeX; 
  G4double phantomSizeY; 
  G4double phantomSizeZ;

  G4double detectorSizeX; 
  G4double detectorSizeY; 
  G4double detectorSizeZ;

  G4ThreeVector phantomPosition, detectorPosition, detectorToPhantomPosition; //  phantom center, detector center, detector to phantom relative position

  G4Material *phantomMaterial, *detectorMaterial;
  G4Region* aRegion;
  
  //Disco0 IORT
  G4Tubs* solidDiscoIORT0;
  G4LogicalVolume* logicDiscoIORT0;
  G4VPhysicalVolume* physiDiscoIORT0;
  G4double AngleDiscoIORT0; 

  // Disco1 IORT
  G4VisAttributes* white;
  G4VisAttributes* gray;
  G4VisAttributes* gray1;
  G4double innerRadiusDiscoIORT;
  G4double OuterRadiusDiscoIORT;
  G4double heightDiscoIORT;
  G4double DiscoXPositionIORT;
  G4double DiscoYPositionIORT;
  G4double DiscoZPositionIORT;
  G4Tubs* solidDiscoIORT; 
  G4LogicalVolume* logicDiscoIORT;
  G4VPhysicalVolume* physiDiscoIORT;
  G4Material* DiscoMaterialIORT;

   // Disco2 IORT
  
  G4double innerRadiusDiscoIORT1;
  G4double OuterRadiusDiscoIORT1;
  G4double heightDiscoIORT1;
  G4double DiscoXPositionIORT1;
  G4Tubs* solidDiscoIORT1; 
  G4LogicalVolume* logicDiscoIORT1;
  G4VPhysicalVolume* physiDiscoIORT1;
  G4Material* DiscoMaterialIORT1;
};
#endif



