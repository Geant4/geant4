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
// $Id: HadrontherapyDetectorConstruction.hh; 
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyDetectorConstruction_H
#define HadrontherapyDetectorConstruction_H 1

#include "G4Box.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class HadrontherapyDetectorROGeometry;
class HadrontherapyDetectorMessenger;
class HadrontherapyDetectorSD;
class HadrontherapyMatrix;

class HadrontherapyDetectorConstruction 
{
public:

  HadrontherapyDetectorConstruction(G4VPhysicalVolume*);

  ~HadrontherapyDetectorConstruction();


private: 

  void ConstructPhantom();
  void ConstructDetector();
  void ConstructSensitiveDetector(G4ThreeVector position_respect_to_WORLD);
  
public: 
// Get detector position relative to WORLD
inline G4ThreeVector GetDetectorToWorldPosition()
  {
    return phantomPosition + detectorPosition;
  }
/////////////////////////////////////////////////////////////////////////////
// Get displacement between phantom and detector by detector position, phantom and detector sizes
inline G4ThreeVector GetDetectorToPhantomPosition()
{
    return G4ThreeVector(phantomSizeX - detectorSizeX + detectorPosition.getX(),
                         phantomSizeY - detectorSizeY + detectorPosition.getY(),
                         phantomSizeZ - detectorSizeZ + detectorPosition.getZ()
		          );
}

/////////////////////////////////////////////////////////////////////////////
// Calculate (and set) detector position by displacement, phantom and detector sizes
inline void SetDetectorPosition()
  {
	  // Adjust detector position
	  detectorPosition.setX(detectorToPhantomPosition.getX() - phantomSizeX + detectorSizeX);
	  detectorPosition.setY(detectorToPhantomPosition.getY() - phantomSizeY + detectorSizeY);
	  detectorPosition.setZ(detectorToPhantomPosition.getZ() - phantomSizeZ + detectorSizeZ);
     
      if (detectorPhysicalVolume) detectorPhysicalVolume -> SetTranslation(detectorPosition); 
  }
/////////////////////////////////////////////////////////////////////////////
// Check whether detector is inside phantom
inline bool IsInside(G4double detectorHalfX,
		     G4double detectorHalfY,
		     G4double detectorHalfZ,
		     G4double phantomHalfX,
		     G4double phantomHalfY,
		     G4double phantomHalfZ,
		     G4ThreeVector detectorToPhantomPosition)
{
// Dimensions check... X Y and Z
// Firstly check what dimension we are modifying
	if (detectorHalfX > 0. && phantomHalfX > 0. && detectorToPhantomPosition.getX() >=0.)
	{
	    if (detectorHalfX > phantomHalfX) 
		 {
		    G4cout << "Error: Detector X dimension must be smaller or equal to the corrispondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( 2*(phantomHalfX - detectorHalfX) < detectorToPhantomPosition.getX()) 
	         {
		    G4cout << "Error: X dimension doesn't fit with detector to phantom relative position" << G4endl;
		    return false;
	         }
	}

	if (detectorHalfY > 0. && phantomHalfY > 0.&& detectorToPhantomPosition.getY() >=0.)
	{
	    if (detectorHalfY > phantomHalfY) 
		 {
		    G4cout << "Error: Detector Y dimension must be smaller or equal to the corrispondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( 2*(phantomHalfY - detectorHalfY) < detectorToPhantomPosition.getY()) 
	     {
		   G4cout << "Error: Y dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}			 

	if (detectorHalfZ > 0. && phantomHalfZ > 0.&& detectorToPhantomPosition.getZ() >=0.)
	{
	    if (detectorHalfZ > phantomHalfZ) 
		 {
		    G4cout << "Error: Detector Z dimension must be smaller or equal to the corrispondent of the phantom" << G4endl;
		    return false;
		 }
	    if ( 2*(phantomHalfZ - detectorHalfZ) < detectorToPhantomPosition.getZ()) 
	     {
		   G4cout << "Error: Z dimension doesn't fit with detector to phantom relative position" << G4endl;
		   return false;
	     }
	}

    G4cout << "Displacement between Phantom and Detector is: "; 
    G4cout << "DX= "<< G4BestUnit(detectorToPhantomPosition.getX(),"Length") << 
              "DY= "<< G4BestUnit(detectorToPhantomPosition.getY(),"Length") << 
              "DZ= "<< G4BestUnit(detectorToPhantomPosition.getZ(),"Length") << G4endl;

	return true;
}
/////////////////////////////////////////////////////////////////////////////

  G4bool SetNumberOfVoxelBySize(G4double sizeX, G4double sizeY, G4double sizeZ);
  G4bool SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  G4bool SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ);
  G4bool SetPhantomPosition(G4ThreeVector);
  G4bool SetDetectorToPhantomPosition(G4ThreeVector DetectorToPhantomPosition);
  G4LogicalVolume* GetDetectorLogicalVolume(){ return detectorLogicalVolume;}


private:

  HadrontherapyDetectorMessenger* detectorMessenger; 

  G4VisAttributes* skyBlue;
  G4VisAttributes* red;

  G4VPhysicalVolume* motherPhys;

  HadrontherapyDetectorSD*         detectorSD; // Pointer to sensitive detector
  HadrontherapyDetectorROGeometry* detectorROGeometry; // Pointer to ROGeometry 
  HadrontherapyMatrix*             matrix;

  G4VPhysicalVolume* phantomPhysicalVolume;
  G4LogicalVolume*   phantomLogicalVolume; 
  G4LogicalVolume*   detectorLogicalVolume;
  G4VPhysicalVolume* detectorPhysicalVolume;
  
  G4double phantomSizeX; 
  G4double phantomSizeY; 
  G4double phantomSizeZ;

  G4double detectorSizeX; 
  G4double detectorSizeY; 
  G4double detectorSizeZ;

  G4ThreeVector phantomPosition, detectorPosition, detectorToPhantomPosition; //  phantom center, detector center, detector to phantom relative position

  G4double sizeOfVoxelAlongX; 
  G4double sizeOfVoxelAlongY; 
  G4double sizeOfVoxelAlongZ; 

  G4int numberOfVoxelsAlongX; 
  G4int numberOfVoxelsAlongY;
  G4int numberOfVoxelsAlongZ;  

  G4Box* phantom;
  G4Box* detector;

};
#endif
