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
//
// $Id: BrachyDetectorConstruction.hh,v 1.14 2003-05-09 16:52:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//  
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************
// this class manages the geometry of the simulation set up
//

#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class BrachyPhantomSD;
class BrachyDetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Tubs;
class G4Box;
class G4Sphere;
class G4Tubs;
class G4Colour;
class G4VPhysicalVolume;
class BrachyPhantomSD;
class BrachyPhantomROGeometry;
class G4VPhysicalVolume;
class BrachyMaterial;
class BrachyFactory;
class BrachyVoxelParameterisation;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  BrachyDetectorConstruction(G4String&);
  ~BrachyDetectorConstruction();

  G4VPhysicalVolume*   Construct();  

  void PrintDetectorParameters(); 
  void SetAbsorberMaterial(G4String);

  const G4double VoxelWidth_X() {return m_BoxDimX/NumVoxelX;}
  const G4double VoxelWidth_Z() {return m_BoxDimZ/NumVoxelZ;}
  const G4int   GetNumVoxelX()  {return NumVoxelX;}
  const G4int   GetNumVoxelZ()  {return NumVoxelZ;}
  const G4double GetDimX()      {return m_BoxDimX;}
  const G4double GetBoxDim_Z()  {return  m_BoxDimZ;}

  void SelectDetector(G4String val);
  void SwitchDetector();

  void ConstructPhantom();
  void ConstructSensitiveDetector();
  void ComputeDimVoxel() { dimVoxel=m_BoxDimX/NumVoxelX; }

private:

  G4int* pVoxel;
  G4double m_BoxDimX;
  G4double m_BoxDimY;
  G4double m_BoxDimZ;
  G4int NumVoxelX;
  G4int NumVoxelZ;
  G4double dimVoxel;
  G4String m_SDName;
  G4int detectorChoice;
  G4double ExpHall_x ;
  G4double ExpHall_y ;
  G4double ExpHall_z ;

  BrachyPhantomSD* pPhantomSD;
  BrachyPhantomROGeometry*   pPhantomROGeometry;

  BrachyFactory* factory;

  BrachyMaterial* pMat;

private:

  BrachyDetectorMessenger* detectorMessenger; 

  G4Box*             ExpHall;        //pointer to the solid World 
  G4LogicalVolume*   ExpHallLog;     //pointer to the logical World
  G4VPhysicalVolume* ExpHallPhys;    //pointer to the physical World

  G4Box*              Phantom;
  G4LogicalVolume*    PhantomLog;
  G4VPhysicalVolume*  PhantomPhys;

  G4Material*         AbsorberMaterial;
};

#endif
