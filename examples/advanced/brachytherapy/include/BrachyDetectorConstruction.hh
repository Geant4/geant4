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
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************


#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class BrachyWaterBoxSD;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
 public:
	BrachyDetectorConstruction(G4String &SDName,G4int NumVoxelX,G4int NumVoxelZ);
	~BrachyDetectorConstruction();

 public:

  G4double GetBoxDim_X() {return  m_BoxDimX;}; 
  G4double GetBoxDim_Z() {return  m_BoxDimZ;};
  
 public:
	G4VPhysicalVolume* Construct();

private:
 const  G4int          m_NumVoxelX;
 const  G4int          m_NumVoxelZ;
  G4double        m_BoxDimX;
  G4double        m_BoxDimY;
  G4double        m_BoxDimZ;
  G4String m_SDName;
 };
#endif

