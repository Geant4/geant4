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
//    ************************************
//    *                                  *
//    *    BrachyWaterBoxROGeometry.hh   *
//    *                                  *
//    ************************************


#ifndef BrachyWaterBoxROGeometry_h
#define BrachyWaterBoxROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class BrachyWaterBoxROGeometry : public G4VReadOutGeometry
{
  public:
	BrachyWaterBoxROGeometry(G4String aString,G4double DetDimX,G4double DetDimZ,G4int NumVoxelX,G4int NumVoxelZ);
  	~BrachyWaterBoxROGeometry();

  public:
	const G4double m_DetDimX;
	const G4double m_DetDimZ;
	const G4int m_NumVoxelX;
	const G4int m_NumVoxelZ;

  private:
  	G4VPhysicalVolume* Build();
};

#endif
