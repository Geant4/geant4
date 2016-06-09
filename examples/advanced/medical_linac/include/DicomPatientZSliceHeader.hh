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
// Author: P. Arce
// History: 30.11.07  First version
//*******************************************************
//
// DicomPatientZSliceHeader.hh :
//	- Contains the meta data information corresponding to one or several Z slices (number of voxels, dimension)
//*******************************************************
//
// The code is part of the DICOM extended example and it was modified by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#ifndef DicomPatientZSliceHeader_h
#define DicomPatientZSliceHeader_h 1

#include "globals.hh"
class G4material;
#include <fstream>
#include <vector>


class DicomPatientZSliceHeader 
{
public:

  DicomPatientZSliceHeader( const DicomPatientZSliceHeader& rhs );
  // build object copying an existing one (except Z dimensions)

  DicomPatientZSliceHeader( std::ifstream& fin );
  // build object reading data from a file

  ~DicomPatientZSliceHeader(){};

  // Get and set methods
  G4int GetNoVoxelX() const { return fNoVoxelX; };
  G4int GetNoVoxelY() const { return fNoVoxelY; };
  G4int GetNoVoxelZ() const { return fNoVoxelZ; };
  G4int GetNoVoxels() const { return fNoVoxelX*fNoVoxelY*fNoVoxelZ; };

  G4double GetMinX() const { return fMinX; };
  G4double GetMinY() const { return fMinY; };
  G4double GetMinZ() const { return fMinZ; };
  G4double GetMaxX() const { return fMaxX; };
  G4double GetMaxY() const { return fMaxY; };
  G4double GetMaxZ() const { return fMaxZ; };

  G4double GetVoxelHalfX() const { return (fMaxX-fMinX)/fNoVoxelX/2.; };
  G4double GetVoxelHalfY() const { return (fMaxY-fMinY)/fNoVoxelY/2.; };
  G4double GetVoxelHalfZ() const { return (fMaxZ-fMinZ)/fNoVoxelZ/2.; };

  std::vector<G4String> GetMaterialNames() const { return fMaterialNames; };
 

  void SetNoVoxelX(const G4int val) { fNoVoxelX = val; }
  void SetNoVoxelY(const G4int val) { fNoVoxelY = val; }
  void SetNoVoxelZ(const G4int val) { fNoVoxelZ = val; }

  void SetMinX(const G4double val) { fMinX = val; };
  void SetMaxX(const G4double val) { fMaxX = val; };
  void SetMinY(const G4double val) { fMinY = val; };
  void SetMaxY(const G4double val) { fMaxY = val; };
  void SetMinZ(const G4double val) { fMinZ = val; };
  void SetMaxZ(const G4double val) { fMaxZ = val; };

  void SetMaterialNames(std::vector<G4String>& mn ){ fMaterialNames = mn; }


  void operator+=( const DicomPatientZSliceHeader& rhs );  
  DicomPatientZSliceHeader operator+( const DicomPatientZSliceHeader& rhs );
  // add two slices that have the same dimensions, merging them in Z 

private:
  G4bool CheckMaterialExists( const G4String& mateName );
  // check that material read exists as a G4Material

private:
  G4int fNoVoxelX, fNoVoxelY, fNoVoxelZ;  // number of voxels in each dimensions
  G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position of wall)
  G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position of wall)

  std::vector<G4String> fMaterialNames; // list of material names
 
};

#endif
