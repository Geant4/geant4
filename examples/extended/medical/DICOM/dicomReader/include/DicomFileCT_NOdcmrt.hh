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
#ifndef DicomFileCT_NOdcmrt__HH
#define DicomFileCT_NOdcmrt__HH

#include "DicomVFile.hh"
#include "DicomFileMgr.hh"
#include "G4ThreeVector.hh"

class DicomFileCT_NOdcmrt : public DicomVFile
{ 
public:
  DicomFileCT_NOdcmrt();
  DicomFileCT_NOdcmrt(DcmDataset* dset);
  ~DicomFileCT_NOdcmrt(){};

public:
  virtual void ReadData();

  void operator+=( const DicomFileCT_NOdcmrt& rhs );
  DicomFileCT_NOdcmrt operator+( const DicomFileCT_NOdcmrt& rhs );
  // add two slices that have the same dimensions, merging them in Z

  void BuildMaterials();
  void DumpHeaderToTextFile(std::ofstream& fout);
  void DumpMateIDsToTextFile(std::ofstream& fout);
  void DumpDensitiesToTextFile(std::ofstream& fout);
  void BuildStructureIDs();
  void DumpStructureIDsToTextFile(std::ofstream& fout);

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
  
  void SetNoVoxelX(const G4int& val) { fNoVoxelX = val; }
  void SetNoVoxelY(const G4int& val) { fNoVoxelY = val; }
  void SetNoVoxelZ(const G4int& val) { fNoVoxelZ = val; }
  
  void SetMinX(const G4double& val) { fMinX = val; };
  void SetMaxX(const G4double& val) { fMaxX = val; };
  void SetMinY(const G4double& val) { fMinY = val; };
  void SetMaxY(const G4double& val) { fMaxY = val; };
  void SetMinZ(const G4double& val) { fMinZ = val; };
  void SetMaxZ(const G4double& val) { fMaxZ = val; };
    
  const G4double& GetLocation() const { return fLocation; }
  
  void SetLocation(const G4double& val) { fLocation = val; }

  G4ThreeVector GetOrientationRows() const { return fOrientationRows; }
  G4ThreeVector GetOrientationColumns() const { return fOrientationColumns; }
  
  void DumpToTextFile();
  void DumpToBinaryFile();

  void ReadDataFromFile();

private:
  template <typename T> inline bool CheckConsistency(const T&, const T&, G4String);

  void ReadPixelData();
  void Print( std::ostream& out );
  
private:
  G4double fLocation;
  G4double fBitAllocated;
  G4double fRescaleSlope;
  G4double fRescaleIntercept;

  G4int fNoVoxelX, fNoVoxelY, fNoVoxelZ;  // number of voxels in each dimensions
  G4double fMinX,fMinY,fMinZ; // minimum extension of voxels (position of wall)
  G4double fMaxX,fMaxY,fMaxZ; // maximum extension of voxels (position of wall)
  G4double fVoxelDimX,fVoxelDimY,fVoxelDimZ; // maximum extension of voxels (position of wall)

  G4ThreeVector fOrientationRows;
  G4ThreeVector fOrientationColumns;
  
  std::vector<int> fHounsfieldV;
  std::vector<size_t> fMateIDs;
  std::vector<G4double> fDensities;
  std::vector<G4int> fStructure;
  //  G4int* fStructure;
  
  DicomFileMgr* theFileMgr;
};

//============================================================================
template <typename T>
inline bool DicomFileCT_NOdcmrt::CheckConsistency(const T& val1, const T& val2, 
                                                       G4String category) {
  if(val1 != val2) {
    G4Exception("DicomFileCT_NOdcmrtr::CheckConsistency", 
                "Consistency Mismatch : Keeping previous value if nonzero",
                JustWarning, category.c_str());
        return false;
  }
  return true;
}

#endif
