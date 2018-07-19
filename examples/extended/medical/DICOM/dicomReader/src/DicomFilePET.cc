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
#include "DicomFilePET.hh"
#include "DicomFileStructure.hh"
#include "DicomROI.hh"

#include "G4GeometryTolerance.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcpixel.h"
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmdata/dcpixseq.h"
#include "dcmtk/dcmrt/drtimage.h"

#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomFilePET::DicomFilePET()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomFilePET::DicomFilePET(DcmDataset* dset) : DicomVFileImage(dset)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePET::BuildActivities()
{
  G4int fCompress = theFileMgr->GetCompression();
  if( fNoVoxelX%fCompress != 0 || fNoVoxelY%fCompress != 0 ) {
    G4Exception("DicompFileMgr.:BuildMaterials",
                "DFC004",
                FatalException,
                ("Compression factor = " + std::to_string(fCompress) 
                 + " has to be a divisor of Number of voxels X = " + std::to_string(fNoVoxelX) 
                 + " and Y " + std::to_string(fNoVoxelY)).c_str());
  }

  //  if( DicomVerb(debugVerb) ) G4cout << " BuildMaterials " << fFileName << G4endl;
  double meanHV = 0.;
  for( int ir = 0; ir < fNoVoxelY; ir += fCompress ) {
    for( int ic = 0; ic < fNoVoxelX; ic += fCompress ) {
      meanHV = 0.;
      int isumrMax = std::min(ir+fCompress,fNoVoxelY);
      int isumcMax = std::min(ic+fCompress,fNoVoxelX);
      for( int isumr = ir; isumr < isumrMax; isumr ++ ) {
        for( int isumc = ic; isumc < isumcMax; isumc ++ ) {
          meanHV += fHounsfieldV[isumc+isumr*fNoVoxelX];
        }
      }
      meanHV /= (isumrMax-ir)*(isumcMax-ic);
      fActivities.push_back(meanHV);
    }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePET::DumpActivitiesToTextFile(std::ofstream& fout)
{
  G4int fCompress = theFileMgr->GetCompression();
  if( DicomFileMgr::verbose >= warningVerb ) G4cout << fLocation << " DumpDensitiesToTextFile " 
          << fFileName << " " << fActivities.size() << G4endl;
  
  G4int copyNo = 0;
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      fout << fActivities[ic+ir*fNoVoxelX/fCompress];
      if( ic != fNoVoxelX/fCompress-1) fout << " ";
      if( copyNo%8 == 7 ) fout << G4endl;
      copyNo++;
    }
    if( copyNo%8 != 0 ) fout << G4endl;
  }

}

