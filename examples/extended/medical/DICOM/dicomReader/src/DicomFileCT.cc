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
#include "DicomFileCT.hh"
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
DicomFileCT::DicomFileCT()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomFileCT::DicomFileCT(DcmDataset* dset) : DicomVFileImage(dset)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileCT::BuildMaterials()
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
 // G4cout << isumr << " " << isumc << " GET mean " << meanHV << G4endl;
        }
      }
      meanHV /= (isumrMax-ir)*(isumcMax-ic);
      G4double meanDens = theFileMgr->Hounsfield2density(meanHV);
      //      G4cout << ir << " " << ic << " FINAL mean " << meanDens << G4endl;

      fDensities.push_back(meanDens);
      size_t mateID;
      if( theFileMgr->IsMaterialsDensity() ) {
        mateID = theFileMgr->GetMaterialIndexByDensity(meanDens);
      } else {
        mateID = theFileMgr->GetMaterialIndex(meanHV);
      }
      fMateIDs.push_back(mateID);
    }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileCT::DumpMateIDsToTextFile(std::ofstream& fout)
{
  G4int fCompress = theFileMgr->GetCompression();
  if( DicomFileMgr::verbose >= warningVerb ) G4cout << fLocation << " DumpMateIDsToTextFile " 
            << fFileName << " " << fMateIDs.size() << G4endl;
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      fout << fMateIDs[ic+ir*fNoVoxelX/fCompress];
      if( ic != fNoVoxelX/fCompress-1) fout << " ";
    }
    fout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileCT::DumpDensitiesToTextFile(std::ofstream& fout)
{
  G4int fCompress = theFileMgr->GetCompression();
   if( DicomFileMgr::verbose >= warningVerb ) G4cout << fLocation << " DumpDensitiesToTextFile " 
          << fFileName << " " << fDensities.size() << G4endl;
  
   G4int copyNo = 0;
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      fout << fDensities[ic+ir*fNoVoxelX/fCompress];
      if( ic != fNoVoxelX/fCompress-1) fout << " ";
      if( copyNo%8 == 7 ) fout << G4endl;
      copyNo++;
    }
    if( copyNo%8 != 0 ) fout << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileCT::BuildStructureIDs()
{
  G4int fCompress = theFileMgr->GetCompression();
  std::vector<DicomFileStructure*> dfs = theFileMgr->GetStructFiles();
  if( dfs.size() == 0 ) return;

  G4int NMAXROI = DicomFileMgr::GetInstance()->GetStructureNMaxROI();
  G4int NMAXROI_real = std::log(INT_MAX)/std::log(NMAXROI);
    
  //  fStructure = new G4int(fNoVoxelX*fNoVoxelY);
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      //      fStructure[ic+ir*fNoVoxelX] = 0;
      fStructure.push_back(0);
    }
  }

  std::set<double> distInters;

  //  std::fill_n(fStructure,fNoVoxelX*fNoVoxelY,0);
  //
  for( size_t ii = 0; ii < dfs.size(); ii++ ){
    std::vector<DicomROI*> rois = dfs[ii]->GetROIs();
    for( size_t jj = 0; jj < rois.size(); jj++ ){
      if( DicomFileMgr::verbose >= debugVerb ) std::cout << " BuildStructureIDs checking ROI "
            << rois[jj]->GetNumber() << " " << rois[jj]->GetName() << G4endl;
      std::vector<DicomROIContour*> contours = rois[jj]->GetContours();
      //      G4int roi = jj+1;
      G4int roiID = rois[jj]->GetNumber();
      for( size_t kk = 0; kk < contours.size(); kk++ ){
        // check contour corresponds to this CT slice
        DicomROIContour* roic = contours[kk];
      // if( DicomVerb(-debugVerb) ) G4cout << jj << " " << kk << " INTERS CONTOUR " << roic 
        //                << " " << fLocation << G4endl;
        if( DicomFileMgr::verbose >= debugVerb ) G4cout << " " << roiID << " " << kk 
            << " INTERS CONTOUR " << roic->GetZ() << " SLICE Z " << fMinZ << " " << fMaxZ << G4endl;
        // Check Contour correspond to slice

        if( roic->GetZ() > fMaxZ || roic->GetZ() < fMinZ ) continue;
        if( DicomFileMgr::verbose >= debugVerb ) G4cout << " INTERS CONTOUR WITH Z SLIZE " 
            << fMinZ << " < " << roic->GetZ() << " < " << fMaxZ << G4endl;
        if( roic->GetGeomType() == "CLOSED_PLANAR" ){
          // Get min and max X and Y, and corresponding slice indexes
          std::vector<G4ThreeVector> points = roic->GetPoints();
           if( DicomFileMgr::verbose >= debugVerb ) G4cout << jj << " " << kk << " NPOINTS " 
                 << points.size() << G4endl;
          std::vector<G4ThreeVector> dirs = roic->GetDirections();
          double minXc = DBL_MAX;
          double maxXc = -DBL_MAX;
          double minYc = DBL_MAX;
          double maxYc = -DBL_MAX;
          for( size_t ll = 0; ll < points.size(); ll++ ){
            minXc = std::min(minXc,points[ll].x());
            maxXc = std::max(maxXc,points[ll].x());
            minYc = std::min(minYc,points[ll].y());
            maxYc = std::max(maxYc,points[ll].y());
          }
          if( minXc < fMinX || maxXc > fMaxX || minYc < fMinY || maxYc > fMaxY ){
            G4cerr << " minXc " << minXc << " < " << fMinX
                   << " maxXc " << maxXc << " > " << fMaxX
                   << " minYc " << minYc << " < " << fMinY
                   << " maxYc " << maxYc << " > " << fMaxY << G4endl;
              G4Exception("DicomFileCT::BuildStructureIDs",
                          "DFCT001",
                          JustWarning,
                          "Contour limits exceed Z slice extent");
          }
          int idMinX = std::max(0,int((minXc-fMinX)/fVoxelDimX/fCompress));
          int idMaxX = std::min(fNoVoxelX/fCompress-1,int((maxXc-fMinX)/fVoxelDimX/fCompress+1));
          int idMinY = std::max(0,int((minYc-fMinY)/fVoxelDimY/fCompress));
          int idMaxY = std::min(fNoVoxelY/fCompress-1,int((maxYc-fMinY)/fVoxelDimY/fCompress+1));
          if( DicomFileMgr::verbose >= debugVerb )
            G4cout << " minXc " << minXc << " < " << fMinX
                 << " maxXc " << maxXc << " > " << fMaxX
                 << " minYc " << minYc << " < " << fMinY
                 << " maxYc " << maxYc << " > " << fMaxY << G4endl;
          if( DicomFileMgr::verbose >= debugVerb )
            G4cout << " idMinX " << idMinX 
                 << " idMaxX " << idMaxX 
                 << " idMinY " << idMinY 
                 << " idMaxY " << idMaxY << G4endl;
          //for each voxel: build 4 lines from the corner towards the center 
          // and check how many contour segments it crosses, and the minimum distance to a segment
          for( int ix = idMinX; ix <= idMaxX; ix++ ) {
            for( int iy = idMinY; iy <= idMaxY; iy++ ) {
              G4bool bOK = false;
              G4int bOKs;
               if( DicomFileMgr::verbose >= debugVerb ) G4cout << "@@@@@ TRYING POINT (" 
                    <<  fMinX + fVoxelDimX*fCompress*(ix+0.5) << "," 
                    <<  fMinY + fVoxelDimY*fCompress*(iy+0.5) << ")" << G4endl;
              // four corners
              for( int icx = 0; icx <= 1; icx++ ){
                for( int icy = 0; icy <= 1; icy++ ){
                  bOKs = 0;
                  if( bOK ) continue;
                  double p0x = fMinX + fVoxelDimX*fCompress * (ix+icx);
                  double p0y = fMinY + fVoxelDimY*fCompress * (iy+icy);
                  double v0x = 1.;
                  if( icx == 1 ) v0x = -1.;
                  double v0y = 0.99*fVoxelDimY/fVoxelDimX*std::pow(-1.,icy);
                  if( DicomFileMgr::verbose >= testVerb ) G4cout << ix << " + " << icx << " " 
                      << iy << " + " << icy << " CORNER (" << p0x << "," << p0y << ") " 
                      << " DIR= (" << v0x << "," << v0y << ") " << G4endl;
                  int NTRIES = theFileMgr->GetStructureNCheck();
                 for( int ip = 0; ip < NTRIES; ip++) {
                   distInters.clear();
                   v0y -= ip*0.1*v0y;
                   G4double halfDiagonal = 0.5*fVoxelDimX*fCompress;
                   if( DicomFileMgr::verbose >= testVerb ) G4cout << ip 
                            << " TRYING WITH DIRECTION (" << " DIR= (" << v0x << "," 
                            << v0y << ") " << bOKs << G4endl;
                  for( size_t ll = 0; ll < points.size(); ll++ ){
                    double d0x = points[ll].x() - p0x;
                    double d0y = points[ll].y() - p0y;
                    double w0x = dirs[ll].x();
                    double w0y = dirs[ll].y();
                    double fac1 = w0x*v0y - w0y*v0x;
                    if( fac1 == 0 ) { // parallel lines
                      continue;
                    }
                    double fac2 = d0x*v0y - d0y*v0x;
                    double fac3 = d0y*w0x - d0x*w0y;
                    double lambdaq = -fac2/fac1;
                    if( lambdaq < 0. || lambdaq >= 1. ) continue; 
                              // intersection further than segment length
                    double lambdap = fac3/fac1;
                    if( lambdap > 0. ) {
                      distInters.insert(lambdap);
                      if( DicomFileMgr::verbose >= testVerb ) G4cout << " !! GOOD INTERS " 
                          <<lambdaq << "  (" << d0x+p0x+lambdaq*w0x << "," << d0y+p0y+lambdaq*w0y 
                          << ")  =  (" << p0x+lambdap*v0x << "," << p0y+lambdap*v0y << ") " 
                          << " N " << distInters.size() << G4endl;                      
                    }
                    if( DicomFileMgr::verbose >= testVerb ) G4cout << " INTERS LAMBDAQ " 
                          << lambdaq << " P " << lambdap << G4endl;

                    if( DicomFileMgr::verbose >= debugVerb ) G4cout << " INTERS POINT (" 
                          << d0x+p0x+lambdaq*w0x << "," << d0y+p0y+lambdaq*w0y << ")  =  (" 
                          << p0x+lambdap*v0x << "," << p0y+lambdap*v0y << ") " << G4endl;
                  }
                  if( distInters.size() % 2 == 1 ) {
                    if( *(distInters.begin() ) > halfDiagonal ) {
                      //                      bOK = true;
                      bOKs += 1;
                      if( DicomFileMgr::verbose >= debugVerb ) G4cout << " OK= " << bOKs 
                            << " N INTERS " << distInters.size() << " " << *(distInters.begin()) 
                            << " > " << halfDiagonal <<  G4endl;
                    }
                  }
                 }
                 if(bOKs == NTRIES ) {
                    bOK = true;
                    if( DicomFileMgr::verbose >= debugVerb ) G4cout << "@@@@@ POINT OK (" 
                           <<  fMinX + fVoxelDimX*fCompress*(ix+0.5) << "," 
                           <<  fMinY + fVoxelDimY*fCompress*(iy+0.5) << ")" << G4endl;
                 }

                }
              } // loop to four corners 
              if( bOK ) {
                // extract previous ROI value
                int roival = fStructure[ix+iy*fNoVoxelX/fCompress];
                //                roival = 2 + NMAXROI*3 + NMAXROI*NMAXROI*15;
                if(roival != 0 && roival != int(roiID) ) {
                  std::set<G4int> roisVoxel;
                  int nRois = std::log10(roival)/std::log10(NMAXROI)+1;
                  if( nRois > NMAXROI_real ) { // another one cannot be added
                    G4Exception("DicomFileCT:BuildStructureIDs",
                                "DFC0004",
                                FatalException,
                                G4String("Too many ROIs associated to a voxel: \
" + std::to_string(nRois) + " > " + std::to_string(NMAXROI_real) + " TRY CHAN\
GING -NStructureNMaxROI argument to a lower value").c_str());
                  }
                  for( int inr = 0; inr < nRois; inr++ ) {
                    roisVoxel.insert( int(roival/std::pow(NMAXROI,inr))%NMAXROI );
                  }
                  roisVoxel.insert(roiID);
                  roival = 0;
                  size_t inr = 0;
                  for( std::set<G4int>::const_iterator ite = roisVoxel.begin(); 
                   ite != roisVoxel.end(); ite++, inr++ ) {
                    roival += (*ite)*std::pow(NMAXROI,inr);
                  }
                  fStructure[ix+iy*fNoVoxelX/fCompress] = roival;
                  if( DicomFileMgr::verbose >= testVerb ){
                    G4cout << " WITH PREVIOUS ROI IN VOXEL " << roival << G4endl;
                  }
                } else {
                  fStructure[ix+iy*fNoVoxelX/fCompress] = roiID;
                }
              } 
              
            }
          }
        }
      }
    }
  }

  if( DicomFileMgr::verbose >= 1 ) {
  //@@@@ PRINT structures
  //@@@ PRINT points of structures in this Z slice
   if( DicomFileMgr::verbose >= 0 ) G4cout << " STRUCTURES FOR SLICE " << fLocation << G4endl;
   for( size_t ii = 0; ii < dfs.size(); ii++ ){
    std::vector<DicomROI*> rois = dfs[ii]->GetROIs();
    for( size_t jj = 0; jj < rois.size(); jj++ ){
      int roi = rois[jj]->GetNumber(); 
      std::vector<DicomROIContour*> contours = rois[jj]->GetContours();
      for( size_t kk = 0; kk < contours.size(); kk++ ){
        DicomROIContour* roic = contours[kk];
        // check contour corresponds to this CT slice
         if( roic->GetZ() > fMaxZ || roic->GetZ() < fMinZ ) continue;
        if( roic->GetGeomType() == "CLOSED_PLANAR" ){
          if( DicomFileMgr::verbose >= testVerb ) G4cout << " PRINTING CONTOUR? " << roi << " " 
                << kk << " INTERS CONTOUR " << roic->GetZ() << " SLICE Z " 
                << fMinZ << " " << fMaxZ << G4endl;
          if( roic->GetZ() > fMaxZ || roic->GetZ() < fMinZ ) continue;
          std::vector<G4ThreeVector> points = roic->GetPoints();
          std::vector<G4ThreeVector> dirs = roic->GetDirections();
           if( DicomFileMgr::verbose >= 0 ) std::cout << " STRUCTURE Z " << roic->GetZ() 
                << std::endl;
          for( size_t ll = 0; ll < points.size(); ll++ ){
            if( DicomFileMgr::verbose >= 0 ) G4cout << roi << " : " << ll 
                << " STRUCTURE POINT (" << points[ll].x() << "," << points[ll].y() << ") " 
                << " (" << dirs[ll].x() << "," << dirs[ll].y() << ") = " << roi << G4endl;
          }
        }
      }
    }
  }
  //@@@ PRINT points in slice inside structure
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      if( fStructure[ic+ir*fNoVoxelX/fCompress] != 0 ) {
         if( DicomFileMgr::verbose >= 0 ) G4cout << ic+ir*fNoVoxelX/fCompress << " = " << ic 
                << " " << ir << " STRUCTURE VOXEL (" << fMinX + fVoxelDimX*fCompress * (ic+0.5) 
                << "," << fMinY + fVoxelDimY*fCompress * (ir+0.5) << ") = " 
                << fStructure[ic+ir*fNoVoxelX/fCompress] << G4endl;
      }
    }
  }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileCT::DumpStructureIDsToTextFile(std::ofstream& fout)
{
  G4int fCompress = theFileMgr->GetCompression();
   if( DicomFileMgr::verbose >= 0 ) G4cout << fLocation << " DumpStructureIDsToTextFile " 
          << fFileName << " " << fStructure.size() << G4endl;
  std::vector<DicomFileStructure*> dfs = theFileMgr->GetStructFiles();
  if( dfs.size() == 0 ) return;
  
  for( int ir = 0; ir < fNoVoxelY/fCompress; ir++ ) {
    for( int ic = 0; ic < fNoVoxelX/fCompress; ic++ ) {
      fout << fStructure[ic+ir*fNoVoxelX/fCompress];
      if( ic != fNoVoxelX/fCompress-1) fout << " ";
    }
    fout << G4endl;
  }
}

