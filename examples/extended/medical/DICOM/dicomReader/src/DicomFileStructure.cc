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
#include "DicomFileStructure.hh"

#include "G4ThreeVector.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmdata/dcpixel.h"
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmdata/dcpixseq.h"
#include "dcmtk/dcmrt/drtstrct.h"
#include "dcmtk/dcmrt/seq/drtrfors.h"  // for ReferencedFrameOfReferenceSequence
#include "dcmtk/dcmrt/seq/drtssrs.h"   // for StructureSetROISequence
#include "dcmtk/dcmrt/seq/drtrcs.h"      // for ROIContourSequence
#include "dcmtk/dcmrt/seq/drtcs.h"      // for ContourSequence
#include "dcmtk/dcmrt/seq/drtcis.h"      // for ContourImageSequence
#include "dcmtk/config/osconfig.h"   // make sure OS specific configuration is included 
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomFileStructure::DicomFileStructure(DcmDataset* dset) : DicomVFile(dset)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFileStructure::ReadData()
{
  DRTStructureSetIOD rtstruct;
  OFCondition result = rtstruct.read(*theDataset);
  //  DCMRT_INFO("Read RT Structure Set: " << status.text());
  if (!result.good()) {
    G4Exception("DicomFileStructure::ReadData",
                "DFS001",
                FatalException,
                result.text());
  }


  //@@@@@@@@@@@@ DRTReferencedFrameOfReferenceSequence
  DRTReferencedFrameOfReferenceSequence refSeq = rtstruct.getReferencedFrameOfReferenceSequence();
  if( refSeq.isEmpty() ) {
    G4Exception("DicomFileStructure::ReadData",
                "DFS002",
                JustWarning,
                "DRTReferencedFrameOfReferenceSequence is empty");
  }

  G4cout << "@@@@@ NUMBER OF ReferenceSequences " << refSeq.getNumberOfItems() << G4endl;
  refSeq.gotoFirstItem();
  for( size_t i1 = 0; i1 < refSeq.getNumberOfItems(); i1++ ) {
    DRTReferencedFrameOfReferenceSequence::Item &item = refSeq.getCurrentItem();
    OFString uid;
    item.getFrameOfReferenceUID(uid);
    G4cout << " FrameOfReferenceUID " << uid << G4endl;
    DRTRTReferencedStudySequence &reference_study_sequence_ref = 
     item.getRTReferencedStudySequence();   
    G4cout << "@@@@ NUMBER OF ReferenceStudySequences " 
           << reference_study_sequence_ref.getNumberOfItems() << G4endl;
    reference_study_sequence_ref.gotoFirstItem();
    for( size_t i2 = 0; i2 < reference_study_sequence_ref.getNumberOfItems(); i2++ ) {
      DRTRTReferencedStudySequence::Item &rss_item = reference_study_sequence_ref.getCurrentItem();
      DRTRTReferencedSeriesSequence &series_seq_ref = rss_item.getRTReferencedSeriesSequence();
      G4cout << "@@@ NUMBER OF SeriesSequences " << series_seq_ref.getNumberOfItems() << G4endl;
      series_seq_ref.gotoFirstItem();
      for( size_t i3 = 0; i3 < series_seq_ref.getNumberOfItems(); i3++ ) {
        DRTRTReferencedSeriesSequence::Item &ref_series_seq_item = series_seq_ref.getCurrentItem();
        DRTContourImageSequence &image_sequence_seq_ref = 
         ref_series_seq_item.getContourImageSequence();
        G4cout << "@@ NUMBER OF ContourImageSequences "
               << image_sequence_seq_ref.getNumberOfItems() << G4endl;
        image_sequence_seq_ref.gotoFirstItem();
        for( size_t i4 = 0; i4 < image_sequence_seq_ref.getNumberOfItems(); i4++ ) {
          DRTContourImageSequence::Item &image_contour_item = 
           image_sequence_seq_ref.getCurrentItem();
          OFString refSOPInstUID;
          image_contour_item.getReferencedSOPInstanceUID(refSOPInstUID);
          std::cout <<"ReferencedSOPInstanceUID= " << refSOPInstUID << std::endl;
          image_sequence_seq_ref.gotoNextItem().good();
        } // end if image_sequence_seq_ref
        series_seq_ref.gotoNextItem();
      } // end if series_seq_ref good
      reference_study_sequence_ref.gotoNextItem();
    } // end if reference_study_sequence_ref good
    refSeq.gotoNextItem();
  } // end if refSeq.first item
  
  //@@@@@@@@@@@@   DRTROISequence 
  DRTStructureSetROISequence ROISeq = rtstruct.getStructureSetROISequence();
  G4cout << "@@@@@ NUMBER OF ROISequences " << ROISeq.getNumberOfItems() << G4endl;
  for( size_t i1 = 0; i1 < ROISeq.getNumberOfItems(); i1++ ) {
    DRTStructureSetROISequence::Item &item = ROISeq.getCurrentItem();
    OFString ROIName, ROINumber, ROIGenerationAlgorithm;
    item.getROINumber(ROINumber);
    item.getROIName(ROIName);
    item.getROIGenerationAlgorithm(ROIGenerationAlgorithm);
    if( ROINumber != "" ) {
      DicomROI* roi = new DicomROI(G4UIcommand::ConvertToInt(ROINumber.c_str()), ROIName.c_str());
      theROIs.push_back( roi );
      G4cout << " ROI: " << roi->GetNumber() << " " << roi->GetName() << " "
             << ROIGenerationAlgorithm << G4endl;
    }

    ROISeq.gotoNextItem().good();
  } // end if ROISeq.first item


  //@@@@@@@@@@@@   DRTROIContourSequence
  DRTROIContourSequence ROIContourSeq = rtstruct.getROIContourSequence();
  if( ROISeq.getNumberOfItems() != ROIContourSeq.getNumberOfItems() ) {
    G4Exception("DicomFileStructure",
                "DCS0001",
                FatalException,
                ("Different number of ROIs and ROI Contours " 
                 + std::to_string(ROISeq.getNumberOfItems()) + " <> " 
                 + std::to_string(ROIContourSeq.getNumberOfItems())).c_str());
  } 

  ROIContourSeq.gotoFirstItem();
  for( size_t i1 = 0; i1 < ROIContourSeq.getNumberOfItems(); i1++ ) {
    DRTROIContourSequence::Item &item = ROIContourSeq.getCurrentItem();
    OFString displayColor;
    item.getROIDisplayColor(displayColor);
    //    G4cout << " ROIDisplayColor " << displayColor << G4endl;
    
    DRTContourSequence contour_seq = item.getContourSequence();
    //    G4cout << "@@@@ NUMBER OF ContourSequences " << contour_seq.getNumberOfItems() << G4endl;
    contour_seq.gotoFirstItem();
    for( size_t i2 = 0; i2 < contour_seq.getNumberOfItems(); i2++ ) {
      //      if (contour_seq.gotoFirstItem().good()) {
      //        do {
      DRTContourSequence::Item &cs_item = contour_seq.getCurrentItem();
      
      DicomROIContour* roiC = new DicomROIContour();
      
      DRTContourImageSequence &contour_image_seq = cs_item.getContourImageSequence();
      
      contour_image_seq.gotoFirstItem();
      for( size_t i3 = 0; i3 < contour_image_seq.getNumberOfItems(); i3++ ) {
        DRTContourImageSequence::Item cis_item = contour_image_seq.getCurrentItem();
        OFString refSOPCUID;
        cis_item.getReferencedSOPClassUID(refSOPCUID);
        OFString refSOPIUID;
        cis_item.getReferencedSOPInstanceUID(refSOPIUID);
        if( refSOPIUID != "") roiC->AddImageIUID(refSOPIUID.c_str());
        contour_image_seq.gotoNextItem();
      } // end if contour_image_seq
      
        //@@@
      OFString geomType;
      cs_item.getContourGeometricType(geomType);
      Sint32 nPoints;
      cs_item.getNumberOfContourPoints(nPoints);
      roiC->SetGeomType(geomType);
      OFVector<Float64> data;
      cs_item.getContourData(data);
      std::vector<G4ThreeVector> dataV;
      for( Sint32 ii = 0; ii < nPoints*3; ii++ ) {
        if( ii%3 == 2 ) dataV.push_back( G4ThreeVector( data[ii-2], data[ii-1], data[ii] ) );
      }
      roiC->SetData(dataV);
      theROIs[i1]->AddContour(roiC);

      contour_seq.gotoNextItem();
    }
    ROIContourSeq.gotoNextItem();
  } // end if ROIContourSeq.first item

  //@@@@ Print ROIs
  G4cout << " @@@@@@@@@@@ ROIs " << G4endl;
  for( size_t ii = 0; ii < theROIs.size(); ii++ ) {
    theROIs[ii]->Print(G4cout);
  }
}

