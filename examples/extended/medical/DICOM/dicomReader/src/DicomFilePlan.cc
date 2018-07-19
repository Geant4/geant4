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
#include "DicomFilePlan.hh"
#include "DicomBeam.hh"
#include "DicomBeamDeviceRef.hh"
#include "DicomBeamDevicePos.hh"
#include "DicomBeamControlPoint.hh"
#include "DicomBeamCompensator.hh"
#include "DicomBeamBlock.hh"
#include "DicomBeamWedge.hh"

#include "G4ThreeVector.hh"

#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcdeftag.h"
#include "dcmtk/dcmrt/drtplan.h"
#include "dcmtk/dcmrt/seq/drtfgs.h" // DRTFractionGroupSequence
#include "dcmtk/dcmrt/seq/drtrbs8.h" // DRTReferencedBeamSequenceInRTFractionSchemeModule
#include "dcmtk/dcmrt/seq/drtbs.h"     // for BeamSequence
#include "dcmtk/dcmrt/seq/drtblds1.h"  // for BeamLimitingDeviceSequence
#include "dcmtk/dcmrt/seq/drtcps.h"    // for ControlPointSequence
#include "dcmtk/dcmrt/seq/drtbldps.h"  // for BeamLimitingDevicePositionSequence
#include "dcmtk/dcmrt/seq/drtcos.h" // for CompensatorSequence
#include "dcmtk/dcmrt/seq/drtbl2.h" // for BlockSequence
#include "dcmtk/dcmrt/seq/drtws.h"     // for WedgeSequence

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomFilePlan::DicomFilePlan(DcmDataset* dset) : DicomVFile(dset)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePlan::ReadData()
{
  DRTPlanIOD rtplan;
  OFCondition result = rtplan.read(*theDataset);
  if (!result.good()) {
    G4Exception("DicomFilePlan::ReadData",
                "DFS001",
                FatalException,
                result.text());
  }
  OFString fstr;
  Sint32 fint;
  Float64 ffloat;
  OFVector<Float64> fvfloat; 
  
  DRTFractionGroupSequence frgSeq = rtplan.getFractionGroupSequence();
  if( frgSeq.isEmpty() ) {
    G4Exception("DicomFilePlan::ReadData",
                "DFS002",
                JustWarning,
                "DRTFractionGroupSequence is empty");
  }
  G4cout << "@@@@@ NUMBER OF FractionSequences " << frgSeq.getNumberOfItems() << G4endl;
  frgSeq.gotoFirstItem();
  for( size_t i1 = 0; i1 < frgSeq.getNumberOfItems(); i1++ ) {
    DRTFractionGroupSequence::Item &rfgItem = frgSeq.getCurrentItem();
    rfgItem.getBeamDoseMeaning(fstr);
    G4cout << " " << i1 << " BeamDoseMeaning " << fstr << G4endl;
    rfgItem.getFractionGroupDescription(fstr);
    G4cout << " " << i1 << " FractionGroupDescription " << fstr << G4endl;
    rfgItem.getFractionGroupNumber(fint);
    G4cout << " " << i1 << " FractionGroupNumber " << fint << G4endl;
    rfgItem.getFractionPattern(fstr);
    G4cout << " " << i1 << " FractionPattern " << fstr << G4endl;
    rfgItem.getNumberOfBeams(fint);
    G4cout << " " << i1 << " NumberOfBeams " << fint << G4endl;
    theNumberOfBeams = fint;
    rfgItem.getNumberOfBrachyApplicationSetups(fint);
    G4cout << " " << i1 << " NumberOfBrachyApplicationSetups " << fint << G4endl;
    CheckData0(" NumberOfBrachyApplicationSetups ", fint);
    rfgItem.getNumberOfFractionPatternDigitsPerDay(fint);
    G4cout << " " << i1 << " NumberOfFractionPatternDigitsPerDay " << fint << G4endl;
    rfgItem.getNumberOfFractionsPlanned(fint);
    G4cout << " " << i1 << " NumberOfFractionsPlanned " << fint << G4endl;
    rfgItem.getRepeatFractionCycleLength(fint);
    G4cout << " " << i1 << " RepeatFractionCycleLength " << fint << G4endl;
    DRTReferencedBeamSequenceInRTFractionSchemeModule refBeamSeq = 
     rfgItem.getReferencedBeamSequence();
    G4cout << " @@@ NUMBER OF ReferencedBeamSequences " << refBeamSeq.getNumberOfItems() << G4endl;

    refBeamSeq.gotoFirstItem();
    for( size_t i2 = 0; i2 < refBeamSeq.getNumberOfItems(); i2++ ) {
      DicomBeam* db = new DicomBeam();
      theBeams.push_back(db);
      DRTReferencedBeamSequenceInRTFractionSchemeModule::Item &rbsItem = 
       refBeamSeq.getCurrentItem();
      rbsItem.getBeamDeliveryDurationLimit(ffloat);
      G4cout << "  " << i2 << " BeamDeliveryDurationLimit " << ffloat << G4endl;
      rbsItem.getBeamDose(ffloat);
      G4cout << "  " << i2 << " BeamDose " << ffloat << G4endl; // dose at dose point
      rbsItem.getBeamDoseSpecificationPoint(fvfloat);
      G4cout << "  " << i2 << " BeamDoseSpecificationPoint (" << fvfloat[0] << "," << fvfloat[1] 
             << "," << fvfloat[2] << ")" << G4endl;
      db->SetDoseSpecificationPoint(G4ThreeVector(fvfloat[0],fvfloat[1],fvfloat[2]));
      rbsItem.getBeamMeterset(ffloat);
      G4cout << "  " << i2 << " BeamMeterset " << ffloat << G4endl;
      db->SetMeterset(ffloat);
      rbsItem.getReferencedBeamNumber(fint);
      G4cout << "  " << i2 << " ReferencedBeamNumber " << fint << G4endl;
    
      refBeamSeq.gotoNextItem();
    }

    frgSeq.gotoNextItem();
  }


  DRTBeamSequence beamSeq = rtplan.getBeamSequence();
  if( beamSeq.isEmpty() ) {
    G4Exception("DicomFilePlan::ReadData",
                "DFS002",
                JustWarning,
                "DRTBeamSequence is empty");
  }
  G4cout << "@@@@@ NUMBER OF BeamSequences " << beamSeq.getNumberOfItems() << G4endl;
  beamSeq.gotoFirstItem();
  for( size_t i1 = 0; i1 < beamSeq.getNumberOfItems(); i1++ ) {
    DicomBeam* db = theBeams[i1];
    DRTBeamSequence::Item &beamItem = beamSeq.getCurrentItem();

    beamItem.getManufacturer(fstr);
    G4cout << " " << i1 << " Manufacturer " << fstr << G4endl;
    beamItem.getManufacturerModelName(fstr);
    G4cout << " " << i1 << " ManufacturerModelName " << fstr << G4endl;
    beamItem.getTreatmentMachineName(fstr);
    G4cout << " " << i1 << " TreatmentMachineName " << fstr << G4endl;
    beamItem.getPrimaryDosimeterUnit(fstr);
    G4cout << " " << i1 << " PrimaryDosimeterUnit " << fstr << G4endl;
    beamItem.getSourceAxisDistance(ffloat);
    G4cout << " " << i1 << " SourceAxisDistance " << ffloat << G4endl;
    db->SetSourceAxisDistance(ffloat);
    
    DRTBeamLimitingDeviceSequenceInRTBeamsModule beamLDS = beamItem.getBeamLimitingDeviceSequence();
    G4cout << " @@@ NUMBER OF BeamLimitingDeviceSequence " << beamLDS.getNumberOfItems() << G4endl;
    beamLDS.gotoFirstItem();
    for( size_t i2 = 0; i2 < beamLDS.getNumberOfItems(); i2++ ) {
      DRTBeamLimitingDeviceSequenceInRTBeamsModule::Item bldsItem = beamLDS.getCurrentItem();
      DicomBeamDeviceRef* dbd = new DicomBeamDeviceRef(bldsItem);
      db->AddDevice(dbd);
      
      beamLDS.gotoNextItem();
    }
 
    beamItem.getBeamNumber(fint);
    G4cout << " " << i1 << " BeamNumber " << fint << G4endl;
    db->SetNumber(fint);
    beamItem.getBeamName(fstr);
    G4cout << " " << i1 << " BeamName " << fstr << G4endl;
    beamItem.getBeamDescription(fstr);
    G4cout << " " << i1 << " BeamDescription " << fstr << G4endl;
    beamItem.getBeamType(fstr);
    G4cout << " " << i1 << " BeamType " << fstr << G4endl;
    beamItem.getRadiationType(fstr);
    G4cout << " " << i1 << " RadiationType " << fstr << G4endl;
    db->SetRadiationType(fstr);
    beamItem.getTreatmentDeliveryType(fstr);
    G4cout << " " << i1 << " TreatmentDeliveryType " << fstr << G4endl;
    beamItem.getNumberOfWedges(fint);
    G4cout << " " << i1 << " NumberOfWedges " << fint << G4endl;
    DRTWedgeSequence beamWedge = beamItem.getWedgeSequence();
    beamWedge.gotoFirstItem();
    for( size_t i2 = 0; i2 < beamWedge.getNumberOfItems(); i2++ ) {
      DRTWedgeSequence::Item bwedItem = beamWedge.getCurrentItem();
      DicomBeamWedge* dbwed = new DicomBeamWedge( bwedItem );
      db->AddWedge( dbwed );
      beamWedge.gotoNextItem();
    }
    
    beamItem.getNumberOfCompensators(fint);
    G4cout << " " << i1 << " NumberOfCompensators " << fint << G4endl;
    DRTCompensatorSequence beamCompens = beamItem.getCompensatorSequence();
    beamCompens.gotoFirstItem();
    for( size_t i2 = 0; i2 < beamCompens.getNumberOfItems(); i2++ ) {
      DRTCompensatorSequence::Item bcompItem = beamCompens.getCurrentItem();
      DicomBeamCompensator* dbcomp = new DicomBeamCompensator( bcompItem );
      db->AddCompensator( dbcomp );
      beamCompens.gotoNextItem();
    }
  
    beamItem.getNumberOfBoli(fint);
    G4cout << " " << i1 << " NumberOfBoli " << fint << G4endl;
    //Bolus has no relevant info (see drtrbos1.h)

    beamItem.getNumberOfBlocks(fint);
    G4cout << " " << i1 << " NumberOfBlocks " << fint << G4endl;
    DRTBlockSequenceInRTBeamsModule beamBlock = beamItem.getBlockSequence();
    beamBlock.gotoFirstItem();
    for( size_t i2 = 0; i2 < beamBlock.getNumberOfItems(); i2++ ) {
      DRTBlockSequenceInRTBeamsModule::Item bblItem = beamBlock.getCurrentItem();
      DicomBeamBlock* dbbl = new DicomBeamBlock( bblItem );
      db->AddBlock( dbbl );
      beamBlock.gotoNextItem();
    }

    beamItem.getFinalCumulativeMetersetWeight(fstr);
    G4cout << " " << i1 << " FinalCumulativeMetersetWeight " << fstr << G4endl;
    beamItem.getDeviceSerialNumber(fstr);
    G4cout << " " << i1 << " DeviceSerialNumber " << fstr << G4endl;
    beamItem.getHighDoseTechniqueType(fstr);
    G4cout << " " << i1 << " HighDoseTechniqueType " << fstr << G4endl;
    beamItem.getInstitutionAddress(fstr);
    G4cout << " " << i1 << " InstitutionAddress " << fstr << G4endl;
    beamItem.getInstitutionName(fstr);
    G4cout << " " << i1 << " InstitutionName " << fstr << G4endl;
    beamItem.getInstitutionalDepartmentName(fstr);
    G4cout << " " << i1 << " InstitutionalDepartmentName " << fstr << G4endl;
    beamItem.getReferencedPatientSetupNumber(fint);
    G4cout << " " << i1 << " ReferencedPatientSetupNumber " << fint << G4endl;
    beamItem.getReferencedToleranceTableNumber(fint);
    G4cout << " " << i1 << " ReferencedToleranceTableNumber " << fint << G4endl;
    beamItem.getTotalBlockTrayFactor(ffloat);
    G4cout << " " << i1 << " TotalBlockTrayFactor " << ffloat << G4endl;
    beamItem.getTotalCompensatorTrayFactor(ffloat);
    G4cout << " " << i1 << " TotalCompensatorTrayFactor " << ffloat << G4endl;
    
    beamItem.getNumberOfControlPoints(fint);
    DRTControlPointSequence controlPSeq = beamItem.getControlPointSequence();
    G4cout << " @@@ NUMBER OF ControlPointSequences " << controlPSeq.getNumberOfItems() << " = " 
           << fint << G4endl;
    controlPSeq.gotoFirstItem();
    DicomBeamControlPoint* dbcp0 = 0; 
    // Only first ControlPoint has some info if it does not change 
    for( size_t i2 = 0; i2 < controlPSeq.getNumberOfItems(); i2++ ) {
      DRTControlPointSequence::Item &cpItem = controlPSeq.getCurrentItem();
      if( db->GetNControlPoints() != 0 ) dbcp0 = db->GetControlPoint(0);
      DicomBeamControlPoint* dbcp = new DicomBeamControlPoint( cpItem, dbcp0 );
      db->AddControlPoint( dbcp );
      controlPSeq.gotoNextItem();
  
    }
    
    beamSeq.gotoNextItem();
  }

  //(300a,0180) SQ (Sequence with explicit length #=1)      #  30, 1 PatientSetupSequence
  //(300c,0060) SQ (Sequence with explicit length #=1)      # 112, 1 ReferencedStructureSetSequence
  //(300e,0002) CS [UNAPPROVED]                             #  10, 1 ApprovalStatus
  //(7fe0,0010) OW (no value available)                     #   0, 1 PixelData

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePlan::CheckData0(OFString title, Sint32 val )
{
  if( val != 0 ){
    G4Exception("DicomFilePlan::CheckData",
                "DFP003",
                FatalException,
                (title + " exists, and code is not ready ").c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePlan::SetControlPointMetersets()
{
  for( size_t ii = 0; ii < theBeams.size(); ii++ ){
    theBeams[ii]->SetControlPointMetersets();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomFilePlan::DumpToFile()
{
  for( size_t ii = 0; ii < theBeams.size(); ii++ ){
    theBeams[ii]->DumpToFile();
  }
}
