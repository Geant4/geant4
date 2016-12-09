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
#include "DicomBeamControlPoint.hh"
#include "DicomBeamDevicePos.hh"

#include "dcmtk/dcmrt/seq/drtcps.h"    // for ControlPointSequence
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeamControlPoint::DicomBeamControlPoint(DRTControlPointSequence::Item cpItem, 
 DicomBeamControlPoint* point0 )
{
  OFString fstr;
  Sint32 fint;
  Float64 ffloat;
  Float32 ffloat32;
  OFVector<Float64> fvfloat;
  
  cpItem.getControlPointIndex(fint);
  G4cout << "  @ DicomBeamControlPoint: " << fint << G4endl;
  G4cout << "   " << " ControlPointIndex " << fint << G4endl;
  SetIndex( fint );
  if( cpItem.getNominalBeamEnergy(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetNominalBeamEnergy();
  }
  G4cout << "   " << " NominalBeamEnergy " << ffloat << G4endl;
  SetNominalBeamEnergy(ffloat);
  cpItem.getDoseRateSet(ffloat); // != EC_Normal ) {
  G4cout << "   " << " DoseRateSet " << ffloat << G4endl;
  
  DRTBeamLimitingDevicePositionSequence beamLDPS = cpItem.getBeamLimitingDevicePositionSequence();
  G4cout << "  @ NUMBER OF BeamLimitingDevicePositionSequence " << beamLDPS.getNumberOfItems()
         << G4endl;
  beamLDPS.gotoFirstItem();
  for( size_t i3 = 0; i3 < beamLDPS.getNumberOfItems(); i3++ ) {
    DRTBeamLimitingDevicePositionSequence::Item bldpsItem = beamLDPS.getCurrentItem();
    DicomBeamDevicePos* dbd = new DicomBeamDevicePos(bldpsItem);
    AddDevice(dbd);
 
    beamLDPS.gotoNextItem();
  }
  
  cpItem.getGantryAngle(ffloat);
  G4cout << "   " << " GantryAngle " << ffloat << G4endl;

  cpItem.getGantryRotationDirection(fstr); //**
  G4cout << "   " << " GantryRotationDirection " << fstr << G4endl;
  if( fstr == "CC" ) { // counter-clockwise
    SetGantryAngle(-ffloat);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetGantryAngle(ffloat);
  }
  if( cpItem.getBeamLimitingDeviceAngle(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetBeamLimitingDeviceAngle();
  }
  G4cout << "   " << " BeamLimitingDeviceAngle " << ffloat << G4endl;

  if( cpItem.getBeamLimitingDeviceRotationDirection(fstr) != EC_Normal ) {
    if( point0 ) fstr = point0->GetBeamLimitingDeviceRotationDirection();
  }
  if( fstr == "CC" ) { // counter-clockwise
    SetBeamLimitingDeviceAngle(-ffloat);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetBeamLimitingDeviceAngle(ffloat);
  }   
  G4cout << "   " << " BeamLimitingDeviceRotationDirection " << fstr << G4endl;
  SetBeamLimitingDeviceRotationDirection(fstr);

  if( cpItem.getPatientSupportAngle(ffloat) != EC_Normal ) {
    if( point0 ) fstr = point0->GetPatientSupportAngle();
  }
  G4cout << "   " << " PatientSupportAngle " << ffloat << G4endl;

  if( cpItem.getPatientSupportRotationDirection(fstr) != EC_Normal ) {
    if( point0 ) fstr = point0->GetPatientSupportRotationDirection();
  }
  G4cout << "   " << " PatientSupportRotationDirection " << fstr << G4endl;
  SetPatientSupportRotationDirection(fstr);
  if( fstr == "CC" ) { // counter-clockwise
    SetPatientSupportAngle(-ffloat);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetPatientSupportAngle(ffloat);
  }

  if( cpItem.getTableTopEccentricAngle(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetTableTopEccentricAngle();
  }
  G4cout << "   " << " TableTopEccentricAngle " << ffloat << G4endl;

  if( cpItem.getTableTopEccentricRotationDirection(fstr) != EC_Normal ) {
    if( point0 ) fstr = point0->GetTableTopEccentricRotationDirection();
  }
  if( fstr == "CC" ) { // counter-clockwise
    SetTableTopEccentricAngle(-ffloat);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetTableTopEccentricAngle(ffloat);
  }
  G4cout << "   " << " TableTopEccentricRotationDirection " << fstr << G4endl;
  SetTableTopEccentricRotationDirection(fstr);

  G4ThreeVector isoCenter;
  if( cpItem.getIsocenterPosition(fvfloat) != EC_Normal ) {
    if( point0 ) isoCenter = point0->GetIsocenterPosition();
  } else {
    isoCenter = G4ThreeVector(fvfloat[0],fvfloat[1],fvfloat[2]);
  }
  G4cout << "   " << " IsocenterPosition " << isoCenter << G4endl;
  SetIsocenterPosition(isoCenter);

  if( cpItem.getSourceToSurfaceDistance(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetSourceToSurfaceDistance();
  }
  G4cout << "   " << " SourceToSurfaceDistance " << ffloat << G4endl;
  SetSourceToSurfaceDistance(ffloat);

  cpItem.getCumulativeMetersetWeight(ffloat);
  G4cout << "   " << " CumulativeMetersetWeight " << ffloat << G4endl;
  SetCumulativeMetersetWeight(ffloat);
  
  if( cpItem.getGantryPitchAngle(ffloat32) != EC_Normal ) {
    if( point0 ) ffloat32 = point0->GetGantryPitchAngle();
  }
  G4cout << "   " << " GantryPitchAngle " << ffloat32 << G4endl;
  SetGantryPitchAngle(ffloat32);

  if( cpItem.getSurfaceEntryPoint(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetSurfaceEntryPoint();
  }
  G4cout << "   " << " SurfaceEntryPoint " << ffloat << G4endl;
  SetSurfaceEntryPoint(ffloat);

  if( cpItem.getTableTopEccentricAxisDistance(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetTableTopEccentricAxisDistance();
  }
  G4cout << "   " << " TableTopEccentricAxisDistance " << ffloat << G4endl;
  SetTableTopEccentricAxisDistance(ffloat);

  if( cpItem.getTableTopLateralPosition(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetTableTopLateralPosition();
  }
  G4cout << "   " << " TableTopLateralPosition " << ffloat << G4endl;
  SetTableTopLateralPosition(ffloat);

  if( cpItem.getTableTopLongitudinalPosition(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetTableTopLongitudinalPosition();
  }
  G4cout << "   " << " TableTopLongitudinalPosition " << ffloat << G4endl;
  SetTableTopLongitudinalPosition(ffloat);

  if( cpItem.getTableTopPitchAngle(ffloat32) != EC_Normal ) {
    if( point0 ) ffloat32 = point0->GetTableTopPitchAngle();
  }
  G4cout << "   " << " TableTopPitchAngle " << ffloat32 << G4endl;

  if( cpItem.getTableTopPitchRotationDirection(fstr) != EC_Normal ) {
    if( point0 ) fstr = point0->GetTableTopPitchRotationDirection();
  }
  if( fstr == "CC" ) { // counter-clockwise
    SetTableTopPitchAngle(-ffloat32);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetTableTopPitchAngle(ffloat32);
  }
  G4cout << "   " << " TableTopPitchRotationDirection " << fstr << G4endl;
  SetTableTopPitchRotationDirection(fstr);

  if( cpItem.getTableTopRollAngle(ffloat32) != EC_Normal ) {
    if( point0 ) ffloat32 = point0->GetTableTopRollAngle();
  }
  G4cout << "   " << " TableTopRollAngle " << ffloat32 << G4endl;

  if( cpItem.getTableTopRollRotationDirection(fstr) != EC_Normal ) {
    if( point0 ) fstr = point0->GetTableTopRollRotationDirection();
  }
  if( fstr == "CC" ) { // counter-clockwise
    SetTableTopRollAngle(-ffloat32);
  } else if( fstr == "CW" || fstr == "NONE" || fstr == "") { // clockwise
    SetTableTopRollAngle(ffloat32);
  }
  G4cout << "   " << " TableTopRollRotationDirection " << fstr << G4endl;
  SetTableTopRollRotationDirection(fstr);

  if( cpItem.getTableTopVerticalPosition(ffloat) != EC_Normal ) {
    if( point0 ) ffloat = point0->GetTableTopVerticalPosition();
  }
  G4cout << "   " << " TableTopVerticalPosition " << ffloat << G4endl;
  SetTableTopVerticalPosition(ffloat);
  
  // --- get DICOM sequence attributes ---
  //t        DRTWedgePositionSequence &getWedgePositionSequence()
  //t        const DRTWedgePositionSequence &getWedgePositionSequence() const
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamControlPoint::Print( std::ostream&  )
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamControlPoint::DumpToFile( std::ofstream& fout )
{
  fout << ":P ControlPointIndex " << theIndex << G4endl;
  fout << ":P NominalBeamEnergy " << theNominalBeamEnergy << G4endl;
  fout << ":P IsocenterPositionX " << theIsocenterPosition.x() << G4endl;
  fout << ":P IsocenterPositionY " << theIsocenterPosition.y() << G4endl;
  fout << ":P IsocenterPositionZ " << theIsocenterPosition.z() << G4endl;

  //  std::string iistr = std::to_string(theControlPointIndex);
  fout << ":P SourceToSurfaceDistance " << theSourceToSurfaceDistance << G4endl;
  fout << ":P GantryAngle " << theGantryAngle << G4endl;
  fout << ":P BeamLimitingDeviceAngle " << theBeamLimitingDeviceAngle << G4endl;
  fout << ":P PatientSupportAngle " << thePatientSupportAngle << G4endl;
  fout << ":P TableTopEccentricAngle " << theTableTopEccentricAngle << G4endl;
  fout << ":P SourceToSurfaceDistance " << theSourceToSurfaceDistance<< G4endl;
  fout << ":P MetersetWeight " << theMetersetWeight<< G4endl;
  fout << ":P GantryPitchAngle " << theGantryPitchAngle << G4endl;
  fout << ":P SurfaceEntryPoint " << theSurfaceEntryPoint << G4endl;
  fout << ":P TableTopEccentricAxisDistance " << theTableTopEccentricAxisDistance << G4endl;
  fout << ":P TableTopLateralPosition " << theTableTopLateralPosition << G4endl;
  fout << ":P TableTopLongitudinalPosition " << theTableTopLongitudinalPosition<< G4endl; 
  fout << ":P TableTopPitchAngle " << theTableTopPitchAngle << G4endl;
  fout << ":P TableTopRollAngle " << theTableTopRollAngle << G4endl;
  fout << ":P TableTopVerticalPosition " << theTableTopVerticalPosition << G4endl;

  for( size_t ii = 0; ii < theDevices.size(); ii++ ){
    theDevices[ii]->DumpToFile(fout);
  }
  
}
