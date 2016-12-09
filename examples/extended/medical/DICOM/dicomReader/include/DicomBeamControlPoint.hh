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
#ifndef DicomBeamControlPoint__HH
#define DicomBeamControlPoint__HH

#include <vector>
#include <iostream>
#include "dcmtk/dcmrt/seq/drtcps.h"    // for ControlPointSequence
class DicomBeamDevicePos;

#include "G4ThreeVector.hh"

class DicomBeamControlPoint 
{ 
public:
  DicomBeamControlPoint(DRTControlPointSequence::Item cpItem, DicomBeamControlPoint* point0 );
  ~DicomBeamControlPoint(){};

public:
  void SetIndex( Sint32 dat ) {
    theIndex = dat;
  }
  Sint32 GetIndex() const {
    return theIndex;
  }
  void SetNominalBeamEnergy(Float64 dat){
    theNominalBeamEnergy = dat;
  }
  Float64 GetNominalBeamEnergy() const {
    return theNominalBeamEnergy;
  }
  void SetGantryAngle(Float64 dat){
    theGantryAngle = dat;
  }
  void SetGantryRotationDirection(OFString dat){
    theGantryRotationDirection = dat;
  }
  void SetBeamLimitingDeviceAngle(Float64 dat){
    theBeamLimitingDeviceAngle = dat;
  }
  void SetBeamLimitingDeviceRotationDirection(OFString dat){
    theBeamLimitingDeviceRotationDirection = dat;
  }
  void SetPatientSupportAngle(Float64 dat){
    thePatientSupportAngle = dat;
  }
  void SetPatientSupportRotationDirection(OFString dat){
    thePatientSupportRotationDirection = dat;
  }
  void SetTableTopEccentricAngle(Float64 dat){
    theTableTopEccentricAngle = dat;
  }
  void SetTableTopEccentricRotationDirection(OFString dat){
    theTableTopEccentricRotationDirection = dat;
  }
  void SetIsocenterPosition(G4ThreeVector dat){
    theIsocenterPosition = dat;
  }
  void SetSourceToSurfaceDistance(Float64 dat){
    theSourceToSurfaceDistance = dat;
  }
  void SetCumulativeMetersetWeight(Float64 dat){
    theCumulativeMetersetWeight = dat;
  }
  void SetMetersetWeight(Float64 dat){
    theMetersetWeight = dat;
  }
  void SetGantryPitchAngle(Float32 dat){
    theGantryPitchAngle = dat;
  }
  void SetSurfaceEntryPoint(Float64 dat){
    theSurfaceEntryPoint = dat;
  }
  void SetTableTopEccentricAxisDistance(Float64 dat){
    theTableTopEccentricAxisDistance = dat;
  }
  void SetTableTopLateralPosition(Float64 dat){
    theTableTopLateralPosition = dat;
  }
  void SetTableTopLongitudinalPosition(Float64 dat){
    theTableTopLongitudinalPosition = dat;
  }
  void SetTableTopPitchAngle(Float32 dat){
    theTableTopPitchAngle = dat;
  }
  void SetTableTopPitchRotationDirection(OFString dat){
    theTableTopPitchRotationDirection = dat;
  }
  void SetTableTopRollAngle(Float32 dat){
    theTableTopRollAngle = dat;
  }
  void SetTableTopRollRotationDirection(OFString dat){
    theTableTopRollRotationDirection = dat;
  }
  void SetTableTopVerticalPosition(Float64 dat){
    theTableTopVerticalPosition = dat;
  }
  OFString GetGantryRotationDirection() const {
    return theGantryRotationDirection;
  }
  Float64 GetBeamLimitingDeviceAngle() const { 
    return theBeamLimitingDeviceAngle;
  }
  OFString GetBeamLimitingDeviceRotationDirection() const {
    return theBeamLimitingDeviceRotationDirection;
  }
  Float64 GetPatientSupportAngle() const {
    return thePatientSupportAngle;
  }
  OFString GetPatientSupportRotationDirection() const {
    return thePatientSupportRotationDirection;
  }
  Float64 GetTableTopEccentricAngle() const {
    return theTableTopEccentricAngle;
  }
  OFString GetTableTopEccentricRotationDirection() const {
    return theTableTopEccentricRotationDirection;
  }
  G4ThreeVector GetIsocenterPosition() const {
    return theIsocenterPosition;
  }
  Float64 GetSourceToSurfaceDistance() const {
    return theSourceToSurfaceDistance;
  }
  Float64 GetCumulativeMetersetWeight() const {
    return theCumulativeMetersetWeight;
  }   
  Float64 GetMetersetWeight() const {
    return theMetersetWeight;
  }   
  Float32 GetGantryPitchAngle() const {
    return theGantryPitchAngle;
  }
  Float64 GetSurfaceEntryPoint() const {
    return theSurfaceEntryPoint;
  }
  Float64 GetTableTopEccentricAxisDistance() const {
    return theTableTopEccentricAxisDistance;
  }
  Float64 GetTableTopLateralPosition() const {
    return theTableTopLateralPosition;
  }
  Float64 GetTableTopLongitudinalPosition() const {
    return theTableTopLongitudinalPosition;
  }
  Float32 GetTableTopPitchAngle() const {
    return theTableTopPitchAngle;
  }
  OFString GetTableTopPitchRotationDirection() const {
    return theTableTopPitchRotationDirection;
  }
  Float32 GetTableTopRollAngle() const {
    return theTableTopRollAngle;
  }
  OFString GetTableTopRollRotationDirection() const {
    return theTableTopRollRotationDirection;
  }
  Float64 GetTableTopVerticalPosition() const {
    return theTableTopVerticalPosition;
  }  
  
  void AddDevice(DicomBeamDevicePos* dbd){
    theDevices.push_back(dbd);
  }

  void DumpToFile( std::ofstream& out );
  
  void Print( std::ostream& out );
  
private:
  Sint32 theIndex;
  Float64 theNominalBeamEnergy;
  Float64 theGantryAngle;
  OFString theGantryRotationDirection;
  Float64 theBeamLimitingDeviceAngle; 
  OFString theBeamLimitingDeviceRotationDirection;
  Float64 thePatientSupportAngle;
  OFString thePatientSupportRotationDirection;
  Float64 theTableTopEccentricAngle;
  OFString theTableTopEccentricRotationDirection;
  G4ThreeVector theIsocenterPosition;
  Float64 theSourceToSurfaceDistance;
  Float64 theCumulativeMetersetWeight;
  Float64 theMetersetWeight;
  Float32 theGantryPitchAngle;
  Float64 theSurfaceEntryPoint;
  Float64 theTableTopEccentricAxisDistance;
  Float64 theTableTopLateralPosition;
  Float64 theTableTopLongitudinalPosition;
  Float32 theTableTopPitchAngle;
  OFString theTableTopPitchRotationDirection;
  Float32 theTableTopRollAngle;
  OFString theTableTopRollRotationDirection;
  Float64 theTableTopVerticalPosition;

  std::vector<DicomBeamDevicePos *> theDevices;

};

#endif  
