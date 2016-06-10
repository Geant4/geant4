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
/*
 * G4ITNavigatorState2.cc
 *
 *  Created on: 25 févr. 2014
 *      Author: kara
 */

#include "G4ITNavigator.hh"
//#include "G4ITNavigator2.hh"

// !>

G4ITNavigator2 ::G4NavigatorState::G4NavigatorState() :
    G4ITNavigatorState_Lock2()
{
  ResetState();
}

void G4ITNavigator2::G4NavigatorState::ResetStack()
{
  fHistory.Reset();
}

void G4ITNavigator2::G4NavigatorState::ResetState()
{
  fCalculatedExitNormal = false;
  fChangedGrandMotherRefFrame = false;
  fLastTriedStepComputation = false;
  fWasLimitedByGeometry = false;
  fEntering = false;
  fExiting = false;
  fLocatedOnEdge = false;
  fLastStepWasZero = false;
  fEnteredDaughter = false;
  fExitedMother = false;
  fPushed = false;

  fValidExitNormal = false;
  fExitNormal = G4ThreeVector(0, 0, 0);

  fPreviousSftOrigin = G4ThreeVector(0, 0, 0);
  fPreviousSafety = 0.0;

  fNumberZeroSteps = 0;

  fStepEndPoint = G4ThreeVector(kInfinity, kInfinity, kInfinity);
  fLastStepEndPointLocal = G4ThreeVector(kInfinity, kInfinity, kInfinity);

  fBlockedPhysicalVolume = 0;
  fBlockedReplicaNo = -1;

  fLastLocatedPointLocal = G4ThreeVector(kInfinity, -kInfinity, 0.0);
  fLocatedOutsideWorld = false;
}

G4ITNavigator2 ::G4NavigatorState::G4NavigatorState(const G4NavigatorState& rhs) :
    G4ITNavigatorState_Lock2()
{
  fExitNormal = rhs.fExitNormal;
  fValidExitNormal = rhs.fValidExitNormal;
  fExiting = rhs.fExiting;
  fEntering = rhs.fEntering;

  fBlockedPhysicalVolume = rhs.fBlockedPhysicalVolume;
  fBlockedReplicaNo = rhs.fBlockedReplicaNo,

  fLastStepWasZero = rhs.fLastStepWasZero;

  fLocatedOutsideWorld = rhs.fLocatedOutsideWorld;
  fLastLocatedPointLocal = rhs.fLastLocatedPointLocal;
  fEnteredDaughter = rhs.fEnteredDaughter;
  fExitedMother = rhs.fExitedMother;
  fWasLimitedByGeometry = rhs.fWasLimitedByGeometry;

  fPreviousSftOrigin = rhs.fPreviousSftOrigin;
  fPreviousSafety = rhs.fPreviousSafety;

  fLastTriedStepComputation = rhs.fLastTriedStepComputation;
  fChangedGrandMotherRefFrame = rhs.fChangedGrandMotherRefFrame;
  fCalculatedExitNormal = rhs.fCalculatedExitNormal;

  fNumberZeroSteps = rhs.fNumberZeroSteps;
  fLocatedOnEdge = rhs.fLocatedOnEdge;
  fPushed = rhs.fPushed;
  fNumberZeroSteps = rhs.fNumberZeroSteps;
}

G4ITNavigator2::G4NavigatorState&
G4ITNavigator2::G4NavigatorState::operator=(const G4NavigatorState& rhs)
{
  if (this == &rhs) return *this;
  fExitNormal = rhs.fExitNormal;
  fValidExitNormal = rhs.fValidExitNormal;
  fExiting = rhs.fExiting;
  fEntering = rhs.fEntering;

  fBlockedPhysicalVolume = rhs.fBlockedPhysicalVolume;
  fBlockedReplicaNo = rhs.fBlockedReplicaNo;
  fCalculatedExitNormal = rhs.fCalculatedExitNormal;

  fLastStepWasZero = rhs.fLastStepWasZero;
  fLastTriedStepComputation = rhs.fLastTriedStepComputation;
  fChangedGrandMotherRefFrame = rhs.fChangedGrandMotherRefFrame;

  fPreviousSftOrigin = rhs.fPreviousSftOrigin;
  fPreviousSafety = rhs.fPreviousSafety;
  fNumberZeroSteps = rhs.fNumberZeroSteps;
  fLocatedOnEdge = rhs.fLocatedOnEdge;
  fWasLimitedByGeometry = rhs.fWasLimitedByGeometry;
  fPushed = rhs.fPushed;
  fNumberZeroSteps = rhs.fNumberZeroSteps;
  fEnteredDaughter = rhs.fEnteredDaughter;
  fExitedMother = rhs.fExitedMother;

  fLastLocatedPointLocal = rhs.fLastLocatedPointLocal;
  fLocatedOutsideWorld = rhs.fLocatedOutsideWorld;

  return *this;
}

G4ITNavigator2 ::G4SaveNavigatorState::G4SaveNavigatorState()
{

  sWasLimitedByGeometry = false;
  sEntering = false;
  sExiting = false;
  sLastStepWasZero = false;
  sEnteredDaughter = false;
  sExitedMother = false;


  sValidExitNormal = false;
  sExitNormal = G4ThreeVector(0, 0, 0);

  sPreviousSftOrigin = G4ThreeVector(0, 0, 0);
  sPreviousSafety = 0.0;

//  sLocatedOnEdge = false;
//  sPushed = false;
//  sNumberZeroSteps = 0;

  spBlockedPhysicalVolume = 0;
  sBlockedReplicaNo = -1;

  sLastLocatedPointLocal = G4ThreeVector(kInfinity, -kInfinity, 0.0);
  sLocatedOutsideWorld = false;
}

G4ITNavigator2 ::G4SaveNavigatorState::G4SaveNavigatorState(G4NavigatorState* rhs)
{
  sExitNormal = rhs->fExitNormal;
  sValidExitNormal = rhs->fValidExitNormal;
  sExiting = rhs->fExiting;
  sEntering = rhs->fEntering;

  spBlockedPhysicalVolume = rhs->fBlockedPhysicalVolume;
  sBlockedReplicaNo = rhs->fBlockedReplicaNo;

  sLastStepWasZero = rhs->fLastStepWasZero;

  sPreviousSftOrigin = rhs->fPreviousSftOrigin;
  sPreviousSafety = rhs->fPreviousSafety;

  sWasLimitedByGeometry = rhs->fWasLimitedByGeometry;

//  sLocatedOnEdge = rhs->fLocatedOnEdge;
//  sPushed = rhs->fPushed;
//  sNumberZeroSteps = rhs->fNumberZeroSteps;

  sEnteredDaughter = rhs->fEnteredDaughter;
  sExitedMother = rhs->fExitedMother;

  sLastLocatedPointLocal = rhs->fLastLocatedPointLocal;
  sLocatedOutsideWorld = rhs->fLocatedOutsideWorld;
}

G4ITNavigator2::G4NavigatorState& G4ITNavigator2::G4NavigatorState::operator=(const G4SaveNavigatorState& rhs)
{
  fExitNormal = rhs.sExitNormal;
  fValidExitNormal = rhs.sValidExitNormal;
  fExiting = rhs.sExiting;
  fEntering = rhs.sEntering;

  fBlockedPhysicalVolume = rhs.spBlockedPhysicalVolume;
  fBlockedReplicaNo = rhs.sBlockedReplicaNo;
//	fCalculatedExitNormal = rhs.sCalculatedExitNormal;

  fLastStepWasZero = rhs.sLastStepWasZero;
//	fLastTriedStepComputation =rhs.sLastTriedStepComputation;
//	fChangedGrandMotherRefFrame = rhs.sChangedGrandMotherRefFrame;

  fPreviousSftOrigin = rhs.sPreviousSftOrigin;
  fPreviousSafety = rhs.sPreviousSafety;
//  fNumberZeroSteps = rhs.sNumberZeroSteps;
//  fLocatedOnEdge = rhs.sLocatedOnEdge;
//  fPushed = rhs.sPushed;
  fWasLimitedByGeometry = rhs.sWasLimitedByGeometry;

  fEnteredDaughter = rhs.sEnteredDaughter;
  fExitedMother = rhs.sExitedMother;

  fLastLocatedPointLocal = rhs.sLastLocatedPointLocal;
  fLocatedOutsideWorld = rhs.sLocatedOutsideWorld;
  return *this;
}
