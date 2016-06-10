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
//
// $Id: G4ParticleChangeForMSC.cc 68795 2013-04-05 13:24:46Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	
//	
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//   Update for model variant of msc                16 Jan  2004  V.Ivanchenko
// --------------------------------------------------------------

#include "G4ParticleChangeForMSC.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

G4ParticleChangeForMSC::G4ParticleChangeForMSC()
  : G4VParticleChange()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForMSC::G4ParticleChangeForMSC() " << G4endl;
  }
#endif
}

G4ParticleChangeForMSC::~G4ParticleChangeForMSC()
{
#ifdef G4VERBOSE
  if (verboseLevel>2) {
    G4cout << "G4ParticleChangeForMSC::~G4ParticleChangeForMSC() " << G4endl;
  }
#endif
}

G4ParticleChangeForMSC::
G4ParticleChangeForMSC(const G4ParticleChangeForMSC &right)
  : G4VParticleChange(right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForMSC::  copy constructor is called " << G4endl;
   }
      theMomentumDirection = right.theMomentumDirection;
      thePosition = right.thePosition;
}

// assignment operator
G4ParticleChangeForMSC & G4ParticleChangeForMSC::operator=(
                                   const G4ParticleChangeForMSC &right)
{
   if (verboseLevel>1) {
    G4cout << "G4ParticleChangeForMSC:: assignment operator is called " << G4endl;
   }
   if (this != &right)
   {
      theListOfSecondaries = right.theListOfSecondaries;
      theSizeOftheListOfSecondaries = right.theSizeOftheListOfSecondaries;
      theNumberOfSecondaries = right.theNumberOfSecondaries;
      theStatusChange = right.theStatusChange;
      theLocalEnergyDeposit = right.theLocalEnergyDeposit;
      theSteppingControlFlag = right.theSteppingControlFlag;
      theTrueStepLength = right.theTrueStepLength;

      theMomentumDirection = right.theMomentumDirection;
      thePosition = right.thePosition;
   }
   return *this;
}

//----------------------------------------------------------------
// methods for updating G4Step
//

G4Step* G4ParticleChangeForMSC::UpdateStepForAlongStep(G4Step* pStep)
{
  //  Update the G4Step specific attributes
  pStep->SetStepLength(theTrueStepLength);
  theStatusChange = pStep->GetTrack()->GetTrackStatus();

  // Multiple scattering calculates the final state of the particle
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint->SetPosition( thePosition );
  return pStep;
}

G4Step* G4ParticleChangeForMSC::UpdateStepForPostStep(G4Step* pStep)
{
  // Multiple scattering calculates the final state of the particle
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint->SetPosition( thePosition );

  //  Update the G4Step specific attributes
  return pStep;
}


//----------------------------------------------------------------
// methods for printing messages
//

void G4ParticleChangeForMSC::DumpInfo() const
{
// use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4int oldprc = G4cout.precision(3);
  G4cout << "        Position - x (mm)   : "
       << std::setw(20) << thePosition.x()/mm
       << G4endl;
  G4cout << "        Position - y (mm)   : "
       << std::setw(20) << thePosition.y()/mm
       << G4endl;
  G4cout << "        Position - z (mm)   : "
       << std::setw(20) << thePosition.z()/mm
       << G4endl;
  G4cout << "        Momentum Direct - x : "
       << std::setw(20) << theMomentumDirection.x()
       << G4endl;
  G4cout << "        Momentum Direct - y : "
       << std::setw(20) << theMomentumDirection.y()
       << G4endl;
  G4cout << "        Momentum Direct - z : "
       << std::setw(20) << theMomentumDirection.z()
       << G4endl;
  G4cout.precision(oldprc);
}


G4bool G4ParticleChangeForMSC::CheckIt(const G4Track& aTrack)
{
  G4bool    itsOK = true;
  G4bool    exitWithError = false;

  G4double  accuracy;

  // check

  // MomentumDirection should be unit vector
  accuracy = std::fabs(theMomentumDirection.mag2()-1.0);
  if (accuracy > accuracyForWarning) {
    itsOK = false;
    exitWithError = (accuracy > accuracyForException);
#ifdef G4VERBOSE
    G4cout << "  G4ParticleChangeForMSC::CheckIt  : ";
    G4cout << "the Momentum Change is not unit vector !!"
	   << "  Difference:  " << accuracy << G4endl;
    G4cout << aTrack.GetDefinition()->GetParticleName()
	   << " E=" << aTrack.GetKineticEnergy()/MeV
	   << " pos=" << aTrack.GetPosition().x()/m
	   << ", " << aTrack.GetPosition().y()/m
	   << ", " << aTrack.GetPosition().z()/m
	   <<G4endl;
#endif
  }

  // dump out information of this particle change
#ifdef G4VERBOSE
  if (!itsOK) DumpInfo();
#endif

  // Exit with error
  if (exitWithError) {
    G4Exception("G4ParticleChangeForMSC::CheckIt",
		"300",
		EventMustBeAborted,
		"momentum direction was illegal");
  }
  //correction
  if (!itsOK) {
    G4double vmag = theMomentumDirection.mag();
    theMomentumDirection = (1./vmag)*theMomentumDirection;
  }

  itsOK = (itsOK) && G4VParticleChange::CheckIt(aTrack);
  return itsOK;
}




