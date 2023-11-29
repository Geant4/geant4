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
// G4ParticleChangeForMSC class implementation
//
// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 16 January 2004
//                                 23 August 2022
// --------------------------------------------------------------------

#include "G4ParticleChangeForMSC.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ExceptionSeverity.hh"

// --------------------------------------------------------------------
G4ParticleChangeForMSC::G4ParticleChangeForMSC() 
{
  // Disable flag that is enabled in G4VParticleChange if G4VERBOSE.
  debugFlag = false;
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForMSC::UpdateStepForAlongStep(G4Step* pStep)
{
  // update the G4Step specific attributes
  pStep->SetStepLength(theTrueStepLength);

  // multiple scattering calculates the final state of the particle
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint->SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint->SetPosition(thePosition);
  return pStep;
}

