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
// $Id: G4VErrorLimitProcess.cc 66892 2013-01-17 10:57:59Z gunter $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4VErrorLimitProcess.hh"
#include "G4ErrorMessenger.hh"

//------------------------------------------------------------------------
G4VErrorLimitProcess::G4VErrorLimitProcess(const G4String& processName)
  : G4VDiscreteProcess (processName) 
{
  theStepLimit = kInfinity;
  theStepLength = kInfinity;
}


//------------------------------------------------------------------------
G4VErrorLimitProcess::~G4VErrorLimitProcess()
{
}


//------------------------------------------------------------------------
G4double G4VErrorLimitProcess::
GetMeanFreePath(const class G4Track &, G4double, enum G4ForceCondition *)
{
  return theStepLength;
}


//------------------------------------------------------------------------
G4VParticleChange* G4VErrorLimitProcess::
PostStepDoIt(const G4Track& aTrack, const G4Step& )
{
  theParticleChange.Initialize(aTrack);
  return &theParticleChange;
}
