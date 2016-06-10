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
// $Id: G4ErrorStepLengthLimitProcess.cc 66892 2013-01-17 10:57:59Z gunter $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorStepLengthLimitProcess.hh"
#include "G4ErrorMessenger.hh"
#include "G4SystemOfUnits.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh"
#endif

//------------------------------------------------------------------------
G4ErrorStepLengthLimitProcess::
G4ErrorStepLengthLimitProcess(const G4String& processName)
  : G4VErrorLimitProcess(processName) 
{
  theStepLimit = 1000.*mm; // kInfinity;
}


//------------------------------------------------------------------------
G4ErrorStepLengthLimitProcess::~G4ErrorStepLengthLimitProcess()
{ }


//------------------------------------------------------------------------
G4double G4ErrorStepLengthLimitProcess::
PostStepGetPhysicalInteractionLength( const G4Track&, G4double,
                                      G4ForceCondition* condition )
{
  *condition = NotForced;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 )
  { 
     G4cout << "G4ErrorStepLengthLimitProcess::PostStepGetPhysicalInteractionLength "
            << theStepLimit << G4endl;
  }
#endif

  return theStepLimit;
}
