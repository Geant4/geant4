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
// $Id: G4ErrorMagFieldLimitProcess.cc 66892 2013-01-17 10:57:59Z gunter $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorMagFieldLimitProcess.hh"
#include "G4ErrorMessenger.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#include "G4Track.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh"
#endif

//------------------------------------------------------------------------
G4ErrorMagFieldLimitProcess::
G4ErrorMagFieldLimitProcess(const G4String& processName)
  : G4VErrorLimitProcess(processName) 
{
  theStepLimit = kInfinity;
}


//------------------------------------------------------------------------
G4ErrorMagFieldLimitProcess::~G4ErrorMagFieldLimitProcess()
{ }


//------------------------------------------------------------------------
G4double G4ErrorMagFieldLimitProcess::
PostStepGetPhysicalInteractionLength( const G4Track& aTrack, G4double ,
                                            G4ForceCondition* condition )
{
  *condition = NotForced;
  const G4Field* field =
    G4TransportationManager::GetTransportationManager()->GetFieldManager()
    ->GetDetectorField();

  theStepLength = kInfinity;
  if( field != 0 ) {
    G4ThreeVector trkPosi = aTrack.GetPosition();
    G4double pos1[3];
       pos1[0] = trkPosi.x(); pos1[1] = trkPosi.y(); pos1[2] = trkPosi.z();

    G4double h1[3];
    field->GetFieldValue( pos1, h1 );

    G4ThreeVector BVec(h1[0],h1[1],h1[2]);
    G4double pmag = aTrack.GetMomentum().mag();
    G4double BPerpMom = BVec.cross( G4ThreeVector(pmag,0.,0.) ).mag() / pmag;

    theStepLength = theStepLimit * pmag / BPerpMom; 
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 ) { 
    G4cout <<  "G4ErrorMagFieldLimitProcess:: stepLength "
           << theStepLength << " B " << BPerpMom << " BVec " << BVec
           << " pmag " << pmag << G4endl;
  }
#endif
  }

  return theStepLength;
}
