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
//
//---------------------------------------------------------------
//
//  G4FastTrack.cc
//
//  Description:
//    Keeps the current track information and special features
//    for Parameterised Simulation Models.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------

#include "G4ios.hh"
#include "G4FastTrack.hh"
#include "G4TransportationManager.hh"
#include "G4TouchableHistoryHandle.hh"

// -----------
// Constructor
// -----------
//
G4FastTrack::G4FastTrack(G4Envelope *anEnvelope, G4bool IsUnique)
  : fTrack                      ( nullptr    ),
    fAffineTransformationDefined( false      ),
    fEnvelope                   ( anEnvelope ),
    fIsUnique                   ( IsUnique   ),
    fEnvelopeLogicalVolume      ( nullptr    ),
    fEnvelopePhysicalVolume     ( nullptr    ),
    fEnvelopeSolid              ( nullptr    )
{}

// -----------
// Destructor:
// -----------
G4FastTrack::~G4FastTrack()
{}

//------------------------------------------------------------
// The parameterised simulation manager uses the SetCurrentTrack
// method to setup the current G4FastTrack object 
//------------------------------------------------------------
void G4FastTrack::SetCurrentTrack(const G4Track& track,
                                  const G4Navigator* theNavigator) 
{

  // -- Register track pointer (used everywhere):
  fTrack = &track;

  //-----------------------------------------------------
  // First time the track enters the volume or if the
  // Logical Volume was placed n-Times in the geometry :
  // 
  // Records the Rotation+Translation for the Envelope !
  // When the particle is inside or on the boundary, the 
  // NavigationHistory IS UP TO DATE.
  //------------------------------------------------------
  if (!fAffineTransformationDefined || !fIsUnique) FRecordsAffineTransformation(theNavigator);
  
  //-------------------------------------------
  // Records local position/momentum/direction
  // of the Track.
  // They are accessible to the user through a
  // set of Get functions and should be useful
  // to decide to trigger or not.
  //-------------------------------------------
  // -- local position:
  fLocalTrackPosition = fAffineTransformation.TransformPoint(fTrack->GetPosition());
  // -- local momentum:
  fLocalTrackMomentum = fAffineTransformation.TransformAxis(fTrack->GetMomentum());
  // -- local direction:
  fLocalTrackDirection = fLocalTrackMomentum.unit();
  // -- local polarization:
  fLocalTrackPolarization = fAffineTransformation.TransformAxis(fTrack->GetPolarization());
}

//------------------------------------
//
// 3D transformation of the envelope
// This is Done only one time.
//
//------------------------------------
void 
G4FastTrack::FRecordsAffineTransformation(const G4Navigator* theNavigator)
{

  //--------------------------------------------------------
  // Get the touchable history which represents the current
  // volume hierachy the particle is in.
  // Note that TouchableHistory allocated by the Navigator
  // must be deleted by G4FastTrack.
  //--------------------------------------------------------
  const G4Navigator* NavigatorToUse;
  if(theNavigator != nullptr ) NavigatorToUse = theNavigator;
  else NavigatorToUse = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  G4TouchableHistoryHandle history = NavigatorToUse->CreateTouchableHistoryHandle();
  
  //-----------------------------------------------------
  // Run accross the hierarchy to find the physical volume
  // associated with the envelope
  //-----------------------------------------------------
  G4int depth = (G4int)history->GetHistory()->GetDepth();
  G4int idepth;
  G4bool Done = false;
  for (idepth = 0; idepth <= depth; ++idepth)
  {
    G4VPhysicalVolume* currPV = history->GetHistory()->GetVolume(idepth);
    G4LogicalVolume* currLV   = currPV->GetLogicalVolume();
    if ( (currLV->GetRegion() == fEnvelope) && (currLV->IsRootRegion()) )
    {
      fEnvelopePhysicalVolume = currPV;
      fEnvelopeLogicalVolume  = currLV;
      fEnvelopeSolid          = currLV->GetSolid();
      Done = true;
      break;
    }
  }
  //---------------------------------------------
  //-- Verification: should be removed in future:
  //---------------------------------------------
  if ( Done == false )
    {
      G4ExceptionDescription ed;
      ed << "Can't find transformation for `" << fEnvelopePhysicalVolume->GetName() << "'" << G4endl;
      G4Exception("G4FastTrack::FRecordsAffineTransformation()",
		  "FastSim011",
		  JustWarning, ed);
    }
  else
    {
      //-------------------------------------------------------
      // Records the transformation and inverse transformation:
      //-------------------------------------------------------
      fAffineTransformation = history->GetHistory()->GetTransform(idepth);
      fInverseAffineTransformation = fAffineTransformation.Inverse();
      
      fAffineTransformationDefined = true;
    }
}
