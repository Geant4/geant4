// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastTrack.cc,v 1.3 1999-12-15 14:53:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//$Id:
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

// -----------
// Constructor
// -----------
//
G4FastTrack::G4FastTrack(G4Envelope *anEnvelope,
			 G4bool IsUnique) :
  fEnvelope(anEnvelope),fEnvelopeSolid(fEnvelope->GetSolid()),
  fIsUnique(IsUnique), fAffineTransformationDefined(false)
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
  if (!fAffineTransformationDefined || !fIsUnique) 
    FRecordsAffineTransformation(theNavigator);
  
  //-------------------------------------------
  // Records local position/momentum/direction
  // of the Track.
  // They are accessible to the user through a
  // set of Get functions and should be useful
  // to decide to trigger or not.
  //-------------------------------------------
  // -- local position:
  fLocalTrackPosition = fAffineTransformation.
    TransformPoint(fTrack->GetPosition());
  // -- local momentum:
  fLocalTrackMomentum = fAffineTransformation.
    TransformAxis(fTrack->GetMomentum());
  // -- local direction:
  fLocalTrackDirection = fLocalTrackMomentum.unit();
  // -- local polarization:
  fLocalTrackPolarization = fAffineTransformation.
    TransformAxis(fTrack->GetPolarization());
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
  if(theNavigator != 0 ) NavigatorToUse=theNavigator;
  else
    NavigatorToUse=
      G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking();
  
  G4TouchableHistory *history =  
    NavigatorToUse->CreateTouchableHistory();
  
  //-----------------------------------------------------
  // Run accross the hierarchy to find the physical volume
  // associated with the envelope
  //-----------------------------------------------------
  int depth = history->GetHistory()->GetDepth();
  int idepth, Done = 0;
  for (idepth = 0; idepth <= depth; idepth++) {
    if (history->GetHistory()->GetVolume(idepth)->GetLogicalVolume() == 
	fEnvelope) {
      fEnvelopePhysicalVolume=history->GetHistory()->GetVolume(idepth);
      Done = 1;
      break;
    }
  }
  //---------------------------------------------
  //-- Verification: should be removed in future:
  //---------------------------------------------
  if ( !Done )
    {
      G4cout << "\n\nERROR !!! can't find Transform for " <<
	fEnvelopePhysicalVolume->GetName() << "\n\n" << G4endl;
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
  
  //------------------------------------------------------
  // Delete the TouchableHistory created by the Navigator:
  //------------------------------------------------------
  delete history;
  history = 0;
}


