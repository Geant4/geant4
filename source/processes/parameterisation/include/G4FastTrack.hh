// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastTrack.hh,v 1.1 1999-01-07 16:14:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $Id:
//---------------------------------------------------------------
//
//  G4FastTrack.hh
//
//  Description:
//    Keeps the current track information and special features
//    for Parameterised Simulation Models.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------


#ifndef G4FastTrack_h
#define G4FastTrack_h

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4AffineTransform.hh"
#include "G4Track.hh"
#include "G4Navigator.hh"

//---------------------------
// For possible future needs:
//---------------------------
typedef G4LogicalVolume G4Envelope;


//-------------------------------------------
//
//        G4FastTrack class
//
//-------------------------------------------
class G4FastTrack
{
public:
  //------------------------
  // Constructor/Destructor
  //------------------------
  // Only one Constructor. By default the envelope can
  // be placed n-Times. If the user is sure that it'll be 
  // placed just one time, the IsUnique flag should be set 
  // TRUE to avoid the G4AffineTransform re-calculations each 
  // time we reach the envelope.
  G4FastTrack(G4Envelope *anEnvelope,
	      G4bool IsUnique);
  ~G4FastTrack();

  //------------------------------------------------------------
  // The fast simulation manager uses the SetCurrentTrack
  // method to setup the current G4FastTrack object 
  //------------------------------------------------------------
  void SetCurrentTrack(const G4Track&, const G4Navigator* a = NULL);

  //------------------------------------------------------------
  // The fast simulation manager uses the OnTheBoundaryButExiting
  // method to test if the particle is leaving the envelope.
  //------------------------------------------------------------
  G4bool OnTheBoundaryButExiting() const;

  //----------------------------------
  // Informations useful to the user :
  // General public get functions.
  //----------------------------------
  const G4Track* GetPrimaryTrack() const;
  G4Envelope* GetEnvelope() const;
  G4VPhysicalVolume* GetEnvelopePhysicalVolume() const;
  G4VSolid* GetEnvelopeSolid() const;

  //-----------------------------------
  // Primary track informations in the
  // Envelope coordinate system.
  //-----------------------------------
  G4ThreeVector GetPrimaryTrackLocalPosition() const;
  G4ThreeVector GetPrimaryTrackLocalMomentum() const;
  G4ThreeVector GetPrimaryTrackLocalDirection() const;
  G4ThreeVector GetPrimaryTrackLocalPolarization() const;
  
  //------------------------------------
  // 3D transformation of the envelope:
  //------------------------------------
  // Global -> Local
  const G4AffineTransform* GetAffineTransformation() const;  
  // Local -> Global
  const G4AffineTransform* GetInverseAffineTransformation() const; 

  //-----------------
  // Private members
  //-----------------
private:

  // Current G4Track pointer
  const G4Track* fTrack;

  //------------------------------------------------
  // Records the Affine/InverseAffine transformation
  // of the envelope.
  //------------------------------------------------
  void FRecordsAffineTransformation(const G4Navigator*);
  G4bool fAffineTransformationDefined;
  G4Envelope* fEnvelope;
  G4bool fIsUnique;
  G4VPhysicalVolume* fEnvelopePhysicalVolume;
  G4VSolid* fEnvelopeSolid;
  G4ThreeVector fLocalTrackPosition, fLocalTrackMomentum, 
    fLocalTrackDirection, fLocalTrackPolarization;
  G4AffineTransform fAffineTransformation, fInverseAffineTransformation;
};

//*******************************************************************
//
//  Inline functions
//
//*******************************************************************

inline G4Envelope* G4FastTrack::GetEnvelope() const
{
  return fEnvelope;
}

inline G4VPhysicalVolume* G4FastTrack::GetEnvelopePhysicalVolume() const
{
  return fEnvelopePhysicalVolume;
}

inline G4VSolid* G4FastTrack::GetEnvelopeSolid() const
{
  return fEnvelopeSolid;
}

inline const G4Track* G4FastTrack::GetPrimaryTrack() const
{
  return fTrack;
}

inline G4ThreeVector G4FastTrack::GetPrimaryTrackLocalPosition() const
{
  return fLocalTrackPosition;
}

inline G4ThreeVector G4FastTrack::GetPrimaryTrackLocalMomentum() const
{
  return fLocalTrackMomentum;
}

inline G4ThreeVector G4FastTrack::GetPrimaryTrackLocalDirection() const
{
  return fLocalTrackDirection;
}

inline G4ThreeVector G4FastTrack::GetPrimaryTrackLocalPolarization() const
{
  return fLocalTrackPolarization;
}

inline const G4AffineTransform* G4FastTrack::GetAffineTransformation() const
{
  return &fAffineTransformation;
}

inline const G4AffineTransform* G4FastTrack::GetInverseAffineTransformation() const
{
  return &fInverseAffineTransformation;
}

inline G4bool G4FastTrack::OnTheBoundaryButExiting() const 
{
  // tests if particle are on the boundary and leaving.
  return GetEnvelopeSolid()->
    DistanceToOut(GetPrimaryTrackLocalPosition(),
		  GetPrimaryTrackLocalDirection())==0.;
}

#endif
