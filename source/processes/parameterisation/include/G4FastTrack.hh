// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastTrack.hh,v 1.4 1999-12-15 14:53:45 gunter Exp $
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

// Class Description:
//  The G4FastTrack provides you access to the current G4Track, gives simple access to envelope 
//  related features (G4LogicalVolume, G4VSolid, G4AffineTransform references between the global 
//  and the envelope local coordinates systems) and simple access to the position, momentum 
//  expressed in the envelope coordinate system. Using those quantities and the G4VSolid methods, 
//  you can for example easily check how far you are from the envelope boundary. 
//


class G4FastTrack
{
public: // without description
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
  void SetCurrentTrack(const G4Track&, const G4Navigator* a = 0);

  //------------------------------------------------------------
  // The fast simulation manager uses the OnTheBoundaryButExiting
  // method to test if the particle is leaving the envelope.
  //------------------------------------------------------------
  G4bool OnTheBoundaryButExiting() const;

  //----------------------------------
  // Informations useful to the user :
  // General public get functions.
  //----------------------------------

public: // with Description

  const G4Track* GetPrimaryTrack() const;
  // Returns the current G4Track.

  G4Envelope* GetEnvelope() const;
  // Returns the Envelope G4LogicalVolume pointer.

  G4VPhysicalVolume* GetEnvelopePhysicalVolume() const;
  // Returns the Envelope G4VPhysicalVolume pointer.

  G4VSolid* GetEnvelopeSolid() const;
  // Returns the Envelope G4VSolid pointer.

  //-----------------------------------
  // Primary track informations in the
  // Envelope coordinate system.
  //-----------------------------------

  G4ThreeVector GetPrimaryTrackLocalPosition() const;
  // Returns the particle position in envelope coordinates.

  G4ThreeVector GetPrimaryTrackLocalMomentum() const;
  // Returns the particle momentum in envelope coordinates.

  G4ThreeVector GetPrimaryTrackLocalDirection() const;
  // Returns the particle direction in envelope coordinates.

  G4ThreeVector GetPrimaryTrackLocalPolarization() const;
  // Returns the particle polarization in envelope coordinates.
  
  //------------------------------------
  // 3D transformation of the envelope:
  //------------------------------------
  // Global -> Local

  const G4AffineTransform* GetAffineTransformation() const;  
  // Returns the envelope Global -> Local G4AffineTransform

  // Local -> Global
  const G4AffineTransform* GetInverseAffineTransformation() const; 
  // Returns the envelope Local -> Global G4AffineTransform

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
