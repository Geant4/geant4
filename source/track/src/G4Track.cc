// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Track.cc,v 1.11 2001-02-17 10:51:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4Track.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------
//   Add copy constructor            Hisaya Feb. 07 01
//   Fix GetVelocity                 Hisaya Feb. 17 01
//
#include "G4Track.hh"

G4Allocator<G4Track> aTrackAllocator;

///////////////////////////////////////////////////////////
G4Track::G4Track(G4DynamicParticle* apValueDynamicParticle,
                 G4double aValueTime,
                 const G4ThreeVector& aValuePosition)
///////////////////////////////////////////////////////////
{
   fpDynamicParticle = apValueDynamicParticle;
   fCurrentStepNumber = 0;
   fGlobalTime = aValueTime;
   fLocalTime = 0.;
   fTrackLength = 0.;
   fPosition = aValuePosition;
   fpTouchable = 0;
   fpNextTouchable = 0; 
   fpStep=0;

   fpLVAtVertex = 0;
   fpCreatorProcess = 0;

   fTrackStatus = fAlive;

   fBelowThreshold = false;
   fGoodForTracking = false;
   fWeight = 1.0;

   fpUserInformation = 0;
}

//////////////////
G4Track::G4Track()
//////////////////
{
   fCurrentStepNumber = 0;
   fGlobalTime = 0.;
   fLocalTime = 0.;
   fTrackLength = 0.;
   fParentID = 0;
   fTrackID = 0;
   fpTouchable = 0;
   fpNextTouchable = 0;
   fpStep=0;

   fpDynamicParticle = 0;
   fpLVAtVertex = 0;
   fpCreatorProcess = 0;

   fTrackStatus = fAlive;
   fBelowThreshold = false;
   fGoodForTracking = false;
   fWeight = 1.0;

   fpUserInformation = 0;
}
//////////////////
G4Track::G4Track(const G4Track& right)
//////////////////
{
  *this = right;
}

///////////////////
G4Track::~G4Track()
///////////////////
{
   delete fpDynamicParticle;
   delete fpUserInformation;
}

//////////////////
G4Track & G4Track::operator=(const G4Track &right)
//////////////////
{
  if (this != &right) {
   fCurrentStepNumber = right.fCurrentStepNumber;
   fGlobalTime = right.fGlobalTime;
   fLocalTime = right.fLocalTime;
   fTrackLength = right.fTrackLength;
   fParentID = right.fParentID;

   // Track ID is not copied and set to zero for new track
   fTrackID = 0;

   // pointers to Touchables or Step are not copied
   fpTouchable = 0;
   fpNextTouchable = 0;
   fpStep=0;

   fpDynamicParticle = new G4DynamicParticle(*(right.fpDynamicParticle));

   fVtxKineticEnergy = right.fVtxKineticEnergy;
   fVtxPosition = right.fVtxPosition;
   fpLVAtVertex = right.fpLVAtVertex;
   fpCreatorProcess = right.fpCreatorProcess;

   fTrackStatus = right.fTrackStatus;
   fBelowThreshold = right.fBelowThreshold;
   fGoodForTracking = right.fGoodForTracking;
   fWeight = right.fWeight;

   fpUserInformation = right.fpUserInformation;
  }
  return *this;
}

///////////////////
G4double G4Track::GetVelocity() const
///////////////////
{ 
  G4double velocity ;
  
  G4double mass = fpDynamicParticle->GetMass();

  // mass less particle  
  if( mass == 0. ){
    velocity = c_light ; 

    // special case for photons
    if(fpDynamicParticle->GetDefinition()->GetParticleName()=="opticalphoton"){
      G4Material*
	mat=fpTouchable->GetVolume()->GetLogicalVolume()->GetMaterial();
 
      if(mat->GetMaterialPropertiesTable() != 0){
	if(mat->GetMaterialPropertiesTable()->GetProperty("RINDEX") != 0 ){ 
          // light velocity = c/reflection-index 
	  velocity /= 
	    mat->GetMaterialPropertiesTable()->GetProperty("RINDEX")->
	    GetProperty(fpDynamicParticle->GetTotalMomentum()) ; 
	}
      }  
    }
  } else {
    G4double T = fpDynamicParticle->GetKineticEnergy();
    velocity = c_light*sqrt(T*(T+2.*mass))/(T+mass) ;
  }
  
  return velocity ; 
}
 



