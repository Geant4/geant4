//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Track.cc,v 1.15 2001-12-10 08:36:54 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4Track.cc
//
//---------------------------------------------------------------
//   Add copy constructor            Hisaya Feb. 07 01
//   Fix GetVelocity                 Hisaya Feb. 17 01
//   Modification for G4TouchableHandle             22 Oct. 2001  R.Chytracek//
#include "G4Track.hh"

G4Allocator<G4Track> aTrackAllocator;

///////////////////////////////////////////////////////////
G4Track::G4Track(G4DynamicParticle* apValueDynamicParticle,
                 G4double aValueTime,
                 const G4ThreeVector& aValuePosition)
///////////////////////////////////////////////////////////
  : fCurrentStepNumber(0),    fPosition(aValuePosition),
    fGlobalTime(aValueTime),  fLocalTime(0.),
    fTrackLength(0.),
    fParentID(0),             fTrackID(0),
    fpDynamicParticle(apValueDynamicParticle),
    fTrackStatus(fAlive),
    fBelowThreshold(false),   fGoodForTracking(false),
    fWeight(1.0),
    fpStep(0),
    fpLVAtVertex(0),          fpCreatorProcess(0),
    fpUserInformation(0)
{    
}

//////////////////
G4Track::G4Track()
//////////////////
  : fCurrentStepNumber(0),    
    fGlobalTime(0),           fLocalTime(0.),
    fTrackLength(0.),
    fParentID(0),             fTrackID(0),
    fpDynamicParticle(0),
    fTrackStatus(fAlive),
    fBelowThreshold(false),   fGoodForTracking(false),
    fWeight(1.0),
    fpStep(0),
    fpLVAtVertex(0),          fpCreatorProcess(0),
    fpUserInformation(0)
{
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
   fPosition = right.fPosition;
   fGlobalTime = right.fGlobalTime;
   fLocalTime = right.fLocalTime;
   fTrackLength = right.fTrackLength;
   fWeight = right.fWeight;

   // Track ID (and Parent ID) is not copied and set to zero for new track
   fTrackID = 0;
   fParentID =0;

   // CurrentStepNumber is set to be 0
   fCurrentStepNumber = 0;

   // dynamic particle information 
   fpDynamicParticle = new G4DynamicParticle(*(right.fpDynamicParticle));
 
   // track status and flags for tracking  
   fTrackStatus = right.fTrackStatus;
   fBelowThreshold = right.fBelowThreshold;
   fGoodForTracking = right.fGoodForTracking;
   
   // Step information (Step Length, Step Number, pointer to the Step,) 
   // are not copied
   fpStep=0;

   // vertex information
   fVtxPosition = right.fVtxPosition;
   fpLVAtVertex = right.fpLVAtVertex;
   fVtxKineticEnergy = right.fVtxKineticEnergy;
   fVtxMomentumDirection = right.fVtxMomentumDirection;

   // CreatorProcess is not copied 
   fpCreatorProcess = 0;
    
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
 



