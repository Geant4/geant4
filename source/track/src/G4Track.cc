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
// G4Track class implementation
//
// Author: Katsuya Amako, KEK - 1996
// Revisions: Hisaya Kurashige, 1998-2011
// --------------------------------------------------------------------

#include "G4Track.hh"
#include "G4PhysicalConstants.hh"
#include "G4VAuxiliaryTrackInformation.hh"

#include <iostream>
#include <iomanip>

// --------------------------------------------------------------------
G4Allocator<G4Track>*& aTrackAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4Track>* _instance = nullptr;
  return _instance;
}

// --------------------------------------------------------------------
G4Track::G4Track(G4DynamicParticle* apValueDynamicParticle,
                 G4double aValueTime,
                 const G4ThreeVector& aValuePosition)
  : fPosition(aValuePosition)
  , fGlobalTime(aValueTime)
  , fVelocity(c_light)
{
  fpDynamicParticle = (apValueDynamicParticle) != nullptr
                    ? apValueDynamicParticle : new G4DynamicParticle();
  // check if the particle type is Optical Photon
  is_OpticalPhoton =
    (fpDynamicParticle->GetDefinition()->GetPDGEncoding() == -22);
}

// --------------------------------------------------------------------
G4Track::G4Track()
  : fVelocity(c_light)
  , fpDynamicParticle(new G4DynamicParticle())
{}

// --------------------------------------------------------------------
G4Track::G4Track(const G4Track& right)
  : fVelocity(c_light)
{
  *this = right;
}

// --------------------------------------------------------------------
G4Track::~G4Track()
{
  delete fpDynamicParticle;
  delete fpUserInformation;
  ClearAuxiliaryTrackInformation();
}

// --------------------------------------------------------------------
G4Track& G4Track::operator=(const G4Track& right)
{
  if(this != &right)
  {
    fPosition    = right.fPosition;
    fGlobalTime  = right.fGlobalTime;
    fLocalTime   = right.fLocalTime;
    fTrackLength = right.fTrackLength;
    fWeight      = right.fWeight;
    fStepLength  = right.fStepLength;

    // additional fields required for geometrical splitting
    fpTouchable = right.fpTouchable;
    fpNextTouchable = right.fpNextTouchable;
    fpOriginTouchable = right.fpOriginTouchable;

    // Track ID (and Parent ID) is not copied and set to zero for new track
    fTrackID  = 0;
    fParentID = 0;

    // CurrentStepNumber is set to be 0
    fCurrentStepNumber = 0;

    // Creator model ID
    fCreatorModelID = right.fCreatorModelID;

    // Parent resonance
    fParentResonanceDef = right.fParentResonanceDef;
    fParentResonanceID  = right.fParentResonanceID;
   
    // velocity information
    fVelocity = right.fVelocity;

    // dynamic particle information
    delete fpDynamicParticle;
    fpDynamicParticle = new G4DynamicParticle(*(right.fpDynamicParticle));

    // track status and flags for tracking
    fTrackStatus     = right.fTrackStatus;
    fBelowThreshold  = right.fBelowThreshold;
    fGoodForTracking = right.fGoodForTracking;

    // Step information (Step Length, Step Number, pointer to the Step,)
    // are not copied
    fpStep = nullptr;

    // vertex information
    fVtxPosition          = right.fVtxPosition;
    fpLVAtVertex          = right.fpLVAtVertex;
    fVtxKineticEnergy     = right.fVtxKineticEnergy;
    fVtxMomentumDirection = right.fVtxMomentumDirection;

    // CreatorProcess and UserInformation are not copied
    fpCreatorProcess = nullptr;
    delete fpUserInformation;
    fpUserInformation = nullptr;

    prev_mat      = right.prev_mat;
    groupvel      = right.groupvel;
    prev_velocity = right.prev_velocity;
    prev_momentum = right.prev_momentum;

    is_OpticalPhoton = right.is_OpticalPhoton;
    useGivenVelocity = right.useGivenVelocity;

    ClearAuxiliaryTrackInformation();
  }
  return *this;
}

// --------------------------------------------------------------------
void G4Track::CopyTrackInfo(const G4Track& right)
{
  *this = right;
}

// --------------------------------------------------------------------
G4double G4Track::CalculateVelocityForOpticalPhoton() const
{
  G4double velocity = c_light;

  G4Material* mat        = nullptr;
  G4bool update_groupvel = false;
  if(fpStep != nullptr)
  {
    mat = this->GetMaterial();  //   Fix for repeated volumes
  }
  else
  {
    if(fpTouchable)
    {
      mat = fpTouchable->GetVolume()->GetLogicalVolume()->GetMaterial();
    }
  }
  // check if previous step is in the same volume
  //  and get new GROUPVELOCITY table if necessary
  if((mat != nullptr) && ((mat != prev_mat) || (groupvel == nullptr)))
  {
    groupvel = nullptr;
    if(mat->GetMaterialPropertiesTable() != nullptr)
      groupvel = mat->GetMaterialPropertiesTable()->GetProperty(kGROUPVEL);
    update_groupvel = true;
  }
  prev_mat = mat;

  if(groupvel != nullptr)
  {
    // light velocity = c/(rindex+d(rindex)/d(log(E_phot)))
    // values stored in GROUPVEL material properties vector
    velocity = prev_velocity;

    // check if momentum is same as in the previous step
    //  and calculate group velocity if necessary
    G4double current_momentum = fpDynamicParticle->GetTotalMomentum();
    if(update_groupvel || (current_momentum != prev_momentum))
    {
      velocity      = groupvel->Value(current_momentum);
      prev_velocity = velocity;
      prev_momentum = current_momentum;
    }
  }

  return velocity;
}

// --------------------------------------------------------------------
void G4Track::SetAuxiliaryTrackInformation(G4int id,
              G4VAuxiliaryTrackInformation* info) const
{
  if(fpAuxiliaryTrackInformationMap == nullptr)
  {
    fpAuxiliaryTrackInformationMap =
      new std::map<G4int, G4VAuxiliaryTrackInformation*>;
  }
  if(G4PhysicsModelCatalog::GetModelIndex(id) < 0)
  {
    G4ExceptionDescription ED;
    ED << "Process/model ID <" << id << "> is invalid.";
    G4Exception("G4VAuxiliaryTrackInformation::G4VAuxiliaryTrackInformation()",
                "TRACK0982", FatalException, ED);
  }
  (*fpAuxiliaryTrackInformationMap)[id] = info;
}

// --------------------------------------------------------------------
G4VAuxiliaryTrackInformation*
G4Track::GetAuxiliaryTrackInformation(G4int id) const
{
  if(fpAuxiliaryTrackInformationMap == nullptr)
    return nullptr;
  
  auto itr = fpAuxiliaryTrackInformationMap->find(id);
  if(itr == fpAuxiliaryTrackInformationMap->cend())
    return nullptr;
  return (*itr).second;
}

// --------------------------------------------------------------------
void G4Track::RemoveAuxiliaryTrackInformation(G4int id)
{
  if(fpAuxiliaryTrackInformationMap != nullptr  &&
     fpAuxiliaryTrackInformationMap->find(id) != fpAuxiliaryTrackInformationMap->cend())
  {
    fpAuxiliaryTrackInformationMap->erase(id);
  }
}

// --------------------------------------------------------------------
void G4Track::RemoveAuxiliaryTrackInformation(G4String& name)
{
  if(fpAuxiliaryTrackInformationMap != nullptr)
  {
    G4int id = G4PhysicsModelCatalog::GetModelID(name);
    RemoveAuxiliaryTrackInformation(id);
  }
}

// --------------------------------------------------------------------
void G4Track::ClearAuxiliaryTrackInformation()
{
  if(fpAuxiliaryTrackInformationMap == nullptr)
    return;
  for(const auto& itr : *fpAuxiliaryTrackInformationMap)
  {
    delete itr.second;
  }
  delete fpAuxiliaryTrackInformationMap;
  fpAuxiliaryTrackInformationMap = nullptr;
}
