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
// $Id: G4OpWLS.cc,v 1.5 2004/12/10 18:53:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS.cc
// Description: Discrete Process -- Wavelength Shifting of Optical Photons
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     
// mail:        gum@triumf.ca
//              jparcham@phys.ualberta.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpWLS.hh"

/////////////////////////
// Class Implementation
/////////////////////////

/////////////////
// Constructors
/////////////////

G4OpWLS::G4OpWLS(const G4String& processName)
  : G4VDiscreteProcess(processName)
{
  theIntegralTable = 0;
 
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  
  BuildThePhysicsTable();
}

////////////////
// Destructors
////////////////

G4OpWLS::~G4OpWLS()
{
  if (theIntegralTable != 0) {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }
}

////////////
// Methods
////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
G4OpWLS::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if (verboseLevel>0) {
    G4cout << "\n** Photon absorbed! **" << G4endl;
  }
  
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  const G4MaterialPropertyVector* WLS_Intensity = 
    aMaterialPropertiesTable->GetProperty("WLSCOMPONENT"); 

  if (!WLS_Intensity)
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

  G4int NumPhotons = 1;

  aParticleChange.SetNumberOfSecondaries(NumPhotons);

  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the WLS Integral for this material
  // new G4PhysicsOrderedFreeVector allocated to hold CII's

  G4double WLSTime = 0.*ns;
  G4PhysicsOrderedFreeVector* WLSIntegral = 0;

  WLSTime   = aMaterialPropertiesTable->
    GetConstProperty("WLSTIMECONSTANT");
  WLSIntegral =
    (G4PhysicsOrderedFreeVector*)((*theIntegralTable)(materialIndex));
   
  // Max WLS Integral
  
  G4double CIImax = WLSIntegral->GetMaxValue();
  
  for (G4int i = 0; i < NumPhotons; i++) {
    
    // Determine photon momentum
    
    G4double CIIvalue = G4UniformRand()*CIImax;
    G4double sampledMomentum = 
      WLSIntegral->GetEnergy(CIIvalue);
    
    if (verboseLevel>1) {
      G4cout << "sampledMomentum = " << sampledMomentum << G4endl;
      G4cout << "CIIvalue =        " << CIIvalue << G4endl;
    }
    
    // Generate random photon direction
    
    G4double cost = 1. - 2.*G4UniformRand();
    G4double sint = std::sqrt((1.-cost)*(1.+cost));

    G4double phi = twopi*G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    
    G4double px = sint*cosp;
    G4double py = sint*sinp;
    G4double pz = cost;
    
    // Create photon momentum direction vector
    
    G4ParticleMomentum photonMomentum(px, py, pz);
    
    // Determine polarization of new photon
    
    G4double sx = cost*cosp;
    G4double sy = cost*sinp;
    G4double sz = -sint;
    
    G4ThreeVector photonPolarization(sx, sy, sz);
    
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);
    
    phi = twopi*G4UniformRand();
    sinp = std::sin(phi);
    cosp = std::cos(phi);
    
    photonPolarization = cosp * photonPolarization + sinp * perp;
    
    photonPolarization = photonPolarization.unit();
    
    // Generate a new photon:
    
    G4DynamicParticle* aWLSPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),
			    photonMomentum);
    aWLSPhoton->SetPolarization
      (photonPolarization.x(),
       photonPolarization.y(),
       photonPolarization.z());
    
    aWLSPhoton->SetKineticEnergy(sampledMomentum);
    
    // Generate new G4Track object:
    
    // Must give position of WLS optical photon
  
    G4double aSecondaryTime = (pPostStepPoint->GetGlobalTime()) + WLSTime;

    G4ThreeVector aSecondaryPosition = pPostStepPoint->GetPosition();

    G4Track* aSecondaryTrack = 
      new G4Track(aWLSPhoton,aSecondaryTime,aSecondaryPosition);
    
    aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);
    
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel>0) {
    G4cout << "\n Exiting from G4OpWLS::DoIt -- NumberOfSecondaries = " 
	   << aParticleChange.GetNumberOfSecondaries() << G4endl;  
  }
  
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the wavelength shifting process
// --------------------------------------------------
//

void G4OpWLS::BuildThePhysicsTable()
{
  if (theIntegralTable) return;
  
  const G4MaterialTable* theMaterialTable = 
    G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();
  
  // create new physics table
  
  if(!theIntegralTable)theIntegralTable = new G4PhysicsTable(numOfMaterials);
  
  // loop for materials
  
  for (G4int i=0 ; i < numOfMaterials; i++)
    {
      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
	new G4PhysicsOrderedFreeVector();
      
      // Retrieve vector of WLS wavelength intensity for
      // the material from the material's optical properties table.
      
      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
	aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {

	G4MaterialPropertyVector* theWLSVector = 
	  aMaterialPropertiesTable->GetProperty("WLSCOMPONENT");

	if (theWLSVector) {
	  
	  // Retrieve the first intensity point in vector
	  // of (photon momentum, intensity) pairs
	  
	  theWLSVector->ResetIterator();
	  ++(*theWLSVector);	// advance to 1st entry 
	  
	  G4double currentIN = theWLSVector->
	    GetProperty();
	  
	  if (currentIN >= 0.0) {

	    // Create first (photon momentum) 
	   
	    G4double currentPM = theWLSVector->
	      GetPhotonMomentum();
	    
	    G4double currentCII = 0.0;
	    
	    aPhysicsOrderedFreeVector->
	      InsertValues(currentPM , currentCII);
	    
	    // Set previous values to current ones prior to loop
	    
	    G4double prevPM  = currentPM;
	    G4double prevCII = currentCII;
	    G4double prevIN  = currentIN;
	    
	    // loop over all (photon momentum, intensity)
	    // pairs stored for this material
	    
	    while(++(*theWLSVector))
	      {
		currentPM = theWLSVector->
		  GetPhotonMomentum();
		
		currentIN=theWLSVector->
		  GetProperty();
		
		currentCII = 0.5 * (prevIN + currentIN);
		
		currentCII = prevCII +
		  (currentPM - prevPM) * currentCII;
		
		aPhysicsOrderedFreeVector->
		  InsertValues(currentPM, currentCII);
		
		prevPM  = currentPM;
		prevCII = currentCII;
		prevIN  = currentIN;
	      }
	  }
	}
      }
	// The WLS integral for a given material
	// will be inserted in the table according to the
	// position of the material in the material table.

	theIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
    }
}

// GetMeanFreePath
// ---------------
//
G4double G4OpWLS::GetMeanFreePath(const G4Track& aTrack,
 				         G4double ,
				         G4ForceCondition* )
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4double thePhotonMomentum = aParticle->GetTotalMomentum();

  G4MaterialPropertiesTable* aMaterialPropertyTable;
  G4MaterialPropertyVector* AttenuationLengthVector;
	
  G4double AttenuationLength = DBL_MAX;

  aMaterialPropertyTable = aMaterial->GetMaterialPropertiesTable();

  if ( aMaterialPropertyTable ) {
    AttenuationLengthVector = aMaterialPropertyTable->
      GetProperty("WLSABSLENGTH");
    if ( AttenuationLengthVector ){
      AttenuationLength = AttenuationLengthVector->
	GetProperty (thePhotonMomentum);
    }
    else {
      //             G4cout << "No WLS absorption length specified" << G4endl;
    }
  }
  else {
    //           G4cout << "No WLS absortion length specified" << G4endl;
  }
  
  return AttenuationLength;
}
