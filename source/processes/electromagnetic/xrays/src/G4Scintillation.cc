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
// $Id: G4Scintillation.cc,v 1.9 2002-05-09 20:35:04 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4Scintillation.cc 
// Description: Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07  
// Author:      Peter Gumplinger
// Updated:     2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma) 
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4Scintillation.hh"

/////////////////////////
// Class Implementation  
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4Scintillation::operator=(const G4Scintillation &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4Scintillation::G4Scintillation(const G4String& processName)
                  : G4VDiscreteProcess(processName)
{
	fTrackSecondariesFirst = false;

        ScintillationYield = 0.0;
        ScintillationTime  = 0.0;
        ResolutionScale    = 1.0;

        thePhysicsTable = NULL;

	if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
	}

	BuildThePhysicsTable();
}

// G4Scintillation::G4Scintillation(const G4Scintillation &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4Scintillation::~G4Scintillation() 
{
	if (thePhysicsTable != NULL) {
	   thePhysicsTable->clearAndDestroy();
           delete thePhysicsTable;
	}
}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
G4Scintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking step of a charged particle
// in a scintillator. A Gaussian-distributed number of photons is generated
// according to the scintillation yield formula, distributed evenly along 
// the track segment and uniformly into 4pi.

{
        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

	G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

	G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
	G4double      t0 = pPreStepPoint->GetGlobalTime();

        G4double TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable)
             return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	const G4MaterialPropertyVector* Intensity = 
                aMaterialPropertiesTable->GetProperty("SCINTILLATION"); 
        if (!Intensity) 
  	     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

	G4double MeanNumPhotons = ScintillationYield * TotalEnergyDeposit;

	G4int NumPhotons = (G4int) MeanNumPhotons +
             int( ResolutionScale * G4RandGauss::shoot(0.0,sqrt(MeanNumPhotons)));

	if (NumPhotons <= 0) {

		// return unchanged particle and no secondaries  

		aParticleChange.SetNumberOfSecondaries(0);
		
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

	if (fTrackSecondariesFirst)
		aParticleChange.SetStatusChange(fSuspend);
	
	////////////////////////////////////////////////////////////////

	G4int materialIndex = aMaterial->GetIndex();

	// Retrieve the Scintillation Integral for this material  
	// new G4PhysicsOrderedFreeVector allocated to hold CII's

	G4PhysicsOrderedFreeVector* ScintillationIntegral =
	(G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));
	
        // Max Scintillation Integral
 
	G4double CIImax = ScintillationIntegral->GetMaxValue();
		
	for (G4int i = 0; i < NumPhotons; i++) {

		// Determine photon momentum

                G4double CIIvalue = G4UniformRand()*CIImax;
		G4double sampledMomentum = 
                              ScintillationIntegral->GetEnergy(CIIvalue);

		if (verboseLevel>1) {
                   G4cout << "sampledMomentum = " << sampledMomentum << G4endl;
		   G4cout << "CIIvalue =        " << CIIvalue << G4endl;
		}

		// Generate random photon direction

                G4double cost = 1. - 2.*G4UniformRand();
                G4double sint = sqrt((1.-cost)*(1.+cost));

		G4double phi = 2*M_PI*G4UniformRand();
		G4double sinp = sin(phi);
		G4double cosp = cos(phi);

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

		phi = 2*M_PI*G4UniformRand();
		sinp = sin(phi);
		cosp = cos(phi);

                photonPolarization = cosp * photonPolarization + sinp * perp;

                photonPolarization = photonPolarization.unit();

                // Generate a new photon:

                G4DynamicParticle* aScintillationPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
  					                 photonMomentum);
		aScintillationPhoton->SetPolarization
				     (photonPolarization.x(),
				      photonPolarization.y(),
				      photonPolarization.z());

		aScintillationPhoton->SetKineticEnergy(sampledMomentum);

                // Generate new G4Track object:

                G4double rand;

                if (aParticle->GetDefinition()->GetPDGCharge() != 0) {
                   rand = G4UniformRand();
                } else {
                   rand = 1.0;
                }

                G4double delta = rand * aStep.GetStepLength();
		G4double deltaTime = delta /
                       ((pPreStepPoint->GetVelocity()+
                         pPostStepPoint->GetVelocity())/2.);

                deltaTime = deltaTime - 
                            ScintillationTime * log( G4UniformRand() );

                G4double aSecondaryTime = t0 + deltaTime;

                G4ThreeVector aSecondaryPosition =
                                    x0 + rand * aStep.GetDeltaPosition();

		G4Track* aSecondaryTrack = 
		new G4Track(aScintillationPhoton,aSecondaryTime,aSecondaryPosition);

                aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

                aSecondaryTrack->SetParentID(aTrack.GetTrackID());

		aParticleChange.AddSecondary(aSecondaryTrack);

	}

	if (verboseLevel>0) {
	G4cout << "\n Exiting from G4Scintillation::DoIt -- NumberOfSecondaries = " 
	     << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}

	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the scintillation process
// --------------------------------------------------
//

void G4Scintillation::BuildThePhysicsTable()
{
	if (thePhysicsTable) return;

	const G4MaterialTable* theMaterialTable = 
                               G4Material::GetMaterialTable();
	G4int numOfMaterials = G4Material::GetNumberOfMaterials();

	// create new physics table
	
	thePhysicsTable = new G4PhysicsTable(numOfMaterials);

	// loop for materials

	for (G4int i=0 ; i < numOfMaterials; i++)
	{
		G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
					new G4PhysicsOrderedFreeVector();

		// Retrieve vector of scintillation wavelength intensity
                // for the material from the material's optical
                // properties table 

		G4Material* aMaterial = (*theMaterialTable)[i];

		G4MaterialPropertiesTable* aMaterialPropertiesTable =
				aMaterial->GetMaterialPropertiesTable();

		if (aMaterialPropertiesTable) {

		   G4MaterialPropertyVector* theScintillationLightVector = 
		   aMaterialPropertiesTable->GetProperty("SCINTILLATION");

		   if (theScintillationLightVector) {
		
		      // Retrieve the first intensity point in vector
		      // of (photon momentum, intensity) pairs 

		      theScintillationLightVector->ResetIterator();
		      ++(*theScintillationLightVector);	// advance to 1st entry 

		      G4double currentIN = theScintillationLightVector->
		  			   GetProperty();

		      if (currentIN >= 0.0) {

			 // Create first (photon momentum, Scintillation 
                         // Integral pair  

			 G4double currentPM = theScintillationLightVector->
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

			 while(++(*theScintillationLightVector))
			 {
				currentPM = theScintillationLightVector->
						GetPhotonMomentum();

				currentIN=theScintillationLightVector->	
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

	// The scintillation integral for a given material
	// will be inserted in thePhysicsTable
	// according to the position of the material in
	// the material table. 

	thePhysicsTable->insertAt(i,aPhysicsOrderedFreeVector); 

	}
}

// GetMeanFreePath
// ---------------
//

G4double G4Scintillation::GetMeanFreePath(const G4Track& aTrack,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = Forced;

	return DBL_MAX;

}
