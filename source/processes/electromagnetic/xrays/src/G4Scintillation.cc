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
// $Id: G4Scintillation.cc,v 1.16 2002-11-26 00:52:13 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4Scintillation.cc 
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07  
// Author:      Peter Gumplinger
// Updated:     2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
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
                  : G4VRestDiscreteProcess(processName)
{
	fTrackSecondariesFirst = false;

        YieldFactor = 1.0;
        ExcitationRatio = 1.0;

        theFastIntegralTable = NULL;
        theSlowIntegralTable = NULL;

	if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
	}

	BuildThePhysicsTable();
}

        ////////////////
        // Destructors
        ////////////////

G4Scintillation::~G4Scintillation() 
{
	if (theFastIntegralTable != NULL) {
	   theFastIntegralTable->clearAndDestroy();
           delete theFastIntegralTable;
	}
        if (theSlowIntegralTable != NULL) {
           theSlowIntegralTable->clearAndDestroy();
           delete theSlowIntegralTable;
        }
}

        ////////////
        // Methods
        ////////////

// AtRestDoIt
// ----------
//
G4VParticleChange*
G4Scintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()

{
        return G4Scintillation::PostStepDoIt(aTrack, aStep);
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
G4Scintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is 
// generated according to the scintillation yield formula, distributed 
// evenly along the track segment and uniformly into 4pi.

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
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

	const G4MaterialPropertyVector* Fast_Intensity = 
                aMaterialPropertiesTable->GetProperty("FASTCOMPONENT"); 
        const G4MaterialPropertyVector* Slow_Intensity =
                aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

        if (!Fast_Intensity && !Slow_Intensity )
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

        G4int nscnt = 1;
        if (Fast_Intensity && Slow_Intensity) nscnt = 2;

        G4double ScintillationYield = aMaterialPropertiesTable->
                                      GetConstProperty("SCINTILLATIONYIELD");
        G4double ResolutionScale    = aMaterialPropertiesTable->
                                      GetConstProperty("RESOLUTIONSCALE");

        ScintillationYield = YieldFactor * ScintillationYield;

	G4double MeanNumPhotons = ScintillationYield * TotalEnergyDeposit;

        G4int NumPhotons;
        if (MeanNumPhotons > 10.) {
           G4double sigma = ResolutionScale * sqrt(MeanNumPhotons);
           NumPhotons = G4int(G4RandGauss::shoot(MeanNumPhotons,sigma)+0.5);
        }
        else {
           NumPhotons = G4int(G4Poisson(MeanNumPhotons));
        }

	if (NumPhotons <= 0) {

	   // return unchanged particle and no secondaries 

	   aParticleChange.SetNumberOfSecondaries(0);

           return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

	if (fTrackSecondariesFirst) {
           if (aTrack.GetTrackStatus() == fAlive )
	  	   aParticleChange.SetStatusChange(fSuspend);
        }
	
	////////////////////////////////////////////////////////////////

	G4int materialIndex = aMaterial->GetIndex();

	// Retrieve the Scintillation Integral for this material  
	// new G4PhysicsOrderedFreeVector allocated to hold CII's

        G4int Num = NumPhotons;

        for (G4int scnt = 1; scnt <= nscnt; scnt++) {

            G4double ScintillationTime = 0.*ns;
            G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;

            if (scnt == 1) {
               if (nscnt == 1) {
                 if(Fast_Intensity){
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("FASTTIMECONSTANT");
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
                 }
                 if(Slow_Intensity){
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("SLOWTIMECONSTANT");
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
                 }
               }
               else {
                 G4double YieldRatio = aMaterialPropertiesTable->
                                          GetConstProperty("YIELDRATIO");
                 if ( ExcitationRatio == 1.0 ) {
                    Num = G4int (G4std::min(YieldRatio,1.0) * NumPhotons);
                 }
                 else {
                    Num = G4int (G4std::min(ExcitationRatio,1.0) * NumPhotons);
                 }
                 ScintillationTime   = aMaterialPropertiesTable->
                                          GetConstProperty("FASTTIMECONSTANT");
                 ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
               }
            }
            else {
               Num = NumPhotons - Num;
               ScintillationTime   =   aMaterialPropertiesTable->
                                          GetConstProperty("SLOWTIMECONSTANT");
               ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
            }

            if (!ScintillationIntegral) continue;
	
            // Max Scintillation Integral
 
	    G4double CIImax = ScintillationIntegral->GetMaxValue();
		
	    for (G4int i = 0; i < Num; i++) {

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
        }

	if (verboseLevel>0) {
	G4cout << "\n Exiting from G4Scintillation::DoIt -- NumberOfSecondaries = " 
	     << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}

	return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the scintillation process
// --------------------------------------------------
//

void G4Scintillation::BuildThePhysicsTable()
{
	if (theFastIntegralTable && theSlowIntegralTable) return;

	const G4MaterialTable* theMaterialTable = 
                               G4Material::GetMaterialTable();
	G4int numOfMaterials = G4Material::GetNumberOfMaterials();

	// create new physics table
	
	if(!theFastIntegralTable)theFastIntegralTable = new G4PhysicsTable(numOfMaterials);
        if(!theSlowIntegralTable)theSlowIntegralTable = new G4PhysicsTable(numOfMaterials);

	// loop for materials

	for (G4int i=0 ; i < numOfMaterials; i++)
	{
		G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
					new G4PhysicsOrderedFreeVector();
                G4PhysicsOrderedFreeVector* bPhysicsOrderedFreeVector =
                                        new G4PhysicsOrderedFreeVector();

		// Retrieve vector of scintillation wavelength intensity for
                // the material from the material's optical properties table.

		G4Material* aMaterial = (*theMaterialTable)[i];

		G4MaterialPropertiesTable* aMaterialPropertiesTable =
				aMaterial->GetMaterialPropertiesTable();

		if (aMaterialPropertiesTable) {

		   G4MaterialPropertyVector* theFastLightVector = 
		   aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");

		   if (theFastLightVector) {
		
		      // Retrieve the first intensity point in vector
		      // of (photon momentum, intensity) pairs 

		      theFastLightVector->ResetIterator();
		      ++(*theFastLightVector);	// advance to 1st entry 

		      G4double currentIN = theFastLightVector->
		  			   GetProperty();

		      if (currentIN >= 0.0) {

			 // Create first (photon momentum, Scintillation 
                         // Integral pair  

			 G4double currentPM = theFastLightVector->
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

			 while(++(*theFastLightVector))
			 {
				currentPM = theFastLightVector->
						GetPhotonMomentum();

				currentIN=theFastLightVector->	
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

                   G4MaterialPropertyVector* theSlowLightVector =
                   aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

                   if (theSlowLightVector) {

                      // Retrieve the first intensity point in vector
                      // of (photon momentum, intensity) pairs

                      theSlowLightVector->ResetIterator();
                      ++(*theSlowLightVector);  // advance to 1st entry

                      G4double currentIN = theSlowLightVector->
                                           GetProperty();

                      if (currentIN >= 0.0) {

                         // Create first (photon momentum, Scintillation
                         // Integral pair

                         G4double currentPM = theSlowLightVector->
                                                 GetPhotonMomentum();

                         G4double currentCII = 0.0;

                         bPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon momentum, intensity)
                         // pairs stored for this material

                         while(++(*theSlowLightVector))
                         {
                                currentPM = theSlowLightVector->
                                                GetPhotonMomentum();

                                currentIN=theSlowLightVector->
                                                GetProperty();

                                currentCII = 0.5 * (prevIN + currentIN);

                                currentCII = prevCII +
                                             (currentPM - prevPM) * currentCII;

                                bPhysicsOrderedFreeVector->
                                    InsertValues(currentPM, currentCII);

                                prevPM  = currentPM;
                                prevCII = currentCII;
                                prevIN  = currentIN;
                         }

                      }
                   }
		}

	// The scintillation integral(s) for a given material
	// will be inserted in the table(s) according to the
	// position of the material in the material table.

	theFastIntegralTable->insertAt(i,aPhysicsOrderedFreeVector); 
        theSlowIntegralTable->insertAt(i,bPhysicsOrderedFreeVector);

	}
}

// GetMeanFreePath
// ---------------
//

G4double G4Scintillation::GetMeanFreePath(const G4Track& aTrack,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;

	return DBL_MAX;

}

// GetMeanLifeTime
// ---------------
//

G4double G4Scintillation::GetMeanLifeTime(const G4Track& aTrack,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;

        return DBL_MAX;

}
