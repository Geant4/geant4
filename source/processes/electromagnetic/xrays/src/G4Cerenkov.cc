// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Cerenkov.cc,v 1.5 1999-11-05 01:56:49 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4Cerenkov.cc 
// Description: Continuous Process -- Generation of Cerenkov Photons
// Version:     2.1
// Created:     1996-02-21  
// Author:      Juliet Armstrong
// Updated:     1999-10-29 by Peter Gumplinger
//              > change: == into <= in GetContinuousStepLimit
//              1997-08-08 by Peter Gumplinger
//              > add protection against /0
//              > G4MaterialPropertiesTable; new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4Cerenkov.hh"

/////////////////////////
// Class Implementation  
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4Cerenkov::operator=(const G4Cerenkov &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4Cerenkov::G4Cerenkov(const G4String& processName)
           : G4VContinuousProcess(processName)
{
	fTrackSecondariesFirst = false;
	fMaxPhotons = 0;

        thePhysicsTable = NULL;

	if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << endl;
	}

	BuildThePhysicsTable();
}

// G4Cerenkov::G4Cerenkov(const G4Cerenkov &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4Cerenkov::~G4Cerenkov() 
{
	if (thePhysicsTable!= NULL) {
	   thePhysicsTable->clearAndDestroy();
           delete thePhysicsTable;
	}
}

        ////////////
        // Methods
        ////////////

// AlongStepDoIt
// -------------
//
G4VParticleChange*
G4Cerenkov::AlongStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The 
// parameters are then transformed into the Master Reference System, and 
// they are added to the particle change. 

{
	//////////////////////////////////////////////////////
	// Should we ensure that the material is dispersive?
	//////////////////////////////////////////////////////

        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

	G4StepPoint* pPreStepPoint = aStep.GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

	G4ThreeVector x0 = pPreStepPoint->GetPosition();
        G4ThreeVector p0 = pPreStepPoint->GetMomentumDirection();
	G4double t0 = pPreStepPoint->GetGlobalTime();

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable)
           return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);

	const G4MaterialPropertyVector* Rindex = 
                aMaterialPropertiesTable->GetProperty("RINDEX"); 
        if (!Rindex) 
	   return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);

	G4double MeanNumPhotons = 
                 GetAverageNumberOfPhotons(aParticle,aMaterial,Rindex);

        G4double step_length;
        step_length = aStep.GetStepLength();

	MeanNumPhotons = MeanNumPhotons * step_length;

	// RandPoisson is a utility class.  It provides functions
	// that act on HepRandom

	G4int NumPhotons = (G4int) RandPoisson::shoot(MeanNumPhotons);

	if (NumPhotons <= 0) {

		// return unchanged particle and no secondaries  

		aParticleChange.SetNumberOfSecondaries(0);
		
                return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
	}

	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

	if (fTrackSecondariesFirst)
		aParticleChange.SetStatusChange(fSuspend);
	
	////////////////////////////////////////////////////////////////

	G4double Pmin = Rindex->GetMinPhotonMomentum();
	G4double Pmax = Rindex->GetMaxPhotonMomentum();
	G4double dp = Pmax - Pmin;

	G4double nMax = Rindex->GetMaxProperty();

        G4double BetaInverse = aParticle->GetTotalEnergy() /
	   	 	       aParticle->GetTotalMomentum();

	G4double maxCos = BetaInverse / nMax; 
	G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

	for (G4int i = 0; i < NumPhotons; i++) {

		// Determine photon momentum

		G4double rand;
		G4double sampledMomentum, sampledRI; 
		G4double cosTheta, sin2Theta;
		
		// sample a momentum 

		do {
			rand = G4UniformRand();	
			sampledMomentum = Pmin + rand * dp; 
			sampledRI = Rindex->GetProperty(sampledMomentum);
			cosTheta = BetaInverse / sampledRI;  

			sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
			rand = G4UniformRand();	

		} while (rand*maxSin2 > sin2Theta);

		// Generate random position of photon on cone surface 
		// defined by Theta 

		rand = G4UniformRand();

		G4double phi = 2*M_PI*rand;
		G4double sinPhi = sin(phi);
		G4double cosPhi = cos(phi);

		// calculate x,y, and z components of photon momentum 
		// (in coord system with primary particle direction 
		//  aligned with the z axis)

		G4double sinTheta = sqrt(sin2Theta); 
		G4double px = sinTheta*cosPhi;
		G4double py = sinTheta*sinPhi;
		G4double pz = cosTheta;

		// Create photon momentum direction vector 
		// The momentum direction is still with respect
	 	// to the coordinate system where the primary
		// particle direction is aligned with the z axis  

		G4ParticleMomentum photonMomentum(px, py, pz);

		// Rotate momentum direction back to global reference
		// system 

                photonMomentum.rotateUz(p0);

		// Determine polarization of new photon 

		G4double sx = cosTheta*cosPhi;
		G4double sy = cosTheta*sinPhi; 
		G4double sz = -sinTheta;

		G4ThreeVector photonPolarization(sx, sy, sz);

		// Rotate back to original coord system 

                photonPolarization.rotateUz(p0);
		
                // Generate a new photon:

                G4DynamicParticle* aCerenkovPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
  					                 photonMomentum);
		aCerenkovPhoton->SetPolarization
				     (photonPolarization.x(),
				      photonPolarization.y(),
				      photonPolarization.z());

		aCerenkovPhoton->SetKineticEnergy(sampledMomentum);

                // Generate new G4Track object:

		rand = G4UniformRand();

		G4double delta = rand * aStep.GetStepLength();
		G4ThreeVector aSecondaryPosition = x0 + delta * p0;

		G4double deltaTime = delta /
                       ((pPreStepPoint->GetVelocity()+
                         pPostStepPoint->GetVelocity())/2.);

                G4double aSecondaryTime = t0 + deltaTime;

		G4Track* aSecondaryTrack = 
		new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);

                aSecondaryTrack->SetTouchable(pPreStepPoint->
                                           GetTouchable());

                aSecondaryTrack->SetParentID(aTrack.GetTrackID());

		aParticleChange.AddSecondary(aSecondaryTrack);
	}

	if (verboseLevel>0) {
	G4cout << "\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = " 
	     << aParticleChange.GetNumberOfSecondaries() << endl;
	}

	return G4VContinuousProcess::AlongStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the Cerenkov process
// ---------------------------------------------
//

void G4Cerenkov::BuildThePhysicsTable()
{
	if (thePhysicsTable) return;

	const G4MaterialTable* theMaterialTable=
	 		       G4Material::GetMaterialTable();
	G4int numOfMaterials = theMaterialTable->length();

	// create new physics table
	
	thePhysicsTable = new G4PhysicsTable(numOfMaterials);

	// loop for materials

	for (G4int i=0 ; i < numOfMaterials; i++)
	{
		G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
					new G4PhysicsOrderedFreeVector();

		// Retrieve vector of refraction indices for the material
		// from the material's optical properties table 

		G4Material* aMaterial = (*theMaterialTable)(i);

		G4MaterialPropertiesTable* aMaterialPropertiesTable =
				aMaterial->GetMaterialPropertiesTable();

		if (aMaterialPropertiesTable) {

		   G4MaterialPropertyVector* theRefractionIndexVector = 
		    	   aMaterialPropertiesTable->GetProperty("RINDEX");

		   if (theRefractionIndexVector) {
		
		      // Retrieve the first refraction index in vector
		      // of (photon momentum, refraction index) pairs 

		      theRefractionIndexVector->ResetIterator();
		      ++(*theRefractionIndexVector);	// advance to 1st entry 

		      G4double currentRI = theRefractionIndexVector->
		  			   GetProperty();

		      if (currentRI > 1.0) {

			 // Create first (photon momentum, Cerenkov Integral)
			 // pair  

			 G4double currentPM = theRefractionIndexVector->
			 			 GetPhotonMomentum();
			 G4double currentCAI = 0.0;

			 aPhysicsOrderedFreeVector->
			 	 InsertValues(currentPM , currentCAI);

			 // Set previous values to current ones prior to loop

			 G4double prevPM  = currentPM;
			 G4double prevCAI = currentCAI;
                	 G4double prevRI  = currentRI;

			 // loop over all (photon momentum, refraction index)
			 // pairs stored for this material  

			 while(++(*theRefractionIndexVector))
			 {
				currentRI=theRefractionIndexVector->	
						GetProperty();

				currentPM = theRefractionIndexVector->
						GetPhotonMomentum();

				currentCAI = 0.5*(1.0/(prevRI*prevRI) +
					          1.0/(currentRI*currentRI));

				currentCAI = prevCAI + 
					     (currentPM - prevPM) * currentCAI;

				aPhysicsOrderedFreeVector->
				    InsertValues(currentPM, currentCAI);

				prevPM  = currentPM;
				prevCAI = currentCAI;
				prevRI  = currentRI;
			 }

		      }
		   }
		}

	// The Cerenkov integral for a given material
	// will be inserted in thePhysicsTable
	// according to the position of the material in
	// the material table. 

	thePhysicsTable->insertAt(i,aPhysicsOrderedFreeVector); 

	}
}

// GetContinuousStepLimit
// ----------------------
//

G4double 
G4Cerenkov::GetContinuousStepLimit(const G4Track& aTrack,
				   G4double  ,
				   G4double  ,
                                   G4double& )
{
	// If user has defined an average maximum number of photons to
	// be generated in a Step, then return the Step length for that 
	// number of photons. 
 
	if (fMaxPhotons <= 0) return DBL_MAX;

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

        G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
        if (!aMaterialPropertiesTable) return DBL_MAX;

        const G4MaterialPropertyVector* Rindex =
                aMaterialPropertiesTable->GetProperty("RINDEX");
        if (!Rindex) return DBL_MAX;

	G4double MeanNumPhotons = 
                 GetAverageNumberOfPhotons(aParticle,aMaterial,Rindex);

        if(MeanNumPhotons <= 0.0) return DBL_MAX;

	G4double StepLimit = fMaxPhotons / MeanNumPhotons;

	return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium. 
//             ^^^^^^^^^^

G4double 
G4Cerenkov::GetAverageNumberOfPhotons(const G4DynamicParticle* aParticle, 
			      const G4Material* aMaterial,
			      const G4MaterialPropertyVector* Rindex) const
{
	const G4double Rfact = 369.81/(eV * cm);

        if(aParticle->GetTotalMomentum() <= 0.0)return 0.0;

	G4double BetaInverse = aParticle->GetTotalEnergy() /
			       aParticle->GetTotalMomentum();

	// Vectors used in computation of Cerenkov Angle Integral:
	// 	- Refraction Indices for the current material
	//	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
 
	G4int materialIndex = G4Material::GetMaterialTable()->index(aMaterial);

	// Retrieve the Cerenkov Angle Integrals for this material  

	G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
	(G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));

	// Min and Max photon momenta  
	G4double Pmin = Rindex->GetMinPhotonMomentum();
	G4double Pmax = Rindex->GetMaxPhotonMomentum();

	// Min and Max Refraction Indices 
	G4double nMin = Rindex->GetMinProperty();	
	G4double nMax = Rindex->GetMaxProperty();

	// Max Cerenkov Angle Integral 
	G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

	G4double dp, ge;

	// If n(Pmax) < 1/Beta -- no photons generated 

	if (nMax < BetaInverse) {
		dp = 0;
		ge = 0;
	} 

	// otherwise if n(Pmin) >= 1/Beta -- photons generated  

	else if (nMin > BetaInverse) {
		dp = Pmax - Pmin;	
		ge = CAImax; 
	} 

	// If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
	// we need to find a P such that the value of n(P) == 1/Beta.
	// Interpolation is performed by the GetPhotonMomentum() and
	// GetProperty() methods of the G4MaterialPropertiesTable and
	// the GetValue() method of G4PhysicsVector.  

	else {
		Pmin = Rindex->GetPhotonMomentum(BetaInverse);
		dp = Pmax - Pmin;

		// need boolean for current implementation of G4PhysicsVector
		// ==> being phased out
		G4bool isOutRange;
		G4double CAImin = CerenkovAngleIntegrals->
                                  GetValue(Pmin, isOutRange);
		ge = CAImax - CAImin;

		if (verboseLevel>0) {
			G4cout << "CAImin = " << CAImin << endl;
			G4cout << "ge = " << ge << endl;
		}
	}
	
	// particle charge 
	G4double charge = aParticle->GetDefinition()->GetPDGCharge();

	// Calculate number of photons 
	G4double NumPhotons =
                 Rfact * charge*charge * (dp - ge * BetaInverse*BetaInverse);

	return NumPhotons;		
}
