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
//
// $Id: G4BirkScintillation.cc,v 1.1 2008-02-01 14:56:06 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// BirkScintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4BirkScintillation.cc 
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     01.02.2008  
// Author:      Vladimir Grichine vbased on G4Scintillation of Peter Gumplinger
// Updated:     
//
// mail:        Vladimir.Grichine@cern.ch
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4BirkScintillation.hh"
#include "G4Scintillation.hh"

using namespace std;


////////////////////////////////////////////////////////////////////////
//
// Constructor

G4BirkScintillation::G4BirkScintillation(const G4String& processName,
                                       G4ProcessType type)
                  : G4Scintillation(processName, type)
{
}


////////////////////////////////////////////////////////////////////////
//
// Distructor

G4BirkScintillation::~G4BirkScintillation() 
{
	if (theFastIntegralTable != NULL) 
        {
	   theFastIntegralTable->clearAndDestroy();
           delete theFastIntegralTable;
	}
        if (theSlowIntegralTable != NULL) 
        {
           theSlowIntegralTable->clearAndDestroy();
           delete theSlowIntegralTable;
        }
}


//////////////////////////////////////////////////////////////////////////
//
// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()

G4VParticleChange*
G4BirkScintillation::AtRestDoIt(const G4Track& aTrack, const G4Step& aStep)


{
        return G4BirkScintillation::PostStepDoIt(aTrack, aStep);
}


////////////////////////////////////////////////////////////////////////////
//
// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is 
// generated according to the scintillation yield formula, distributed 
// evenly along the track segment and uniformly into 4pi.

G4VParticleChange*
G4BirkScintillation::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)


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
        ScintillationYield  *= YieldFactor;

        G4double ResolutionScale    = aMaterialPropertiesTable->
                                      GetConstProperty("RESOLUTIONSCALE");

	// Birks law saturation:

        G4double constBirks    = aMaterialPropertiesTable->
                                      GetConstProperty("BIRKSCONSTANT");

	G4double MeanNumberOfPhotons; 

	if( constBirks > DBL_MIN)       
	{
          G4double length      = aStep.GetStepLength();

          if(length > DBL_MIN) 
	  { 
	    MeanNumberOfPhotons  = ScintillationYield*TotalEnergyDeposit;
            MeanNumberOfPhotons /= 1 + constBirks*TotalEnergyDeposit/length;
	  }
          else
	  {
            MeanNumberOfPhotons = ScintillationYield*length/constBirks;
	  }  
	}
        else
	{
	  MeanNumberOfPhotons  = ScintillationYield*TotalEnergyDeposit;
	}
        G4int NumPhotons;

        if (MeanNumberOfPhotons > 10.) 
        {
          G4double sigma = ResolutionScale * sqrt(MeanNumberOfPhotons);
          NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);
        }
        else 
        {
          NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
        }
	if (NumPhotons <= 0) 
        {

	   // return unchanged particle and no secondaries 

	   aParticleChange.SetNumberOfSecondaries(0);

           return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

	////////////////////////////////////////////////////////////////

	aParticleChange.SetNumberOfSecondaries(NumPhotons);

	if (fTrackSecondariesFirst) 
        {
           if (aTrack.GetTrackStatus() == fAlive )
	  	   aParticleChange.ProposeTrackStatus(fSuspend);
        }
	
	////////////////////////////////////////////////////////////////

	G4int materialIndex = aMaterial->GetIndex();

	// Retrieve the Scintillation Integral for this material  
	// new G4PhysicsOrderedFreeVector allocated to hold CII's

        G4int Num = NumPhotons;

        for (G4int scnt = 1; scnt <= nscnt; scnt++) 
        {
            G4double ScintillationTime = 0.*ns;
            G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;

            if (scnt == 1) 
            {
               if (nscnt == 1) 
               {
                 if(Fast_Intensity)
                 {
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("FASTTIMECONSTANT");
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
                 }
                 if(Slow_Intensity)
                 {
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("SLOWTIMECONSTANT");
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
                 }
               }
               else 
               {
                 G4double YieldRatio = aMaterialPropertiesTable->
                                          GetConstProperty("YIELDRATIO");
                 if ( ExcitationRatio == 1.0 ) {
                    Num = G4int (min(YieldRatio,1.0) * NumPhotons);
                 }
                 else {
                    Num = G4int (min(ExcitationRatio,1.0) * NumPhotons);
                 }
                 ScintillationTime   = aMaterialPropertiesTable->
                                          GetConstProperty("FASTTIMECONSTANT");
                 ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
               }
            }
            else 
            {
               Num = NumPhotons - Num;
               ScintillationTime   =   aMaterialPropertiesTable->
                                          GetConstProperty("SLOWTIMECONSTANT");
               ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
            }

            if (!ScintillationIntegral) continue;
	
            // Max Scintillation Integral
 
	    G4double CIImax = ScintillationIntegral->GetMaxValue();
		
	    for (G4int i = 0; i < Num; i++) 
            {

		// Determine photon momentum

                G4double CIIvalue = G4UniformRand()*CIImax;
		G4double sampledMomentum = 
                              ScintillationIntegral->GetEnergy(CIIvalue);

		if (verboseLevel>1) 
                {
                   G4cout << "sampledMomentum = " << sampledMomentum << G4endl;
		   G4cout << "CIIvalue =        " << CIIvalue << G4endl;
		}

		// Generate random photon direction

                G4double cost = 1. - 2.*G4UniformRand();
                G4double sint = sqrt((1.-cost)*(1.+cost));

		G4double phi = twopi*G4UniformRand();
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

		phi = twopi*G4UniformRand();
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

                if (aParticle->GetDefinition()->GetPDGCharge() != 0) 
                {
                   rand = G4UniformRand();
                } 
                else 
                {
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

	if (verboseLevel>0) 
        {
	  G4cout << "\n Exiting from G4BirkScintillation::DoIt -- NumberOfSecondaries = " 
	     << aParticleChange.GetNumberOfSecondaries() << G4endl;
	}

	return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

