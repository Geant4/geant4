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
// $Id: G4OpRayleigh.cc,v 1.19 2010-10-29 23:18:35 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpRayleigh.cc 
// Description: Discrete Process -- Rayleigh scattering of optical 
//		photons  
// Version:     1.0
// Created:     1996-05-31  
// Author:      Juliet Armstrong
// Updated:     2010-06-11 - Fix Bug 207; Thanks to Xin Qian
//              (Kellogg Radiation Lab of Caltech)
//              2005-07-28 - add G4ProcessType to constructor
//              2001-10-18 by Peter Gumplinger
//              eliminate unused variable warning on Linux (gcc-2.95.2)
//              2001-09-18 by mma
//		>numOfMaterials=G4Material::GetNumberOfMaterials() in BuildPhy
//              2001-01-30 by Peter Gumplinger
//              > allow for positiv and negative CosTheta and force the
//              > new momentum direction to be in the same plane as the
//              > new and old polarization vectors
//              2001-01-29 by Peter Gumplinger
//              > fix calculation of SinTheta (from CosTheta)
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpProcessSubType.hh"

#include "G4OpRayleigh.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4OpRayleigh::operator=(const G4OpRayleigh &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4OpRayleigh::G4OpRayleigh(const G4String& processName, G4ProcessType type)
           : G4VDiscreteProcess(processName, type)
{
        SetProcessSubType(fOpRayleigh);

        thePhysicsTable = 0;

        DefaultWater = false;

        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        BuildThePhysicsTable();
}

// G4OpRayleigh::G4OpRayleigh(const G4OpRayleigh &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4OpRayleigh::~G4OpRayleigh()
{
        if (thePhysicsTable!= 0) {
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
G4OpRayleigh::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

        if (verboseLevel>0) {
                G4cout << "Scattering Photon!" << G4endl;
                G4cout << "Old Momentum Direction: "
                       << aParticle->GetMomentumDirection() << G4endl;
                G4cout << "Old Polarization: "
                       << aParticle->GetPolarization() << G4endl;
        }

        G4double cosTheta;
        G4ThreeVector OldMomentumDirection, NewMomentumDirection;
        G4ThreeVector OldPolarization, NewPolarization;

        do {
           // Try to simulate the scattered photon momentum direction
           // w.r.t. the initial photon momentum direction

           G4double CosTheta = G4UniformRand();
           G4double SinTheta = std::sqrt(1.-CosTheta*CosTheta);
           // consider for the angle 90-180 degrees
           if (G4UniformRand() < 0.5) CosTheta = -CosTheta;

           // simulate the phi angle
           G4double rand = twopi*G4UniformRand();
           G4double SinPhi = std::sin(rand);
           G4double CosPhi = std::cos(rand);

           // start constructing the new momentum direction
	   G4double unit_x = SinTheta * CosPhi; 
	   G4double unit_y = SinTheta * SinPhi;  
	   G4double unit_z = CosTheta; 
	   NewMomentumDirection.set (unit_x,unit_y,unit_z);

           // Rotate the new momentum direction into global reference system
           OldMomentumDirection = aParticle->GetMomentumDirection();
           OldMomentumDirection = OldMomentumDirection.unit();
           NewMomentumDirection.rotateUz(OldMomentumDirection);
           NewMomentumDirection = NewMomentumDirection.unit();

           // calculate the new polarization direction
           // The new polarization needs to be in the same plane as the new
           // momentum direction and the old polarization direction
           OldPolarization = aParticle->GetPolarization();
           G4double constant = -1./NewMomentumDirection.dot(OldPolarization);

           NewPolarization = NewMomentumDirection + constant*OldPolarization;
           NewPolarization = NewPolarization.unit();

           // There is a corner case, where the Newmomentum direction
           // is the same as oldpolariztion direction:
           // random generate the azimuthal angle w.r.t. Newmomentum direction
           if (NewPolarization.mag() == 0.) {
              rand = G4UniformRand()*twopi;
              NewPolarization.set(std::cos(rand),std::sin(rand),0.);
              NewPolarization.rotateUz(NewMomentumDirection);
           } else {
              // There are two directions which are perpendicular
              // to the new momentum direction
              if (G4UniformRand() < 0.5) NewPolarization = -NewPolarization;
           }
	  
	   // simulate according to the distribution cos^2(theta)
           cosTheta = NewPolarization.dot(OldPolarization);
        } while (std::pow(cosTheta,2) < G4UniformRand());

        aParticleChange.ProposePolarization(NewPolarization);
        aParticleChange.ProposeMomentumDirection(NewMomentumDirection);

        if (verboseLevel>0) {
                G4cout << "New Polarization: " 
                     << NewPolarization << G4endl;
                G4cout << "Polarization Change: "
                     << *(aParticleChange.GetPolarization()) << G4endl;  
                G4cout << "New Momentum Direction: " 
                     << NewMomentumDirection << G4endl;
                G4cout << "Momentum Change: "
                     << *(aParticleChange.GetMomentumDirection()) << G4endl; 
        }

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// BuildThePhysicsTable for the Rayleigh Scattering process
// --------------------------------------------------------
//
void G4OpRayleigh::BuildThePhysicsTable()
{
//      Builds a table of scattering lengths for each material

        if (thePhysicsTable) return;

        const G4MaterialTable* theMaterialTable=
                               G4Material::GetMaterialTable();
        G4int numOfMaterials = G4Material::GetNumberOfMaterials();

        // create a new physics table

        thePhysicsTable = new G4PhysicsTable(numOfMaterials);

        // loop for materials

        for (G4int i=0 ; i < numOfMaterials; i++)
        {
            G4PhysicsOrderedFreeVector* ScatteringLengths = NULL;

            G4MaterialPropertiesTable *aMaterialPropertiesTable =
                         (*theMaterialTable)[i]->GetMaterialPropertiesTable();
                                                                                
            if(aMaterialPropertiesTable){

              G4MaterialPropertyVector* AttenuationLengthVector =
                            aMaterialPropertiesTable->GetProperty("RAYLEIGH");

              if(!AttenuationLengthVector){

                if ((*theMaterialTable)[i]->GetName() == "Water")
                {
		   // Call utility routine to Generate
		   // Rayleigh Scattering Lengths

                   DefaultWater = true;

                   ScatteringLengths =
		   RayleighAttenuationLengthGenerator(aMaterialPropertiesTable);
                }
              }
	    }

	    thePhysicsTable->insertAt(i,ScatteringLengths);
        } 
}

// GetMeanFreePath()
// -----------------
//
G4double G4OpRayleigh::GetMeanFreePath(const G4Track& aTrack,
                                     G4double ,
                                     G4ForceCondition* )
{
        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

        G4double thePhotonEnergy = aParticle->GetTotalEnergy();

        G4double AttenuationLength = DBL_MAX;

        if (aMaterial->GetName() == "Water" && DefaultWater){

           G4bool isOutRange;

           AttenuationLength =
                (*thePhysicsTable)(aMaterial->GetIndex())->
                           GetValue(thePhotonEnergy, isOutRange);
        }
        else {

           G4MaterialPropertiesTable* aMaterialPropertyTable =
                           aMaterial->GetMaterialPropertiesTable();

           if(aMaterialPropertyTable){
             G4MaterialPropertyVector* AttenuationLengthVector =
                   aMaterialPropertyTable->GetProperty("RAYLEIGH");
             if(AttenuationLengthVector){
               AttenuationLength = AttenuationLengthVector ->
                                    GetProperty(thePhotonEnergy);
             }
             else{
//               G4cout << "No Rayleigh scattering length specified" << G4endl;
             }
           }
           else{
//             G4cout << "No Rayleigh scattering length specified" << G4endl; 
           }
        }

        return AttenuationLength;
}

// RayleighAttenuationLengthGenerator()
// ------------------------------------
// Private method to compute Rayleigh Scattering Lengths (for water)
//
G4PhysicsOrderedFreeVector* 
G4OpRayleigh::RayleighAttenuationLengthGenerator(G4MaterialPropertiesTable *aMPT) 
{
        // Physical Constants

        // isothermal compressibility of water
        G4double betat = 7.658e-23*m3/MeV;

        // K Boltzman
        G4double kboltz = 8.61739e-11*MeV/kelvin;

        // Temperature of water is 10 degrees celsius
        // conversion to kelvin:
        // TCelsius = TKelvin - 273.15 => 273.15 + 10 = 283.15
        G4double temp = 283.15*kelvin;

        // Retrieve vectors for refraction index
        // and photon energy from the material properties table

        G4MaterialPropertyVector* Rindex = aMPT->GetProperty("RINDEX");

        G4double refsq;
        G4double e;
        G4double xlambda;
        G4double c1, c2, c3, c4;
        G4double Dist;
        G4double refraction_index;

        G4PhysicsOrderedFreeVector *RayleighScatteringLengths = 
				new G4PhysicsOrderedFreeVector();

        if (Rindex ) {

           Rindex->ResetIterator();

           while (++(*Rindex)) {

                e = (Rindex->GetPhotonEnergy());

                refraction_index = Rindex->GetProperty();
                refsq = refraction_index*refraction_index;
                xlambda = h_Planck*c_light/e;

	        if (verboseLevel>0) {
        	        G4cout << Rindex->GetPhotonEnergy() << " MeV\t";
                	G4cout << xlambda << " mm\t";
		}

                c1 = 1 / (6.0 * pi);
                c2 = std::pow((2.0 * pi / xlambda), 4);
                c3 = std::pow( ( (refsq - 1.0) * (refsq + 2.0) / 3.0 ), 2);
                c4 = betat * temp * kboltz;

                Dist = 1.0 / (c1*c2*c3*c4);

	        if (verboseLevel>0) {
	                G4cout << Dist << " mm" << G4endl;
		}
                RayleighScatteringLengths->
			InsertValues(Rindex->GetPhotonEnergy(), Dist);
           }

        }

	return RayleighScatteringLengths;
}
