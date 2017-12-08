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
////////////////////////////////////////////////////////////////////////
//
// File G4OpMieHG.hh
// Description: Discrete Process -- Mie Scattering of Optical Photons
// Created: 2010-07-03
// Author: Xin Qian
// Based on work from Vlasios Vasileiou
//
// This subroutine will mimic the Mie scattering based on 
// Henyey-Greenstein phase function
// Forward and backward angles are treated separately.
//
// mail:  gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4OpMieHG.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

G4OpMieHG::G4OpMieHG(const G4String& processName, G4ProcessType type)
           : G4VDiscreteProcess(processName, type)
{
        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        SetProcessSubType(fOpMieHG);
}

G4OpMieHG::~G4OpMieHG(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// -------------
//
G4VParticleChange* 
G4OpMieHG::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
	const G4Material* aMaterial = aTrack.GetMaterial();
	G4MaterialPropertiesTable* aMaterialPropertyTable =
	  aMaterial->GetMaterialPropertiesTable();

	G4double forward_g =
              aMaterialPropertyTable->GetConstProperty(kMIEHG_FORWARD);
	G4double backward_g =
              aMaterialPropertyTable->GetConstProperty(kMIEHG_BACKWARD);
	G4double ForwardRatio =
              aMaterialPropertyTable->GetConstProperty(kMIEHG_FORWARD_RATIO);

        if (verboseLevel>0) {
		G4cout << "MIE Scattering Photon!" << G4endl;
		G4cout << "MIE Old Momentum Direction: "
	     	     << aParticle->GetMomentumDirection() << G4endl;
		G4cout << "MIE Old Polarization: "
		     << aParticle->GetPolarization() << G4endl;
	}

        G4double gg;
        G4int direction;
        if (G4UniformRand()<=ForwardRatio){
           gg = forward_g;
           direction = 1;
        } else {
           gg = backward_g;
           direction = -1;
	}

        G4double r = G4UniformRand();

        G4double Theta;
        //sample the direction
        if (gg!=0) {
          Theta = std::acos(2*r*(1+gg)*(1+gg)*(1-gg+gg*r)/((1-gg+2*gg*r)*(1-gg+2*gg*r)) -1);
        } else {
          Theta = std::acos(2*r-1.);
	}
        G4double Phi = G4UniformRand()*2*pi;

        if (direction==-1) Theta = pi - Theta; //backward scattering

        G4ThreeVector NewMomentumDirection, OldMomentumDirection;
        G4ThreeVector OldPolarization, NewPolarization;

        NewMomentumDirection.set
                       (std::sin(Theta)*std::cos(Phi), std::sin(Theta)*std::sin(Phi), std::cos(Theta));
        OldMomentumDirection = aParticle->GetMomentumDirection();
        NewMomentumDirection.rotateUz(OldMomentumDirection);
        NewMomentumDirection = NewMomentumDirection.unit();

        OldPolarization = aParticle->GetPolarization();
        G4double constant = -1./NewMomentumDirection.dot(OldPolarization);

        NewPolarization = NewMomentumDirection + constant*OldPolarization;
        NewPolarization = NewPolarization.unit();

        if (NewPolarization.mag()==0) {
           r = G4UniformRand()*twopi;
           NewPolarization.set(std::cos(r),std::sin(r),0.);
           NewPolarization.rotateUz(NewMomentumDirection);
        } else {
	   // There are two directions which perpendicular
           // new momentum direction
           if (G4UniformRand() < 0.5) NewPolarization = -NewPolarization;
        }

        aParticleChange.ProposePolarization(NewPolarization);
        aParticleChange.ProposeMomentumDirection(NewMomentumDirection);

        if (verboseLevel>0) {
              G4cout << "MIE New Polarization: " 
                     << NewPolarization << G4endl;
              G4cout << "MIE Polarization Change: "
                     << *(aParticleChange.GetPolarization()) << G4endl;  
              G4cout << "MIE New Momentum Direction: " 
                     << NewMomentumDirection << G4endl;
              G4cout << "MIE Momentum Change: "
                     << *(aParticleChange.GetMomentumDirection()) << G4endl; 
        }

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath()
// -----------------
//
G4double G4OpMieHG::GetMeanFreePath(const G4Track& aTrack,
                                    G4double ,
                                    G4ForceCondition* )
{
        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

        G4double thePhotonEnergy = aParticle->GetTotalEnergy();

        G4double AttenuationLength = DBL_MAX;

        G4MaterialPropertiesTable* aMaterialPropertyTable =
          aMaterial->GetMaterialPropertiesTable();

        if (aMaterialPropertyTable) {
           G4MaterialPropertyVector* AttenuationLengthVector =
                                 aMaterialPropertyTable->GetProperty("MIEHG");
           if (AttenuationLengthVector) {
              AttenuationLength = AttenuationLengthVector ->
                                    Value(thePhotonEnergy);
           } else {
//              G4cout << "No Mie scattering length specified" << G4endl;
           }
        } else {
//           G4cout << "No Mie scattering length specified" << G4endl; 
        }

//	G4cout << thePhotonEnergy/GeV << " \t" << AttenuationLength/m << G4endl;

        return AttenuationLength;
}
