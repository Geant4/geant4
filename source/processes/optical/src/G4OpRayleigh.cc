// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpRayleigh.cc,v 1.5 2001-01-31 04:08:12 gum Exp $
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
// Updated:     2001-01-30 by Peter Gumplinger
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

G4OpRayleigh::G4OpRayleigh(const G4String& processName)
           : G4VDiscreteProcess(processName)
{

        thePhysicsTable = NULL;

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
        if (thePhysicsTable!= NULL) {
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
        const G4Material* aMaterial = aTrack.GetMaterial();

        if (verboseLevel>0) {
		G4cout << "Scattering Photon!" << G4endl;
		G4cout << "Old Momentum Direction: "
	     	     << aParticle->GetMomentumDirection() << G4endl;
		G4cout << "Old Polarization: "
		     << aParticle->GetPolarization() << G4endl;
	}

	// find polar angle w.r.t. old polarization vector

	G4double rand = G4UniformRand();

	G4double CosTheta = pow(rand, 1./3.);
	G4double SinTheta = sqrt(1.-CosTheta*CosTheta);

        if(G4UniformRand() < 0.5)CosTheta = -CosTheta;

	// find azimuthal angle w.r.t old polarization vector 

	rand = G4UniformRand();

	G4double Phi = twopi*rand;
	G4double SinPhi = sin(Phi); 
	G4double CosPhi = cos(Phi); 
	
	G4double unit_x = SinTheta * CosPhi; 
	G4double unit_y = SinTheta * SinPhi;  
	G4double unit_z = CosTheta; 
	
        G4ThreeVector NewPolarization (unit_x,unit_y,unit_z);

        // Rotate new polarization direction into global reference system 

	G4ThreeVector OldPolarization = aParticle->GetPolarization();
        OldPolarization = OldPolarization.unit();

	NewPolarization.rotateUz(OldPolarization);
        NewPolarization = NewPolarization.unit();
	
        // -- new momentum direction is normal to the new
        // polarization vector and in the same plane as the
        // old and new polarization vectors --

        G4ThreeVector NewMomentumDirection = 
                              OldPolarization - NewPolarization * CosTheta;

        if(G4UniformRand() < 0.5)NewMomentumDirection = -NewMomentumDirection;
        NewMomentumDirection = NewMomentumDirection.unit();

	aParticleChange.SetPolarizationChange(NewPolarization);

	aParticleChange.SetMomentumChange(NewMomentumDirection);

        if (verboseLevel>0) {
		G4cout << "New Polarization: " 
		     << NewPolarization << G4endl;
		G4cout << "Polarization Change: "
		     << *(aParticleChange.GetPolarizationChange()) << G4endl;  
		G4cout << "New Momentum Direction: " 
		     << NewMomentumDirection << G4endl;
		G4cout << "Momentum Change: "
		     << *(aParticleChange.GetMomentumChange()) << G4endl; 
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
        G4int numOfMaterials = theMaterialTable->length();

        // create a new physics table

        thePhysicsTable = new G4PhysicsTable(numOfMaterials);

        // loop for materials

        for (G4int i=0 ; i < numOfMaterials; i++)
        {
                G4PhysicsOrderedFreeVector* ScatteringLengths =
                                new G4PhysicsOrderedFreeVector();

                if ((*theMaterialTable)[i]->GetName() == "Water")
                {
			G4MaterialPropertiesTable *MaterialPT =
			(*theMaterialTable)[i]->GetMaterialPropertiesTable();
			// Call utility routine to Generate
			// Rayleigh Scattering Lengths
			ScatteringLengths =
				RayleighAttenuationLengthGenerator(MaterialPT);
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

        G4double thePhotonMomentum = aParticle->GetTotalMomentum();

        G4double AttenuationLength = DBL_MAX;

        if (aMaterial->GetName() == "Water") {

           G4bool isOutRange;

           AttenuationLength =
                (*thePhysicsTable)(aMaterial->GetIndex())->
                           GetValue(thePhotonMomentum, isOutRange);
        }
        else {

           G4MaterialPropertiesTable* aMaterialPropertyTable =
                           aMaterial->GetMaterialPropertiesTable();

           if(aMaterialPropertyTable){
             G4MaterialPropertyVector* AttenuationLengthVector =
                   aMaterialPropertyTable->GetProperty("RAYLEIGH");
             if(AttenuationLengthVector){
               AttenuationLength = AttenuationLengthVector ->
                                    GetProperty(thePhotonMomentum);
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
        // and photon momentum from the material properties table

        G4MaterialPropertyVector* Rindex = aMPT->GetProperty("RINDEX");

        G4double refsq;
        G4double e;
        G4double xlambda;
        G4double c1, c2, c3, c4;
        G4double Dist;
        G4double refraction_index;

        G4double no_unit = 1.0;

        G4PhysicsOrderedFreeVector *RayleighScatteringLengths = 
				new G4PhysicsOrderedFreeVector();
        Rindex->ResetIterator();

        while (++(*Rindex)) {

                e = (Rindex->GetPhotonMomentum());

                refraction_index = Rindex->GetProperty();
                refsq = refraction_index*refraction_index;
                xlambda = h_Planck*c_light/e;

	        if (verboseLevel>0) {
        	        G4cout << Rindex->GetPhotonMomentum() << " MeV\t";
                	G4cout << xlambda << " mm\t";
		}

                c1 = 1 / (6.0 * pi);
                c2 = pow((2.0 * pi / xlambda), 4);
                c3 = pow( ( (refsq - 1.0) * (refsq + 2.0) / 3.0 ), 2);
                c4 = betat * temp * kboltz;

                Dist = 1.0 / (c1*c2*c3*c4);

	        if (verboseLevel>0) {
	                G4cout << Dist << " mm" << G4endl;
		}
                RayleighScatteringLengths->
			InsertValues(Rindex->GetPhotonMomentum(), Dist);
        }

	return RayleighScatteringLengths;
}
