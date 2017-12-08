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
// $Id: G4OpRayleigh.cc 106117 2017-09-13 10:23:20Z gcosmo $
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
// Updated:     2014-10-10 -  This version calculates the Rayleigh scattering   
//              length for more materials than just Water (although the Water
//              default is kept). To do this the user would need to specify the
//              ISOTHERMAL_COMPRESSIBILITY as a material property and
//              optionally an RS_SCALE_LENGTH (useful for testing). Code comes
//              from Philip Graham (Queen Mary University of London).
//              2010-06-11 - Fix Bug 207; Thanks to Xin Qian
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

#include "G4OpRayleigh.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"

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

        thePhysicsTable = NULL;

        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }
}

// G4OpRayleigh::G4OpRayleigh(const G4OpRayleigh &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4OpRayleigh::~G4OpRayleigh()
{
        if (thePhysicsTable) {
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

        G4double rand, constant;
        G4double CosTheta, SinTheta, SinPhi, CosPhi, unit_x, unit_y, unit_z;

        do {
           // Try to simulate the scattered photon momentum direction
           // w.r.t. the initial photon momentum direction

           CosTheta = G4UniformRand();
           SinTheta = std::sqrt(1.-CosTheta*CosTheta);
           // consider for the angle 90-180 degrees
           if (G4UniformRand() < 0.5) CosTheta = -CosTheta;

           // simulate the phi angle
           rand = twopi*G4UniformRand();
           SinPhi = std::sin(rand);
           CosPhi = std::cos(rand);

           // start constructing the new momentum direction
	   unit_x = SinTheta * CosPhi; 
	   unit_y = SinTheta * SinPhi;  
	   unit_z = CosTheta; 
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
           constant = -NewMomentumDirection.dot(OldPolarization);

           NewPolarization = OldPolarization + constant*NewMomentumDirection;
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
          // Loop checking, 13-Aug-2015, Peter Gumplinger
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

// BuildPhysicsTable for the Rayleigh Scattering process
// --------------------------------------------------------
void G4OpRayleigh::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (thePhysicsTable) {
     thePhysicsTable->clearAndDestroy();
     delete thePhysicsTable;
     thePhysicsTable = NULL;
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  thePhysicsTable = new G4PhysicsTable( numOfMaterials );
  
  for( G4int iMaterial = 0; iMaterial < numOfMaterials; iMaterial++ )
  {
      G4Material* material = (*theMaterialTable)[iMaterial];
      G4MaterialPropertiesTable* materialProperties = 
                                       material->GetMaterialPropertiesTable();
      G4PhysicsOrderedFreeVector* rayleigh = NULL;
      if ( materialProperties != NULL ) {
         rayleigh = materialProperties->GetProperty( kRAYLEIGH );
         if ( rayleigh == NULL ) rayleigh = 
                                   CalculateRayleighMeanFreePaths( material );
      }
      thePhysicsTable->insertAt( iMaterial, rayleigh );
  }
}

// GetMeanFreePath()
// -----------------
//
G4double G4OpRayleigh::GetMeanFreePath(const G4Track& aTrack,
                                       G4double ,
                                       G4ForceCondition* )
{
  const G4DynamicParticle* particle = aTrack.GetDynamicParticle();
  const G4double photonMomentum = particle->GetTotalMomentum();
  const G4Material* material = aTrack.GetMaterial();

  G4PhysicsOrderedFreeVector* rayleigh = 
                              static_cast<G4PhysicsOrderedFreeVector*>
                              ((*thePhysicsTable)(material->GetIndex()));
  
  G4double rsLength = DBL_MAX;
  if( rayleigh != NULL ) rsLength = rayleigh->Value( photonMomentum );
  return rsLength;
}

// CalculateRayleighMeanFreePaths()
// --------------------------------
// Private method to compute Rayleigh Scattering Lengths
G4PhysicsOrderedFreeVector* 
G4OpRayleigh::CalculateRayleighMeanFreePaths( const G4Material* material ) const
{
  G4MaterialPropertiesTable* materialProperties = 
                                       material->GetMaterialPropertiesTable();

  // Retrieve the beta_T or isothermal compressibility value. For backwards
  // compatibility use a constant if the material is "Water". If the material
  // doesn't have an ISOTHERMAL_COMPRESSIBILITY constant then return
  G4double betat;
  if ( material->GetName() == "Water" )
    betat = 7.658e-23*m3/MeV;
  else if(materialProperties->ConstPropertyExists("ISOTHERMAL_COMPRESSIBILITY"))
    betat = materialProperties->GetConstProperty(kISOTHERMAL_COMPRESSIBILITY);
  else
    return NULL;

  // If the material doesn't have a RINDEX property vector then return
  G4MaterialPropertyVector* rIndex = materialProperties->GetProperty(kRINDEX);
  if ( rIndex == NULL ) return NULL;

  // Retrieve the optional scale factor, (this just scales the scattering length
  G4double scaleFactor = 1.0;
  if( materialProperties->ConstPropertyExists( "RS_SCALE_FACTOR" ) )
    scaleFactor= materialProperties->GetConstProperty(kRS_SCALE_FACTOR );

  // Retrieve the material temperature. For backwards compatibility use a 
  // constant if the material is "Water"
  G4double temperature;
  if( material->GetName() == "Water" )
    temperature = 283.15*kelvin; // Temperature of water is 10 degrees celsius
  else
    temperature = material->GetTemperature();

  G4PhysicsOrderedFreeVector* rayleighMeanFreePaths =
                                             new G4PhysicsOrderedFreeVector();
  // This calculates the meanFreePath via the Einstein-Smoluchowski formula
  const G4double c1 = scaleFactor * betat * temperature * k_Boltzmann / 
                      ( 6.0 * pi );

  for( size_t uRIndex = 0; uRIndex < rIndex->GetVectorLength(); uRIndex++ )
  {
     const G4double energy = rIndex->Energy( uRIndex );
     const G4double rIndexSquared = (*rIndex)[uRIndex] * (*rIndex)[uRIndex];
     const G4double xlambda = h_Planck * c_light / energy;
     const G4double c2 = std::pow(twopi/xlambda,4);
     const G4double c3 = 
                    std::pow(((rIndexSquared-1.0)*(rIndexSquared+2.0 )/3.0),2);

     const G4double meanFreePath = 1.0 / ( c1 * c2 * c3 );

     if( verboseLevel>0 )
       G4cout << energy << "MeV\t" << meanFreePath << "mm" << G4endl;

     rayleighMeanFreePaths->InsertValues( energy, meanFreePath );
  }

  return rayleighMeanFreePaths;
}
