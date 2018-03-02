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
// $Id: G4Cerenkov.cc 108508 2018-02-15 15:54:35Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4Cerenkov.cc
// Description: Discrete Process -- Generation of Cerenkov Photons
// Version:     2.1
// Created:     1996-02-21
// Author:      Juliet Armstrong
// Updated:     2007-09-30 by Peter Gumplinger
//              > change inheritance to G4VDiscreteProcess
//              GetContinuousStepLimit -> GetMeanFreePath (StronglyForced)
//              AlongStepDoIt -> PostStepDoIt
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2001-09-17, migration of Materials to pure STL (mma)
//              2000-11-12 by Peter Gumplinger
//              > add check on CerenkovAngleIntegrals->IsFilledVectorExist()
//              in method GetAverageNumberOfPhotons
//              > and a test for MeanNumberOfPhotons <= 0.0 in DoIt
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              1999-10-29 by Peter Gumplinger
//              > change: == into <= in GetContinuousStepLimit
//              1997-08-08 by Peter Gumplinger
//              > add protection against /0
//              > G4MaterialPropertiesTable; new physics/tracking scheme
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"
#include "G4EmProcessSubType.hh"

#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"

#include "G4Cerenkov.hh"

/////////////////////////
// Class Implementation  
/////////////////////////

//G4bool G4Cerenkov::fTrackSecondariesFirst = false;
//G4double G4Cerenkov::fMaxBetaChange = 0.;
//G4int G4Cerenkov::fMaxPhotons = 0;

  //////////////
  // Operators
  //////////////

// G4Cerenkov::operator=(const G4Cerenkov &right)
// {
// }

  /////////////////
  // Constructors
  /////////////////

G4Cerenkov::G4Cerenkov(const G4String& processName, G4ProcessType type)
           : G4VProcess(processName, type),
             fTrackSecondariesFirst(false),
             fMaxBetaChange(0.0),
             fMaxPhotons(0),
             fStackingFlag(true),
             fNumPhotons(0)
{
  SetProcessSubType(fCerenkov);

  thePhysicsTable = nullptr;

  if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created " << G4endl;
  }
}

// G4Cerenkov::G4Cerenkov(const G4Cerenkov &right)
// {
// }

  ////////////////
  // Destructors
  ////////////////

G4Cerenkov::~G4Cerenkov()
{
  if (thePhysicsTable != nullptr) {
     thePhysicsTable->clearAndDestroy();
     delete thePhysicsTable;
  }
}

  ////////////
  // Methods
  ////////////

G4bool G4Cerenkov::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  G4bool result = false;

  if (aParticleType.GetPDGCharge() != 0.0 &&
      aParticleType.GetPDGMass() != 0.0 &&
      aParticleType.GetParticleName() != "chargedgeantino" &&
      !aParticleType.IsShortLived() ) { result = true; }

  return result;
}

void G4Cerenkov::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
}

void G4Cerenkov::SetMaxBetaChangePerStep(const G4double value)
{
  fMaxBetaChange = value*CLHEP::perCent;
}

void G4Cerenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons)
{
  fMaxPhotons = NumPhotons;
}

void G4Cerenkov::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (!thePhysicsTable) BuildThePhysicsTable();
}

// PostStepDoIt
// -------------
//
G4VParticleChange*
G4Cerenkov::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The
// parameters are then transformed into the Master Reference System, and
// they are added to the particle change.

{
  ////////////////////////////////////////////////////
  // Should we ensure that the material is dispersive?
  ////////////////////////////////////////////////////

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial = aTrack.GetMaterial();

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0 = pPreStepPoint->GetGlobalTime();

  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable) return pParticleChange;

  G4MaterialPropertyVector* Rindex = 
                aMaterialPropertiesTable->GetProperty(kRINDEX); 
  if (!Rindex) return pParticleChange;

  // particle charge
  G4double charge = aParticle->GetDefinition()->GetPDGCharge();

  // particle beta
  G4double beta = (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta())*0.5;

  fNumPhotons = 0;

  G4double MeanNumberOfPhotons = 
                     GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);

  if (MeanNumberOfPhotons <= 0.0) {

     // return unchanged particle and no secondaries

     aParticleChange.SetNumberOfSecondaries(0);
 
     return pParticleChange;

  }

  G4double step_length = aStep.GetStepLength();

  MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;

  fNumPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);

  if ( fNumPhotons <= 0 || !fStackingFlag ) {

     // return unchanged particle and no secondaries  

     aParticleChange.SetNumberOfSecondaries(0);

     return pParticleChange;

  }

  ////////////////////////////////////////////////////////////////

  aParticleChange.SetNumberOfSecondaries(fNumPhotons);

  if (fTrackSecondariesFirst) {
     if (aTrack.GetTrackStatus() == fAlive )
                           aParticleChange.ProposeTrackStatus(fSuspend);
  }

  ////////////////////////////////////////////////////////////////

  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();
  G4double dp = Pmax - Pmin;

  G4double nMax = Rindex->GetMaxValue();

  G4double BetaInverse = 1./beta;

  G4double maxCos = BetaInverse / nMax; 
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

  G4double beta1 = pPreStepPoint ->GetBeta();
  G4double beta2 = pPostStepPoint->GetBeta();

  G4double MeanNumberOfPhotons1 =
                     GetAverageNumberOfPhotons(charge,beta1,aMaterial,Rindex);
  G4double MeanNumberOfPhotons2 =
                     GetAverageNumberOfPhotons(charge,beta2,aMaterial,Rindex);

  for (G4int i = 0; i < fNumPhotons; i++) {

      // Determine photon energy

      G4double rand;
      G4double sampledEnergy, sampledRI; 
      G4double cosTheta, sin2Theta;

      // sample an energy

      do {
         rand = G4UniformRand();	
         sampledEnergy = Pmin + rand * dp; 
         sampledRI = Rindex->Value(sampledEnergy);
         cosTheta = BetaInverse / sampledRI;  

         sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);
         rand = G4UniformRand();	

        // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
      } while (rand*maxSin2 > sin2Theta);

      // Generate random position of photon on cone surface 
      // defined by Theta 

      rand = G4UniformRand();

      G4double phi = twopi*rand;
      G4double sinPhi = std::sin(phi);
      G4double cosPhi = std::cos(phi);

      // calculate x,y, and z components of photon energy
      // (in coord system with primary particle direction 
      //  aligned with the z axis)

      G4double sinTheta = std::sqrt(sin2Theta); 
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
        new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);

      aCerenkovPhoton->SetPolarization(photonPolarization.x(),
                                       photonPolarization.y(),
                                       photonPolarization.z());

      aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

      // Generate new G4Track object:

      G4double NumberOfPhotons, N;

      do {
         rand = G4UniformRand();
         NumberOfPhotons = MeanNumberOfPhotons1 - rand *
                                (MeanNumberOfPhotons1-MeanNumberOfPhotons2);
         N = G4UniformRand() *
                        std::max(MeanNumberOfPhotons1,MeanNumberOfPhotons2);
        // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
      } while (N > NumberOfPhotons);

      G4double delta = rand * aStep.GetStepLength();

      G4double deltaTime = delta / (pPreStepPoint->GetVelocity()+
                                      rand*(pPostStepPoint->GetVelocity()-
                                            pPreStepPoint->GetVelocity())*0.5);

      G4double aSecondaryTime = t0 + deltaTime;

      G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

      G4Track* aSecondaryTrack = 
               new G4Track(aCerenkovPhoton,aSecondaryTime,aSecondaryPosition);

      aSecondaryTrack->SetTouchableHandle(
                               aStep.GetPreStepPoint()->GetTouchableHandle());

      aSecondaryTrack->SetParentID(aTrack.GetTrackID());

      aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel>0) {
     G4cout <<"\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = "
	    << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return pParticleChange;
}

// BuildThePhysicsTable for the Cerenkov process
// ---------------------------------------------
//

void G4Cerenkov::BuildThePhysicsTable()
{
  if (thePhysicsTable) return;

  const G4MaterialTable* theMaterialTable=
  G4Material::GetMaterialTable();
  G4int numOfMaterials = G4Material::GetNumberOfMaterials();

  // create new physics table
	
  thePhysicsTable = new G4PhysicsTable(numOfMaterials);

  // loop for materials

  for (G4int i=0 ; i < numOfMaterials; i++) {

      G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector = 0;

      // Retrieve vector of refraction indices for the material
      // from the material's optical properties table 

      G4Material* aMaterial = (*theMaterialTable)[i];

      G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                      aMaterial->GetMaterialPropertiesTable();

      if (aMaterialPropertiesTable) {
         aPhysicsOrderedFreeVector = new G4PhysicsOrderedFreeVector();
         G4MaterialPropertyVector* theRefractionIndexVector = 
                              aMaterialPropertiesTable->GetProperty(kRINDEX);

         if (theRefractionIndexVector) {

            // Retrieve the first refraction index in vector
            // of (photon energy, refraction index) pairs 

            G4double currentRI = (*theRefractionIndexVector)[0];

            if (currentRI > 1.0) {

               // Create first (photon energy, Cerenkov Integral)
               // pair  

               G4double currentPM = theRefractionIndexVector->Energy(0);
               G4double currentCAI = 0.0;

               aPhysicsOrderedFreeVector->InsertValues(currentPM , currentCAI);

               // Set previous values to current ones prior to loop

               G4double prevPM  = currentPM;
               G4double prevCAI = currentCAI;
               G4double prevRI  = currentRI;

               // loop over all (photon energy, refraction index)
               // pairs stored for this material  

               for (size_t ii = 1;
                           ii < theRefractionIndexVector->GetVectorLength();
                           ++ii) {
                   currentRI = (*theRefractionIndexVector)[ii];
                   currentPM = theRefractionIndexVector->Energy(ii);

                   currentCAI = 0.5*(1.0/(prevRI*prevRI) +
                                     1.0/(currentRI*currentRI));

                   currentCAI = prevCAI + (currentPM - prevPM) * currentCAI;

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

// GetMeanFreePath
// ---------------
//

G4double G4Cerenkov::GetMeanFreePath(const G4Track&,
                                           G4double,
                                           G4ForceCondition*)
{
  return 1.;
}

G4double G4Cerenkov::PostStepGetPhysicalInteractionLength(
                                           const G4Track& aTrack,
                                           G4double,
                                           G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double StepLimit = DBL_MAX;

  const G4Material* aMaterial = aTrack.GetMaterial();
  G4int materialIndex = aMaterial->GetIndex();

  // If Physics Vector is not defined no Cerenkov photons
  //    this check avoid string comparison below

  if(!(*thePhysicsTable)[materialIndex]) { return StepLimit; }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

  G4double kineticEnergy = aParticle->GetKineticEnergy();
  const G4ParticleDefinition* particleType = aParticle->GetDefinition();
  G4double mass = particleType->GetPDGMass();

  // particle beta
  G4double beta = aParticle->GetTotalMomentum() /
                  aParticle->GetTotalEnergy();
  // particle gamma
  G4double gamma = aParticle->GetTotalEnergy()/mass;

  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                      aMaterial->GetMaterialPropertiesTable();

  G4MaterialPropertyVector* Rindex = NULL;

  if (aMaterialPropertiesTable)
                     Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);

  G4double nMax;
  if (Rindex) {
     nMax = Rindex->GetMaxValue();
  } else {
     return StepLimit;
  }

  G4double BetaMin = 1./nMax;
  if ( BetaMin >= 1. ) return StepLimit;

  G4double GammaMin = 1./std::sqrt(1.-BetaMin*BetaMin);

  if (gamma < GammaMin ) return StepLimit;

  G4double kinEmin = mass*(GammaMin-1.);

  G4double RangeMin = G4LossTableManager::Instance()->GetRange(particleType,
                                                               kinEmin,
                                                               couple);
  G4double Range    = G4LossTableManager::Instance()->GetRange(particleType,
                                                               kineticEnergy,
                                                               couple);

  G4double Step = Range - RangeMin;
//  if (Step < 1.*um ) return StepLimit;

  if (Step > 0. && Step < StepLimit) StepLimit = Step; 

  // If user has defined an average maximum number of photons to
  // be generated in a Step, then calculate the Step length for
  // that number of photons. 
 
  if (fMaxPhotons > 0) {

     // particle charge
     const G4double charge = aParticle->GetDefinition()->GetPDGCharge();

     G4double MeanNumberOfPhotons = 
                      GetAverageNumberOfPhotons(charge,beta,aMaterial,Rindex);

     Step = 0.;
     if (MeanNumberOfPhotons > 0.0) Step = fMaxPhotons / MeanNumberOfPhotons;

     if (Step > 0. && Step < StepLimit) StepLimit = Step;
  }

  // If user has defined an maximum allowed change in beta per step
  if (fMaxBetaChange > 0.) {

     G4double dedx = G4LossTableManager::Instance()->GetDEDX(particleType,
                                                             kineticEnergy,
                                                             couple);

     G4double deltaGamma = gamma - 1./std::sqrt(1.-beta*beta*
                                                (1.-fMaxBetaChange)*
                                                (1.-fMaxBetaChange));

     Step = mass * deltaGamma / dedx;

     if (Step > 0. && Step < StepLimit) StepLimit = Step;

  }

  *condition = StronglyForced;
  return StepLimit;
}

// GetAverageNumberOfPhotons
// -------------------------
// This routine computes the number of Cerenkov photons produced per
// GEANT-unit (millimeter) in the current medium.
//             ^^^^^^^^^^

G4double
  G4Cerenkov::GetAverageNumberOfPhotons(const G4double charge,
                                        const G4double beta, 
			                const G4Material* aMaterial,
			                G4MaterialPropertyVector* Rindex) const
{
  const G4double Rfact = 369.81/(eV * cm);

  if(beta <= 0.0)return 0.0;

  G4double BetaInverse = 1./beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  // 	- Refraction Indices for the current material
  //	- new G4PhysicsOrderedFreeVector allocated to hold CAI's
 
  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material  

  G4PhysicsOrderedFreeVector* CerenkovAngleIntegrals =
             (G4PhysicsOrderedFreeVector*)((*thePhysicsTable)(materialIndex));

  if(!(CerenkovAngleIntegrals->IsFilledVectorExist()))return 0.0;

  // Min and Max photon energies 
  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();

  // Min and Max Refraction Indices 
  G4double nMin = Rindex->GetMinValue();	
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral 
  G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

  G4double dp, ge;

  // If n(Pmax) < 1/Beta -- no photons generated 

  if (nMax < BetaInverse) {
     dp = 0.0;
     ge = 0.0;
  } 

  // otherwise if n(Pmin) >= 1/Beta -- photons generated  

  else if (nMin > BetaInverse) {
     dp = Pmax - Pmin;	
     ge = CAImax; 
  } 

  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
  // we need to find a P such that the value of n(P) == 1/Beta.
  // Interpolation is performed by the GetEnergy() and
  // Value() methods of the G4MaterialPropertiesTable and
  // the GetValue() method of G4PhysicsVector.  

  else {
     Pmin = Rindex->GetEnergy(BetaInverse);
     dp = Pmax - Pmin;

     // need boolean for current implementation of G4PhysicsVector
     // ==> being phased out
     G4bool isOutRange;
     G4double CAImin = CerenkovAngleIntegrals->GetValue(Pmin, isOutRange);
     ge = CAImax - CAImin;

     if (verboseLevel>0) {
        G4cout << "CAImin = " << CAImin << G4endl;
        G4cout << "ge = " << ge << G4endl;
     }
  }
	
  // Calculate number of photons 
  G4double NumPhotons = Rfact * charge/eplus * charge/eplus *
                                 (dp - ge * BetaInverse*BetaInverse);

  return NumPhotons;		
}

void G4Cerenkov::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = thePhysicsTable->entries();
  G4PhysicsOrderedFreeVector *v;

  for (G4int i = 0 ; i < PhysicsTableSize ; i++ ) {
      v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
      v->DumpValues();
  }
}

