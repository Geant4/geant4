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
// $Id: G4Scintillation.cc,v 1.38 2010-12-15 07:39:26 gunter Exp $
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
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
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
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4EmProcessSubType.hh"

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

G4Scintillation::G4Scintillation(const G4String& processName,
                                       G4ProcessType type)
                  : G4VRestDiscreteProcess(processName, type)
{
        SetProcessSubType(fScintillation);

        fTrackSecondariesFirst = false;
        fFiniteRiseTime = false;

        YieldFactor = 1.0;
        ExcitationRatio = 1.0;

        scintillationByParticleType = false;

        theFastIntegralTable = NULL;
        theSlowIntegralTable = NULL;

        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        BuildThePhysicsTable();

        emSaturation = NULL;
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

        G4double ScintillationYield = 0.;

        if (scintillationByParticleType) {
           // The scintillation response is a function of the energy
           // deposited by particle types.

           // Get the definition of the current particle
           G4ParticleDefinition *pDef = aParticle->GetDefinition();
           const G4MaterialPropertyVector *Scint_Yield_Vector = NULL;

           // Obtain the G4MaterialPropertyVectory containing the
           // scintillation light yield as a function of the deposited
           // energy for the current particle type

           // Protons
           if(pDef==G4Proton::ProtonDefinition()) 
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("PROTONSCINTILLATIONYIELD");

           // Deuterons
           else if(pDef==G4Deuteron::DeuteronDefinition())
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("DEUTERONSCINTILLATIONYIELD");

           // Tritons
           else if(pDef==G4Triton::TritonDefinition())
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("TRITONSCINTILLATIONYIELD");

           // Alphas
           else if(pDef==G4Alpha::AlphaDefinition())
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("ALPHASCINTILLATIONYIELD");
	  
           // Ions (particles derived from G4VIon and G4Ions)
           // and recoil ions below tracking cut from neutrons after hElastic
           else if(pDef->GetParticleType()== "nucleus" || 
                   pDef==G4Neutron::NeutronDefinition()) 
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("IONSCINTILLATIONYIELD");

           // Electrons (must also account for shell-binding energy
           // attributed to gamma from standard PhotoElectricEffect)
           else if(pDef==G4Electron::ElectronDefinition() ||
                   pDef==G4Gamma::GammaDefinition())
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("ELECTRONSCINTILLATIONYIELD");

           // Default for particles not enumerated/listed above
           else
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("ELECTRONSCINTILLATIONYIELD");
	   
           // If the user has not specified yields for (p,d,t,a,carbon)
           // then these unspecified particles will default to the 
           // electron's scintillation yield
           if(!Scint_Yield_Vector){
             Scint_Yield_Vector = aMaterialPropertiesTable->
               GetProperty("ELECTRONSCINTILLATIONYIELD");
           }

           // Throw an exception if no scintillation yield is found
           if (!Scint_Yield_Vector) {
              G4cerr << "\nG4Scintillation::PostStepDoIt(): "
                     << "Request for scintillation yield for energy deposit and particle type without correct entry in MaterialPropertiesTable\n"
                     << "ScintillationByParticleType requires at minimum that ELECTRONSCINTILLATIONYIELD is set by the user\n"
                     << G4endl;
             G4Exception("G4Scintillation::PostStepDoIt",
                         "No correct entry in MaterialPropertiesTable",
                         FatalException,"Missing MaterialPropertiesTable entry.");
             return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
           }

           if (verboseLevel>1) {
             G4cout << "\n"
                    << "Particle = " << pDef->GetParticleName() << "\n"
                    << "Energy Dep. = " << TotalEnergyDeposit/MeV << "\n"
                    << "Yield = " 
                    << Scint_Yield_Vector->GetProperty(TotalEnergyDeposit) 
                    << "\n" << G4endl;
           }

           // Obtain the scintillation yield using the total energy
           // deposited by the particle in this step.

           // Units: [# scintillation photons]
           ScintillationYield = Scint_Yield_Vector->
                                            GetProperty(TotalEnergyDeposit);
        } else {
           // The default linear scintillation process
           ScintillationYield = aMaterialPropertiesTable->
                                      GetConstProperty("SCINTILLATIONYIELD");

           // Units: [# scintillation photons / MeV]
           ScintillationYield *= YieldFactor;
        }

        G4double ResolutionScale    = aMaterialPropertiesTable->
                                      GetConstProperty("RESOLUTIONSCALE");

        // Birks law saturation:

        //G4double constBirks = 0.0;

        //constBirks = aMaterial->GetIonisation()->GetBirksConstant();

        G4double MeanNumberOfPhotons;

        // Birk's correction via emSaturation and specifying scintillation by
        // by particle type are physically mutually exclusive

        if (scintillationByParticleType)
           MeanNumberOfPhotons = ScintillationYield;
        else if (emSaturation)
           MeanNumberOfPhotons = ScintillationYield*
                              (emSaturation->VisibleEnergyDeposition(&aStep));
        else
           MeanNumberOfPhotons = ScintillationYield*TotalEnergyDeposit;

        G4int NumPhotons;

        if (MeanNumberOfPhotons > 10.)
        {
          G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
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

        if (fTrackSecondariesFirst) {
           if (aTrack.GetTrackStatus() == fAlive )
                  aParticleChange.ProposeTrackStatus(fSuspend);
        }

        ////////////////////////////////////////////////////////////////

        G4int materialIndex = aMaterial->GetIndex();

        // Retrieve the Scintillation Integral for this material  
        // new G4PhysicsOrderedFreeVector allocated to hold CII's

        G4int Num = NumPhotons;

        for (G4int scnt = 1; scnt <= nscnt; scnt++) {

            G4double ScintillationTime = 0.*ns;
            G4double ScintillationRiseTime = 0.*ns;
            G4PhysicsOrderedFreeVector* ScintillationIntegral = NULL;

            if (scnt == 1) {
               if (nscnt == 1) {
                 if(Fast_Intensity){
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("FASTTIMECONSTANT");
                   if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                  GetConstProperty("FASTSCINTILLATIONRISETIME");
                   }
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
                 }
                 if(Slow_Intensity){
                   ScintillationTime   = aMaterialPropertiesTable->
                                           GetConstProperty("SLOWTIMECONSTANT");
                   if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                  GetConstProperty("SLOWSCINTILLATIONRISETIME");
                   }
                   ScintillationIntegral =
                   (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
                 }
               }
               else {
                 G4double YieldRatio = aMaterialPropertiesTable->
                                          GetConstProperty("YIELDRATIO");
                 if ( ExcitationRatio == 1.0 ) {
                    Num = G4int (std::min(YieldRatio,1.0) * NumPhotons);
                 }
                 else {
                    Num = G4int (std::min(ExcitationRatio,1.0) * NumPhotons);
                 }
                 ScintillationTime   = aMaterialPropertiesTable->
                                          GetConstProperty("FASTTIMECONSTANT");
                 if (fFiniteRiseTime) {
                      ScintillationRiseTime = aMaterialPropertiesTable->
                                 GetConstProperty("FASTSCINTILLATIONRISETIME");
                 }
                 ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theFastIntegralTable)(materialIndex));
               }
            }
            else {
               Num = NumPhotons - Num;
               ScintillationTime   =   aMaterialPropertiesTable->
                                          GetConstProperty("SLOWTIMECONSTANT");
               if (fFiniteRiseTime) {
                    ScintillationRiseTime = aMaterialPropertiesTable->
                                 GetConstProperty("SLOWSCINTILLATIONRISETIME");
               }
               ScintillationIntegral =
                  (G4PhysicsOrderedFreeVector*)((*theSlowIntegralTable)(materialIndex));
            }

            if (!ScintillationIntegral) continue;
	
            // Max Scintillation Integral
 
            G4double CIImax = ScintillationIntegral->GetMaxValue();

            for (G4int i = 0; i < Num; i++) {

                // Determine photon energy

                G4double CIIvalue = G4UniformRand()*CIImax;
                G4double sampledEnergy = 
                              ScintillationIntegral->GetEnergy(CIIvalue);

                if (verboseLevel>1) {
                   G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
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

                G4DynamicParticle* aScintillationPhoton =
                  new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
                                                         photonMomentum);
                aScintillationPhoton->SetPolarization
                                     (photonPolarization.x(),
                                      photonPolarization.y(),
                                      photonPolarization.z());

                aScintillationPhoton->SetKineticEnergy(sampledEnergy);

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

                // emission time distribution
                if (ScintillationRiseTime==0.0) {
                   deltaTime = deltaTime - 
                          ScintillationTime * std::log( G4UniformRand() );
                } else {
                   deltaTime = deltaTime +
                          sample_time(ScintillationRiseTime, ScintillationTime);
                }

                G4double aSecondaryTime = t0 + deltaTime;

                G4ThreeVector aSecondaryPosition =
                                    x0 + rand * aStep.GetDeltaPosition();

                G4Track* aSecondaryTrack = 
                new G4Track(aScintillationPhoton,aSecondaryTime,aSecondaryPosition);

                aSecondaryTrack->SetTouchableHandle(
                                 aStep.GetPreStepPoint()->GetTouchableHandle());
                // aSecondaryTrack->SetTouchableHandle((G4VTouchable*)0);

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
                      // of (photon energy, intensity) pairs 

                      theFastLightVector->ResetIterator();
                      ++(*theFastLightVector);  // advance to 1st entry 

                      G4double currentIN = theFastLightVector->
                                               GetProperty();

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation 
                         // Integral pair  

                         G4double currentPM = theFastLightVector->
                                                  GetPhotonEnergy();

                         G4double currentCII = 0.0;

                         aPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material  

                         while(++(*theFastLightVector))
                         {
                                currentPM = theFastLightVector->
                                                GetPhotonEnergy();

                                currentIN = theFastLightVector->	
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
                      // of (photon energy, intensity) pairs

                      theSlowLightVector->ResetIterator();
                      ++(*theSlowLightVector);  // advance to 1st entry

                      G4double currentIN = theSlowLightVector->
                                           GetProperty();

                      if (currentIN >= 0.0) {

                         // Create first (photon energy, Scintillation
                         // Integral pair

                         G4double currentPM = theSlowLightVector->
                                                 GetPhotonEnergy();

                         G4double currentCII = 0.0;

                         bPhysicsOrderedFreeVector->
                                 InsertValues(currentPM , currentCII);

                         // Set previous values to current ones prior to loop

                         G4double prevPM  = currentPM;
                         G4double prevCII = currentCII;
                         G4double prevIN  = currentIN;

                         // loop over all (photon energy, intensity)
                         // pairs stored for this material

                         while(++(*theSlowLightVector))
                         {
                                currentPM = theSlowLightVector->
                                                GetPhotonEnergy();

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

// Called by the user to set the scintillation yield as a function
// of energy deposited by particle type

void G4Scintillation::SetScintillationByParticleType(const G4bool scintType)
{
        if (emSaturation) {
           G4Exception("G4Scintillation::SetScintillationByParticleType", "Redefinition",
                       JustWarning, "Birks Saturation is replaced by ScintillationByParticleType!");
           RemoveSaturation();
        }
        scintillationByParticleType = scintType;
}

// GetMeanFreePath
// ---------------
//

G4double G4Scintillation::GetMeanFreePath(const G4Track&,
                                          G4double ,
                                          G4ForceCondition* condition)
{
        *condition = StronglyForced;

        return DBL_MAX;

}

// GetMeanLifeTime
// ---------------
//

G4double G4Scintillation::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
        *condition = Forced;

        return DBL_MAX;

}

G4double G4Scintillation::sample_time(G4double tau1, G4double tau2)
{
// tau1: rise time and tau2: decay time

        while(1) {
          // two random numbers
          G4double ran1 = G4UniformRand();
          G4double ran2 = G4UniformRand();
          //
          // exponential distribution as envelope function: very efficient
          //
          G4double d = (tau1+tau2)/tau2;
          // make sure the envelope function is 
          // always larger than the bi-exponential
          G4double t = -1.0*tau2*std::log(1-ran1);
          G4double g = d*single_exp(t,tau2);
          if (ran2 <= bi_exp(t,tau1,tau2)/g) return t;
        }
        return -1.0;
}
