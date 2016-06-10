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
// $Id: G4PenelopeAnnihilationModel.cc 74452 2013-10-07 15:08:00Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// --------
// 29 Oct 2008   L Pandola    Migration from process to model 
// 15 Apr 2009   V Ivanchenko Cleanup initialisation and generation of 
//                    secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - do not use G4ElementSelector
// 02 Oct 2013   L.Pandola    Migration to MT

#include "G4PenelopeAnnihilationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeAnnihilationModel::fPielr2 = 0;

G4PenelopeAnnihilationModel::G4PenelopeAnnihilationModel(const G4ParticleDefinition* part,
                                             const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),fParticle(0),isInitialised(false)
{
  fIntrinsicLowEnergyLimit = 0.0;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);

  if (part)
    SetParticle(part);

  //Calculate variable that will be used later on
  fPielr2 = pi*classic_electr_radius*classic_electr_radius;

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeAnnihilationModel::~G4PenelopeAnnihilationModel()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeAnnihilationModel::Initialise(const G4ParticleDefinition* part,
					     const G4DataVector&)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeAnnihilationModel::Initialise()" << G4endl;
  SetParticle(part);

  if (IsMaster() && part == fParticle)
    {

      if(verboseLevel > 0) {
	G4cout << "Penelope Annihilation model is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / keV << " keV - "
	       << HighEnergyLimit() / GeV << " GeV"
	       << G4endl;
      }
    }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PenelopeAnnihilationModel::InitialiseLocal(const G4ParticleDefinition* part,
						  G4VEmModel* masterModel)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeAnnihilationModel::InitialiseLocal()" << G4endl;

  //
  //Check that particle matches: one might have multiple master models (e.g. 
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeAnnihilationModel* theModel = 
        static_cast<G4PenelopeAnnihilationModel*> (masterModel);
 
      //Same verbosity for all workers, as the master
      verboseLevel = theModel->verboseLevel;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeAnnihilationModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopeAnnihilationModel" << 
      G4endl;

  G4double cs = Z*ComputeCrossSectionPerElectron(energy);
  
  if (verboseLevel > 2)
    G4cout << "Annihilation cross Section at " << energy/keV << " keV for Z=" << Z << 
      " = " << cs/barn << " barn" << G4endl;
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeAnnihilationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple*,
					      const G4DynamicParticle* aDynamicPositron,
					      G4double,
					      G4double)
{
  //
  // Penelope model to sample final state for positron annihilation. 
  // Target eletrons are assumed to be free and at rest. Binding effects enabling 
  // one-photon annihilation are neglected.
  // For annihilation at rest, two back-to-back photons are emitted, having energy of 511 keV 
  // and isotropic angular distribution.
  // For annihilation in flight, it is used the theory from 
  //  W. Heitler, The quantum theory of radiation, Oxford University Press (1954)
  // The two photons can have different energy. The efficiency of the sampling algorithm 
  // of the photon energy from the dSigma/dE distribution is practically 100% for 
  // positrons of kinetic energy < 10 keV. It reaches a minimum (about 80%) at energy 
  // of about 10 MeV.
  // The angle theta is kinematically linked to the photon energy, to ensure momentum 
  // conservation. The angle phi is sampled isotropically for the first gamma.
  //
  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeAnnihilationModel" << G4endl;

  G4double kineticEnergy = aDynamicPositron->GetKineticEnergy();

  // kill primary
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  
  if (kineticEnergy == 0.0)
    {
      //Old AtRestDoIt
      G4double cosTheta = -1.0+2.0*G4UniformRand();
      G4double sinTheta = std::sqrt(1.0-cosTheta*cosTheta);
      G4double phi = twopi*G4UniformRand();
      G4ThreeVector direction (sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);
      G4DynamicParticle* firstGamma = new G4DynamicParticle (G4Gamma::Gamma(),
							     direction, electron_mass_c2);
      G4DynamicParticle* secondGamma = new G4DynamicParticle (G4Gamma::Gamma(),
							      -direction, electron_mass_c2);
  
      fvect->push_back(firstGamma);
      fvect->push_back(secondGamma);
      return;
    }

  //This is the "PostStep" case (annihilation in flight)
  G4ParticleMomentum positronDirection = 
    aDynamicPositron->GetMomentumDirection();
  G4double gamma = 1.0 + std::max(kineticEnergy,1.0*eV)/electron_mass_c2;
  G4double gamma21 = std::sqrt(gamma*gamma-1);
  G4double ani = 1.0+gamma;
  G4double chimin = 1.0/(ani+gamma21);
  G4double rchi = (1.0-chimin)/chimin;
  G4double gt0 = ani*ani-2.0;
  G4double test=0.0;
  G4double epsilon = 0;
  do{
    epsilon = chimin*std::pow(rchi,G4UniformRand());
    G4double reject = ani*ani*(1.0-epsilon)+2.0*gamma-(1.0/epsilon);
    test = G4UniformRand()*gt0-reject;
  }while(test>0);
   
  G4double totalAvailableEnergy = kineticEnergy + 2.0*electron_mass_c2;
  G4double photon1Energy = epsilon*totalAvailableEnergy;
  G4double photon2Energy = (1.0-epsilon)*totalAvailableEnergy;
  G4double cosTheta1 = (ani-1.0/epsilon)/gamma21;
  G4double cosTheta2 = (ani-1.0/(1.0-epsilon))/gamma21;
  
  //G4double localEnergyDeposit = 0.; 

  G4double sinTheta1 = std::sqrt(1.-cosTheta1*cosTheta1);
  G4double phi1  = twopi * G4UniformRand();
  G4double dirx1 = sinTheta1 * std::cos(phi1);
  G4double diry1 = sinTheta1 * std::sin(phi1);
  G4double dirz1 = cosTheta1;
  
  G4double sinTheta2 = std::sqrt(1.-cosTheta2*cosTheta2);
  G4double phi2  = phi1+pi;
  G4double dirx2 = sinTheta2 * std::cos(phi2);
  G4double diry2 = sinTheta2 * std::sin(phi2);
  G4double dirz2 = cosTheta2;
  
  G4ThreeVector photon1Direction (dirx1,diry1,dirz1);
  photon1Direction.rotateUz(positronDirection);   
  // create G4DynamicParticle object for the particle1  
  G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Gamma::Gamma(),
							   photon1Direction, 
							   photon1Energy);
  fvect->push_back(aParticle1);
 
  G4ThreeVector photon2Direction(dirx2,diry2,dirz2);
  photon2Direction.rotateUz(positronDirection); 
     // create G4DynamicParticle object for the particle2 
  G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
							   photon2Direction,
							   photon2Energy);
  fvect->push_back(aParticle2);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeAnnihilation" << G4endl;
      G4cout << "Kinetic positron energy: " << kineticEnergy/keV << " keV" << G4endl;
      G4cout << "Total available energy: " << totalAvailableEnergy/keV << " keV " << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Photon energy 1: " << photon1Energy/keV << " keV" << G4endl;
      G4cout << "Photon energy 2: " << photon2Energy/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (photon1Energy+photon2Energy)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {      
      G4double energyDiff = std::fabs(totalAvailableEnergy-photon1Energy-photon2Energy);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeAnnihilation: problem with energy conservation: " << 
	  (photon1Energy+photon2Energy)/keV << 
	  " keV (final) vs. " << 
	  totalAvailableEnergy/keV << " keV (initial)" << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeAnnihilationModel:: ComputeCrossSectionPerElectron(G4double energy)
{
  //
  // Penelope model to calculate cross section for positron annihilation.
  // The annihilation cross section per electron is calculated according 
  // to the Heitler formula
  //  W. Heitler, The quantum theory of radiation, Oxford University Press (1954)
  // in the assumptions of electrons free and at rest.
  //
  G4double gamma = 1.0+std::max(energy,1.0*eV)/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double f2 = gamma2-1.0;
  G4double f1 = std::sqrt(f2);
  G4double crossSection = fPielr2*((gamma2+4.0*gamma+1.0)*std::log(gamma+f1)/f2
			 - (gamma+3.0)/f1)/(gamma+1.0);
  return crossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeAnnihilationModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;  
  }
}
