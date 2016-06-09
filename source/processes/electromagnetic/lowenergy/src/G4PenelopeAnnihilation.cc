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
// Author: Luciano Pandola
//
// History:
// --------
// 02 Jul 2003   L.Pandola    First implementation
// 16 Mar 2004   L.Pandola    Removed unnecessary calls to std::pow(a,b)

#include "G4PenelopeAnnihilation.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"

// constructor
 
G4PenelopeAnnihilation::G4PenelopeAnnihilation(const G4String& processName)
  : G4VRestDiscreteProcess (processName),
    lowEnergyLimit(250*eV),
    highEnergyLimit(100*GeV),
    nBins(200),
    cutForLowEnergySecondaryPhotons(250.0*eV)
{
  meanFreePathTable = 0;
  if (verboseLevel > 0)
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: "
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV"
	     << G4endl;
    }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4PenelopeAnnihilation is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

// destructor
G4PenelopeAnnihilation::~G4PenelopeAnnihilation()
{
  if (meanFreePathTable) {
      meanFreePathTable->clearAndDestroy();
      delete meanFreePathTable;
   }
}
 

void G4PenelopeAnnihilation::BuildPhysicsTable(const G4ParticleDefinition& )
{
  G4double lowEdgeEnergy, value;
  G4PhysicsLogVector* dataVector;

  // Build mean free path table for the e+e- annihilation

  if (meanFreePathTable) {
    meanFreePathTable->clearAndDestroy(); delete meanFreePathTable;}

  meanFreePathTable = 
    new G4PhysicsTable(G4Material::GetNumberOfMaterials());
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4Material* material;

  for (size_t j=0;j<G4Material::GetNumberOfMaterials();j++)  
    { 
      //create physics vector then fill it ....
      dataVector = new G4PhysicsLogVector(lowEnergyLimit, 
					 highEnergyLimit,nBins);
      material = (*theMaterialTable)[j];
      const G4ElementVector* theElementVector = material->GetElementVector();
      const G4double* NbOfAtomsPerVolume = 
	material->GetVecNbOfAtomsPerVolume(); 
      for (G4int i=0;i<nBins;i++)      
        {
          lowEdgeEnergy = dataVector->GetLowEdgeEnergy(i);
	  G4double sigma=0.0 ;
	  for (size_t elm=0;elm<material->GetNumberOfElements();elm++)
	    {      
	      G4double Z = (*theElementVector)[elm]->GetZ();
	      sigma += NbOfAtomsPerVolume[elm] * Z * 
		calculateCrossSectionPerElectron(lowEdgeEnergy);
	    }       
	  if (sigma > DBL_MIN) 
	    {
	      value = 1.0/sigma;
	    }
	  else
	    {
	      value = DBL_MAX;
	    }
	  dataVector->PutValue(i,value);
        }
      meanFreePathTable->insertAt(j,dataVector);
    }
  PrintInfoDefinition();
}

G4double G4PenelopeAnnihilation::calculateCrossSectionPerElectron
                              (G4double ene)
{
  //Heitler dcs formula for annihilation with free electrons at rest
  G4double crossSection=0.0;
  G4double gamma = 1.0+std::max(ene,1.0*eV)/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double f2 = gamma2-1.0;
  G4double f1 = std::sqrt(f2);
  G4double pielr2 = pi*classic_electr_radius*classic_electr_radius;
  crossSection = pielr2*((gamma2+4.0*gamma+1.0)*std::log(gamma+f1)/f2
			 - (gamma+3.0)/f1)/(gamma+1.0);
  return crossSection;
}

G4VParticleChange* G4PenelopeAnnihilation::PostStepDoIt(const G4Track& aTrack,
                                                     const G4Step& )
{
   aParticleChange.Initialize(aTrack);
     
   const G4DynamicParticle* incidentPositron = aTrack.GetDynamicParticle();
   G4double kineticEnergy = incidentPositron->GetKineticEnergy();
   G4ParticleMomentum positronDirection = 
     incidentPositron->GetMomentumDirection();

   // Do not make anything if particle is stopped, the annihilation then
   // should be performed by the AtRestDoIt!
   if (aTrack.GetTrackStatus() == fStopButAlive) return &aParticleChange;

   //G4cout << "Sono nel PostStep" << G4endl;

   //Annihilation in flight
   G4double gamma = 1.0 + std::max(kineticEnergy,1.0*eV)/electron_mass_c2;
   G4double gamma21 = std::sqrt(gamma*gamma-1);
   G4double ani = 1.0+gamma;
   G4double chimin = 1.0/(ani+gamma21);
   G4double rchi = (1.0-chimin)/chimin;
   G4double gt0 = ani*ani-2.0;
   G4double epsilon=0.0, reject=0.0, test=0.0;
   do{
     epsilon = chimin*std::pow(rchi,G4UniformRand());
     reject = ani*ani*(1.0-epsilon)+2.0*gamma-(1.0/epsilon);
     test = G4UniformRand()*gt0-reject;
   }while(test>0);
   
   G4double totalAvailableEnergy = kineticEnergy + 2.0*electron_mass_c2;
   G4double photon1Energy = epsilon*totalAvailableEnergy;
   G4double photon2Energy = (1.0-epsilon)*totalAvailableEnergy;
   G4double cosTheta1 = (ani-1.0/epsilon)/gamma21;
   G4double cosTheta2 = (ani-1.0/(1.0-epsilon))/gamma21;

   aParticleChange.SetNumberOfSecondaries(2);
   G4double localEnergyDeposit = 0.; 

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

   if (photon1Energy > cutForLowEnergySecondaryPhotons) {
     G4ThreeVector photon1Direction (dirx1,diry1,dirz1);
     photon1Direction.rotateUz(positronDirection);   
     // create G4DynamicParticle object for the particle1  
     G4DynamicParticle* aParticle1= new G4DynamicParticle (G4Gamma::Gamma(),
							   photon1Direction, 
							   photon1Energy);
     aParticleChange.AddSecondary(aParticle1);
   }
   else  localEnergyDeposit += photon1Energy;  

   if (photon2Energy > cutForLowEnergySecondaryPhotons) {
     G4ThreeVector photon2Direction(dirx2,diry2,dirz2);
     photon2Direction.rotateUz(positronDirection); 
     // create G4DynamicParticle object for the particle2 
     G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
							   photon2Direction,
							   photon2Energy);
     aParticleChange.AddSecondary(aParticle2);
   }   
   else  localEnergyDeposit += photon2Energy;
     
   aParticleChange.ProposeLocalEnergyDeposit(localEnergyDeposit);

   aParticleChange.ProposeMomentumDirection( 0., 0., 0. );
   aParticleChange.ProposeEnergy(0.); 
   aParticleChange.ProposeTrackStatus(fStopAndKill);

   return &aParticleChange;
}

 
G4VParticleChange* G4PenelopeAnnihilation::AtRestDoIt(const G4Track& aTrack,
                                                  const G4Step& )
{
   aParticleChange.Initialize(aTrack);
   aParticleChange.SetNumberOfSecondaries(2); 
   G4double cosTheta = -1.0+2.0*G4UniformRand();
   G4double sinTheta = std::sqrt(1.0-cosTheta*cosTheta);
   G4double phi = twopi*G4UniformRand();
   //G4cout << "cosTheta: " << cosTheta << " sinTheta: " << sinTheta << G4endl;
   //G4cout << "phi: " << phi << G4endl;
   G4ThreeVector direction (sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);   
   aParticleChange.AddSecondary(new G4DynamicParticle (G4Gamma::Gamma(),
                                            direction, electron_mass_c2) );
   aParticleChange.AddSecondary(new G4DynamicParticle (G4Gamma::Gamma(),
                                           -direction, electron_mass_c2) ); 

   aParticleChange.ProposeLocalEnergyDeposit(0.);

   // Kill the incident positron 
   //
   aParticleChange.ProposeTrackStatus(fStopAndKill);
      
   return &aParticleChange;
}

void G4PenelopeAnnihilation::PrintInfoDefinition()
{
  G4String comments = "Total cross section from Heilter formula" 
                      "(annihilation into 2 photons).";
  comments += "\n      Gamma energies sampled according Heitler"; 
  comments += "\n      It can be used for positrons";
  comments += " in the energy range [250eV,100GeV].";
  G4cout << G4endl << GetProcessName() << ":  " << comments
	 << G4endl;
}         

G4bool G4PenelopeAnnihilation::IsApplicable(
					 const G4ParticleDefinition& particle)
{
   return ( &particle == G4Positron::Positron() ); 
}

 
G4double G4PenelopeAnnihilation::GetMeanFreePath(const G4Track& aTrack,
                                                     G4double,
                                                     G4ForceCondition*)
{
  const G4DynamicParticle* incidentPositron = aTrack.GetDynamicParticle();
  G4double kineticEnergy = incidentPositron->GetKineticEnergy();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();
  G4int materialIndex = material->GetIndex();
  if (kineticEnergy<lowEnergyLimit)
    {
      kineticEnergy = lowEnergyLimit;
    }

  G4double meanFreePath;
  G4bool isOutRange ;

  if (kineticEnergy>highEnergyLimit) 
    {
      meanFreePath = DBL_MAX;
    }
  else 
    {
      meanFreePath = (*meanFreePathTable)(materialIndex)->
                    GetValue(kineticEnergy,isOutRange);
    }

  return meanFreePath; 
} 

G4double G4PenelopeAnnihilation::GetMeanLifeTime(const G4Track&,
                                                     G4ForceCondition*)
 
{
  return 0.0; 
} 
