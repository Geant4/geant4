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
//	G4IonCoulombScatteringModel.cc
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4IonCoulombScatteringModel
//
// Author:      Cristina Consolandi
//
// Creation date: 05.10.2010 from G4eCoulombScatteringModel 
//                               & G4CoulombScatteringModel
//
// Class Description:
//      Single Scattering Model for
//      for protons, alpha and heavy Ions
//
// Reference:
//      M.J. Boschini et al. "Nuclear and Non-Ionizing Energy-Loss 
//      for Coulomb ScatteredParticles from Low Energy up to Relativistic 
//      Regime in Space Radiation Environment"
//      Accepted for publication in the Proceedings of  the  ICATPP Conference
//      on Cosmic Rays for Particle and Astroparticle Physics, Villa  Olmo, 7-8
//      October,  2010, to be published by World Scientific (Singapore).
//
//      Available for downloading at:
//	http://arxiv.org/abs/1011.4822
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4IonCoulombScatteringModel.hh"
#include "Randomize.hh"
//#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4IonCoulombScatteringModel::G4IonCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),

    cosThetaMin(1.0),
    isInitialised(false)
{
  	fNistManager = G4NistManager::Instance();
  	theParticleTable = G4ParticleTable::GetParticleTable();
  	theProton   = G4Proton::Proton();

	pCuts=0;
  	currentMaterial = 0;
  	currentElement  = 0;
        currentCouple = 0;

        lowEnergyLimit  = 100*eV;
  	recoilThreshold = 0.*eV;
	heavycorr =0;
  	particle = 0;
	mass=0;
	currentMaterialIndex = -1;

  	ioncross = new G4IonCoulombCrossSection(); 

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4IonCoulombScatteringModel::~G4IonCoulombScatteringModel()
{ delete  ioncross;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
                                           const G4DataVector&  )
{
  	SetupParticle(p);
  	currentCouple = 0;
	currentMaterialIndex = -1;
  	cosThetaMin = cos(PolarAngleLimit());
  	ioncross->Initialise(p,cosThetaMin);
 
        pCuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);


  	if(!isInitialised) {
	  isInitialised = true;
	  fParticleChange = GetParticleChangeForGamma();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4IonCoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double, 
				G4double cutEnergy,
				G4double)
{

  	SetupParticle(p);
 
  	G4double xsec =0.0;
	if(kinEnergy < lowEnergyLimit) return xsec;

  	DefineMaterial(CurrentCouple());

        G4int iz          = G4int(Z);

	//from lab to pCM & mu_rel of effective particle
        ioncross->SetupKinematic(kinEnergy, cutEnergy,iz);


        ioncross->SetupTarget(Z, kinEnergy, heavycorr);

  	xsec = ioncross->NuclearCrossSection();

//cout<< "..........xsec "<<G4BestUnit(xsec,"Surface") <<endl;
  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4IonCoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cutEnergy, 
			       G4double)
{
  	G4double kinEnergy = dp->GetKineticEnergy();

  	if(kinEnergy < lowEnergyLimit) return;
	
  	DefineMaterial(couple);

  	SetupParticle(dp->GetDefinition());

	// Choose nucleus
  	currentElement = SelectRandomAtom(couple,particle,
                                    kinEnergy,cutEnergy,kinEnergy);

	G4double Z  = currentElement->GetZ();
  	G4int iz          = G4int(Z);
  	G4int ia = SelectIsotopeNumber(currentElement);
  	G4double m2 = G4NucleiProperties::GetNuclearMass(ia, iz);



	G4double xsec= ComputeCrossSectionPerAtom(particle,kinEnergy, Z,
                                kinEnergy, cutEnergy, kinEnergy) ;
	if(xsec == 0.0)return;
    
	//scattering angle, z1 == (1-cost)
  	G4double z1 = ioncross->SampleCosineTheta(); 

  	if(z1 <= 0.0) { return; }

  	G4double cost = 1.0 - z1;
  	G4double sint = sqrt(z1*(1.0 + cost));
  	G4double phi  = twopi * G4UniformRand();


  	// kinematics in the Lab system
  	G4double etot = kinEnergy + mass;
  	G4double mom2=  kinEnergy*(kinEnergy+2.0*mass);
  	G4double ptot = sqrt(mom2);

	//CM particle 1
  	G4double bet  = ptot/(etot + m2);
  	G4double gam  = 1.0/sqrt((1.0 - bet)*(1.0 + bet));

	//CM 	
  	G4double momCM2= ioncross->GetMomentum2();
  	G4double momCM =std::sqrt(momCM2); 
        //energy & momentum after scattering of incident particle
	G4double pxCM = momCM*sint*cos(phi);
        G4double pyCM = momCM*sint*sin(phi);
        G4double pzCM = momCM*cost;
  	G4double eCM  = sqrt(momCM2 + mass*mass);

	//CM--->Lab
  	G4ThreeVector v1(pxCM , pyCM, gam*(pzCM + bet*eCM));
  	G4ThreeVector dir = dp->GetMomentumDirection(); 

  	G4ThreeVector newDirection = v1.unit();
  	newDirection.rotateUz(dir);   
  
	fParticleChange->ProposeMomentumDirection(newDirection);   
  
	// recoil.......................................
	G4double trec =(1.0 - cost)* m2*(etot*etot - mass*mass )/
				(mass*mass + m2*m2+ 2.*m2*etot);

        G4double finalT = kinEnergy - trec;


  	if(finalT <= lowEnergyLimit) { 
    		trec = kinEnergy;  
    		finalT = 0.0;
  		} 
    
  	fParticleChange->SetProposedKineticEnergy(finalT);

  	G4double tcut = recoilThreshold;
        if(pCuts) { tcut= std::max(tcut,(*pCuts)[currentMaterialIndex]); 

			//G4cout<<" tcut eV "<<tcut/eV<<endl;
			}
 
  	if(trec > tcut) {
    		G4ParticleDefinition* ion = theParticleTable->FindIon(iz, ia, 0, iz);
    		G4double plab = sqrt(finalT*(finalT + 2.0*mass));
    		G4ThreeVector p2 = (ptot*dir - plab*newDirection).unit();
    		G4DynamicParticle* newdp  = new G4DynamicParticle(ion, p2, trec);
    			fvect->push_back(newdp);
  	} else if(trec > 0.0) {
    		fParticleChange->ProposeLocalEnergyDeposit(trec);
    		fParticleChange->ProposeNonIonizingEnergyDeposit(trec);
  	}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
		
