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
//	G4eSingleCoulombScatteringModel.cc
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4eSingleCoulombScatteringModel
//
// Author:      Cristina Consolandi
//
// Creation date: 20.10.2012  
//                           
//	Class Description:
//	Single Scattering model for electron-nuclei interaction.
//	Suitable for high energy electrons and low scattering angles.
//
//
//	 Reference:
//	M.J. Boschini et al.
//	"Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//	Proc. of the 13th International Conference on Particle Physics and Advanced Technology 
//	(13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//	Available at: http://arxiv.org/abs/1111.4042v4
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4eSingleCoulombScatteringModel.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NucleiProperties.hh"

#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eSingleCoulombScatteringModel::G4eSingleCoulombScatteringModel(const G4String& nam)
  : G4VEmModel(nam),

    cosThetaMin(1.0),
    isInitialised(false)
{
  	fNistManager = G4NistManager::Instance();
  	theParticleTable = G4ParticleTable::GetParticleTable();
	fParticleChange = 0;

	pCuts=0;
  	currentMaterial = 0;
  	currentElement  = 0;
        currentCouple = 0;

        lowEnergyLimit  = 0*eV;
  	recoilThreshold = 0.*eV;
  	particle = 0;
	mass=0;
	currentMaterialIndex = -1;

  	Mottcross = new G4ScreeningMottCrossSection(); 

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eSingleCoulombScatteringModel::~G4eSingleCoulombScatteringModel()
{ delete  Mottcross;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eSingleCoulombScatteringModel::Initialise(const G4ParticleDefinition* p,
                                           const G4DataVector&  )
{
  	SetupParticle(p);
  	currentCouple = 0;
	currentMaterialIndex = -1;
  	cosThetaMin = cos(PolarAngleLimit());
  	Mottcross->Initialise(p,cosThetaMin);
 
        pCuts = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(3);


  	if(!isInitialised) {
	  isInitialised = true;
	  fParticleChange = GetParticleChangeForGamma();
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eSingleCoulombScatteringModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition* p,
				G4double kinEnergy, 
				G4double Z, 
				G4double , 
				G4double, 
				G4double )
{

  	SetupParticle(p);
 
  	G4double cross =0.0;
	if(kinEnergy < lowEnergyLimit) return cross;

  	DefineMaterial(CurrentCouple());

	//Total Cross section
        Mottcross->SetupKinematic(kinEnergy, Z);
  	cross = Mottcross->NuclearCrossSection();

	//cout<< "....cross "<<G4BestUnit(cross,"Surface") << " cm2 "<< cross/cm2 <<endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eSingleCoulombScatteringModel::SampleSecondaries(
			       std::vector<G4DynamicParticle*>* fvect,
			       const G4MaterialCutsCouple* couple,
			       const G4DynamicParticle* dp,
			       G4double cutEnergy, 
			       G4double)
{
  	G4double kinEnergy = dp->GetKineticEnergy();
	//cout<<"--- kinEnergy "<<kinEnergy<<endl;


  	if(kinEnergy < lowEnergyLimit) return;
	
  	DefineMaterial(couple);
  	SetupParticle(dp->GetDefinition());

	// Choose nucleus
  	currentElement = SelectRandomAtom(couple,particle,
                                    kinEnergy,cutEnergy,kinEnergy);//last two :cutEnergy= min e kinEnergy=max

	G4double Z  = currentElement->GetZ();
  	G4int iz    = G4int(Z);
  	G4int ia = SelectIsotopeNumber(currentElement);

	//cout<<"Element "<<currentElement->GetName()<<endl;;	

	G4double cross= Mottcross->GetTotalCross();


	if(cross == 0.0)return;
    		
  	G4ThreeVector dir = dp->GetMomentumDirection(); //old direction
	G4ThreeVector newDirection=Mottcross->GetNewDirection();//new direction
  	newDirection.rotateUz(dir);   
  
	fParticleChange->ProposeMomentumDirection(newDirection);   
  
	//Recoil energy
	G4double trec= Mottcross->GetTrec();
	//Energy after scattering	
        G4double finalT = kinEnergy - trec;


  	if(finalT <= lowEnergyLimit) { 
    		trec = kinEnergy;  
    		finalT = 0.0;
  		} 
    
  	fParticleChange->SetProposedKineticEnergy(finalT);

  	G4double tcut = recoilThreshold;
        if(pCuts) { tcut= std::min(tcut,(*pCuts)[currentMaterialIndex]); 
		}
 

  	if(trec > tcut) {

		//cout<<"Trec "<<trec/eV<<endl;
    		G4ParticleDefinition* ion = theParticleTable->GetIon(iz, ia, 0.0);

		//incident before scattering
		G4double ptot=sqrt(Mottcross->GetMom2Lab());
		//incident after scattering
    		G4double plab = sqrt(finalT*(finalT + 2.0*mass));
    		G4ThreeVector p2 = (ptot*dir - plab*newDirection).unit();
		//secondary particle
    		G4DynamicParticle* newdp  = new G4DynamicParticle(ion, p2, trec);
    			fvect->push_back(newdp);
  			}

	else if(trec > 0.0) {
                fParticleChange->ProposeNonIonizingEnergyDeposit(trec);
   		if(trec< tcut) fParticleChange->ProposeLocalEnergyDeposit(trec);
  		}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
		
