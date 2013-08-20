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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class
// File name:     G4PAIModel.cc
//
// Author: Vladimir.Grichine@cern.ch on base of V.Ivanchenko model interface
//
// Creation date: 05.10.2003
//
// Modifications:
//
// 17.08.04 V.Grichine, bug fixed for Tkin<=0 in SampleSecondary
// 16.08.04 V.Grichine, bug fixed in massRatio for DEDX, CrossSection, 
//          SampleSecondary
// 08.04.05 Major optimisation of internal interfaces (V.Ivantchenko)
// 26.07.09 Fixed logic to work with several materials (V.Ivantchenko)
// 21.11.10 V. Grichine verbose flag for protons and G4PAYySection to 
//          check sandia table 
// 12.06.13 V. Grichine Bug fixed in SampleSecondaries for scaled Tkin 
//          (fMass -> proton_mass_c2)
// 19.08.13 V.Ivanchenko extract data handling to G4PAIModelData class 
//          added sharing of internal data between threads (MT migration)
//

#include "G4PAIModel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Region.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"
#include "G4OrderedTable.hh"

#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4PAIModelData.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIModel::G4PAIModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
    fVerbose(0),
    fModelData(0),
    fParticle(0)
{  
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  fParticleChange = 0;

  if(p) { SetParticle(p); }
  else  { SetParticle(fElectron); }

  fLowKinEnergy = 0.1*MeV;

  isInitialised = false;
}

////////////////////////////////////////////////////////////////////////////

G4PAIModel::~G4PAIModel()
{
  if(IsMaster()) { delete fModelData; }
}

////////////////////////////////////////////////////////////////////////////

void G4PAIModel::Initialise(const G4ParticleDefinition* p,
			    const G4DataVector& cuts)
{
  if(fVerbose > 0) {
    G4cout<<"G4PAIModel::Initialise for "<<p->GetParticleName()<<G4endl;
  }

  if(isInitialised) { return; }
  isInitialised = true;

  SetParticle(p);
  fParticleChange = GetParticleChangeForLoss();

  if(IsMaster()) { 
    if(!fModelData) {

      G4double tmin = fRatio*std::max(fLowKinEnergy, LowEnergyLimit());
      G4double tmax =  HighEnergyLimit()*fRatio;
      fModelData = new G4PAIModelData(tmin, tmax, fVerbose);
    }
    // Prepare initialization
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    size_t numOfMat   = G4Material::GetNumberOfMaterials();
    size_t numRegions = fPAIRegionVector.size();

    for(size_t iReg = 0; iReg < numRegions; ++iReg) {
      const G4Region* curReg = fPAIRegionVector[iReg];
      G4Region* reg = const_cast<G4Region*>(curReg);

      for(size_t jMat = 0; jMat < numOfMat; ++jMat) {
	G4Material* mat = (*theMaterialTable)[jMat];
	const G4MaterialCutsCouple* cutCouple = reg->FindCouple(mat);
	//G4cout << "Couple <" << fCutCouple << G4endl;
	if(cutCouple) {
	  /*
	    G4cout << "Reg <" <<curReg->GetName() << ">  mat <" 
	    << fMaterial->GetName() << ">  fCouple= " 
	    << fCutCouple << " idx= " << fCutCouple->GetIndex()
	    <<"  " << p->GetParticleName() <<G4endl;
	    // G4cout << cuts.size() << G4endl;
	    */
	  // check if this couple is not already initialized
	  size_t n = fMaterialCutsCoupleVector.size();
	  if(0 < n) {
	    for(size_t i=0; i<fMaterialCutsCoupleVector.size(); ++i) {
	      if(cutCouple == fMaterialCutsCoupleVector[i]) {
		break;
	      }
	    }
	  }
	  // initialise data banks
	  fMaterialCutsCoupleVector.push_back(cutCouple);
	  G4double deltaCutInKinEnergy = cuts[cutCouple->GetIndex()];
	  fModelData->Initialise(cutCouple, deltaCutInKinEnergy, this);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::ComputeDEDXPerVolume(const G4Material*,
					  const G4ParticleDefinition* p,
					  G4double kineticEnergy,
					  G4double cutEnergy)
{
  G4int coupleIndex = FindCoupleIndex(CurrentCouple());
  if(0 > coupleIndex) { return 0.0; }

  G4double cut = std::min(MaxSecondaryEnergy(p, kineticEnergy), cutEnergy);

  G4double scaledTkin = kineticEnergy*fRatio;
 
  return fChargeSquare*fModelData->DEDXPerVolume(coupleIndex, scaledTkin, cut);
}

/////////////////////////////////////////////////////////////////////////

void G4PAIModel::InitialiseLocal(const G4ParticleDefinition*, 
				 G4VEmModel* masterModel)
{
  fModelData = static_cast<G4PAIModel*>(masterModel)->GetPAIModelData();
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::CrossSectionPerVolume( const G4Material*,
					    const G4ParticleDefinition* p,
					    G4double kineticEnergy,
					    G4double cutEnergy,
					    G4double maxEnergy  ) 
{
  G4int coupleIndex = FindCoupleIndex(CurrentCouple());
  if(0 > coupleIndex) { return 0.0; }

  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(tmax <= cutEnergy) { return 0.0; }

  G4double scaledTkin = kineticEnergy*fRatio;

  return fChargeSquare*fModelData->CrossSectionPerVolume(coupleIndex, 
							 scaledTkin, 
							 cutEnergy, tmax);
}

///////////////////////////////////////////////////////////////////////////
//
// It is analog of PostStepDoIt in terms of secondary electron.
//
 
void G4PAIModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
				   const G4MaterialCutsCouple* matCC,
				   const G4DynamicParticle* dp,
				   G4double tmin,
				   G4double maxEnergy)
{
  G4int coupleIndex = FindCoupleIndex(matCC);
  if(0 > coupleIndex) { return; }

  SetParticle(dp->GetDefinition());
  G4double kineticEnergy = dp->GetKineticEnergy();

  G4double tmax = MaxSecondaryEnergy(fParticle, kineticEnergy);
  if(maxEnergy < tmax) { tmax = maxEnergy; }
  if(tmin >= tmax) { return; }

  G4ThreeVector direction= dp->GetMomentumDirection();
  G4double scaledTkin    = kineticEnergy*fRatio;
  G4double totalEnergy   = kineticEnergy + fMass;
  G4double pSquare       = kineticEnergy*(totalEnergy+fMass);
 
  G4double deltaTkin = 
    fModelData->SamplePostStepTransfer(coupleIndex, scaledTkin);

  // G4cout<<"G4PAIModel::SampleSecondaries; deltaKIn = "<<deltaTkin/keV
  // <<" keV "<<G4endl;

  if( deltaTkin <= 0. && fVerbose > 0) 
  {
    G4cout<<"G4PAIModel::SampleSecondary e- deltaTkin = "<<deltaTkin<<G4endl;
  }
  if( deltaTkin <= 0.) { return; }

  if( deltaTkin > tmax) { deltaTkin = tmax; }

  G4double deltaTotalMomentum = 
    sqrt(deltaTkin*(deltaTkin + 2.*electron_mass_c2 ));
  G4double totalMomentum      = sqrt(pSquare);
  G4double costheta           = deltaTkin*(totalEnergy + electron_mass_c2)
                                /(deltaTotalMomentum * totalMomentum);

  if( costheta > 1.0 ) { costheta = 1.0; }
  G4double sintheta = 0.0;
  G4double sin2 = 1. - costheta*costheta;
  if( sin2 > 0.) { sintheta = sqrt(sin2); }

  //  direction of the delta electron
  G4double phi  = twopi*G4UniformRand(); 

  G4ThreeVector deltaDirection(sintheta*cos(phi),sintheta*sin(phi),costheta);
  deltaDirection.rotateUz(direction);

  // primary change
  kineticEnergy -= deltaTkin;
  G4ThreeVector dir = 
    totalMomentum*direction - deltaTotalMomentum*deltaDirection;
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  // create G4DynamicParticle object for e- delta ray 
  G4DynamicParticle* deltaRay = new G4DynamicParticle();
  deltaRay->SetDefinition(fElectron);
  deltaRay->SetKineticEnergy( deltaTkin ); //!!! trick for last steps/2.0???
  deltaRay->SetMomentumDirection(deltaDirection); 

  vdp->push_back(deltaRay);
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations( const G4MaterialCutsCouple* matCC,
                                         const G4DynamicParticle* aParticle,
					 G4double, G4double step,
					 G4double eloss)
{
  G4int coupleIndex = FindCoupleIndex(matCC);
  if(0 > coupleIndex) { return eloss; }

  SetParticle(aParticle->GetDefinition());

  /*
  G4cout << "G4PAIModel::SampleFluctuations step(mm)= "<< step/mm
	 << "  Eloss(keV)= " << eloss/keV  << " in " 
	 << matCC->Getmaterial()->GetName() << G4endl;
  */

  G4double Tkin       = aParticle->GetKineticEnergy();
  G4double scaledTkin = Tkin*fRatio;

  G4double loss = fModelData->SampleAlongStepTransfer(coupleIndex, Tkin,
						      scaledTkin,
						      step*fChargeSquare);
  
  // G4cout<<"PAIModel AlongStepLoss = "<<loss/keV<<" keV, on step = "
  //<<step/mm<<" mm"<<G4endl; 
  return loss;

}

//////////////////////////////////////////////////////////////////////
//
// Returns the statistical estimation of the energy loss distribution variance
//


G4double G4PAIModel::Dispersion( const G4Material* material, 
                                 const G4DynamicParticle* aParticle,
 				       G4double tmax, 
			               G4double step       )
{
  G4double particleMass  = aParticle->GetMass();
  G4double electronDensity = material->GetElectronDensity();
  G4double kineticEnergy = aParticle->GetKineticEnergy();
  G4double q = aParticle->GetCharge()/eplus;
  G4double etot = kineticEnergy + particleMass;
  G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*particleMass)/(etot*etot);
  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * step
                 * electronDensity * q * q;

  return siga;
  /*
  G4double loss, sumLoss=0., sumLoss2=0., sigma2, meanLoss=0.;
  for(G4int i = 0; i < fMeanNumber; i++)
  {
    loss      = SampleFluctuations(material,aParticle,tmax,step,meanLoss);
    sumLoss  += loss;
    sumLoss2 += loss*loss;
  }
  meanLoss = sumLoss/fMeanNumber;
  sigma2   = meanLoss*meanLoss + (sumLoss2-2*sumLoss*meanLoss)/fMeanNumber;
  return sigma2;
  */
}

/////////////////////////////////////////////////////////////////////

G4double G4PAIModel::MaxSecondaryEnergy( const G4ParticleDefinition* p,
					 G4double kinEnergy) 
{
  SetParticle(p);
  G4double tmax = kinEnergy;
  if(p == fElectron) { tmax *= 0.5; }
  else if(p != fPositron) { 
    G4double ratio= electron_mass_c2/fMass;
    G4double gamma= kinEnergy/fMass + 1.0;
    tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  }
  return tmax;
}

///////////////////////////////////////////////////////////////

void G4PAIModel::DefineForRegion(const G4Region* r) 
{
  fPAIRegionVector.push_back(r);
}

//
//
/////////////////////////////////////////////////
