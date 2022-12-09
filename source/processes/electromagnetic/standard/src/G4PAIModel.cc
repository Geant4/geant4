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
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialTable.hh"
#include "G4RegionStore.hh"

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
#include "G4DeltaAngle.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIModel::G4PAIModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),G4VEmFluctuationModel(nam),
    fVerbose(0),
    fModelData(nullptr),
    fParticle(nullptr)
{  
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  fParticleChange = nullptr;

  if(p) { SetParticle(p); }
  else  { SetParticle(fElectron); }

  // default generator
  SetAngularDistribution(new G4DeltaAngle());
  fLowestTcut = 12.5*CLHEP::eV;
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
  if(fVerbose > 1) {
    G4cout<<"G4PAIModel::Initialise for "<<p->GetParticleName()<<G4endl;
  }
  SetParticle(p);
  fParticleChange = GetParticleChangeForLoss();

  if(IsMaster()) { 

    delete fModelData;      
    fMaterialCutsCoupleVector.clear(); 
    if(fVerbose > 1) {
      G4cout << "G4PAIModel instantiates data for  " << p->GetParticleName()
	     << G4endl;
    }  
    G4double tmin = LowEnergyLimit()*fRatio;
    G4double tmax = HighEnergyLimit()*fRatio;
    fModelData = new G4PAIModelData(tmin, tmax, fVerbose);
  
    // Prepare initialization
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    std::size_t numOfMat   = G4Material::GetNumberOfMaterials();
    std::size_t numRegions = fPAIRegionVector.size();

    // protect for unit tests
    if(0 == numRegions) {
      G4Exception("G4PAIModel::Initialise()","em0106",JustWarning,
                  "no G4Regions are registered for the PAI model - World is used");
      fPAIRegionVector.push_back(G4RegionStore::GetInstance()
				 ->GetRegion("DefaultRegionForTheWorld", false));
      numRegions = 1;
    }

    if(fVerbose > 1) {
      G4cout << "G4PAIModel is defined for " << numRegions << " regions "   
	     << "; number of materials " << numOfMat << G4endl;
    } 
    for(std::size_t iReg = 0; iReg<numRegions; ++iReg) {
      const G4Region* curReg = fPAIRegionVector[iReg];
      G4Region* reg = const_cast<G4Region*>(curReg);

      for(std::size_t jMat = 0; jMat<numOfMat; ++jMat) {
	G4Material* mat = (*theMaterialTable)[jMat];
	const G4MaterialCutsCouple* cutCouple = reg->FindCouple(mat);
	std::size_t n = fMaterialCutsCoupleVector.size();
	/*
	G4cout << "Region: " << reg->GetName() << "  " << reg
	       << " Couple " << cutCouple 
	       << " PAI defined for " << n << " couples"
	       << " jMat= " << jMat << "  " << mat->GetName()
	       << G4endl;
	*/
	if(nullptr != cutCouple) {
	  if(fVerbose > 1) {
	    G4cout << "Region <" << curReg->GetName() << ">  mat <" 
		   << mat->GetName() << ">  CoupleIndex= " 
		   << cutCouple->GetIndex()
		   << "  " << p->GetParticleName()
		   << " cutsize= " << cuts.size() << G4endl;
	  }
	  // check if this couple is not already initialized
	  G4bool isnew = true;
	  if(0 < n) {
	    for(std::size_t i=0; i<n; ++i) {
	      G4cout << i << G4endl;
	      if(cutCouple == fMaterialCutsCoupleVector[i]) {
		isnew = false;
		break;
	      }
	    }
	  }
	  // initialise data banks
	  // G4cout << "   isNew: " << isnew << "  " << cutCouple << G4endl;
	  if(isnew) { 
	    fMaterialCutsCoupleVector.push_back(cutCouple); 
	    fModelData->Initialise(cutCouple, this);
	  }
	}
      }
    }
    InitialiseElementSelectors(p, cuts);
  }
}

/////////////////////////////////////////////////////////////////////////

void G4PAIModel::InitialiseLocal(const G4ParticleDefinition* p, 
				 G4VEmModel* masterModel)
{
  SetParticle(p);
  fModelData = static_cast<G4PAIModel*>(masterModel)->GetPAIModelData();
  fMaterialCutsCoupleVector = 
    static_cast<G4PAIModel*>(masterModel)->GetVectorOfCouples();
  SetElementSelectors(masterModel->GetElementSelectors());
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::MinEnergyCut(const G4ParticleDefinition*,
				  const G4MaterialCutsCouple*)
{
  return fLowestTcut;
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::ComputeDEDXPerVolume(const G4Material*,
					  const G4ParticleDefinition* p,
					  G4double kineticEnergy,
					  G4double cutEnergy)
{  
  //G4cout << "===1=== " << CurrentCouple() 
  //	 << "  idx= " << CurrentCouple()->GetIndex()
  //	 << "   " <<  fMaterialCutsCoupleVector[0]
  //	 << G4endl; 
  G4int coupleIndex = FindCoupleIndex(CurrentCouple());
  //G4cout << "===2=== " << coupleIndex << G4endl;
  if(0 > coupleIndex) { return 0.0; }

  G4double cut = std::min(MaxSecondaryEnergy(p, kineticEnergy), cutEnergy);

  G4double scaledTkin = kineticEnergy*fRatio;
 
  return fChargeSquare*fModelData->DEDXPerVolume(coupleIndex, scaledTkin, 
						 cut);
}

/////////////////////////////////////////////////////////////////////////

G4double G4PAIModel::CrossSectionPerVolume( const G4Material*,
					    const G4ParticleDefinition* p,
					    G4double kineticEnergy,
					    G4double cutEnergy,
					    G4double maxEnergy  ) 
{
  //G4cout << "===3=== " << CurrentCouple() 
  //	 << "  idx= " << CurrentCouple()->GetIndex()
  //	 << "   " <<  fMaterialCutsCoupleVector[0]
  //	 << G4endl; 
  G4int coupleIndex = FindCoupleIndex(CurrentCouple());
  //G4cout << "===4=== " << coupleIndex << G4endl;
  if(0 > coupleIndex) { return 0.0; }

  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);
  if(tmax <= cutEnergy) { return 0.0; }

  G4double scaledTkin = kineticEnergy*fRatio;

  return fChargeSquare*fModelData->CrossSectionPerVolume(coupleIndex, 
							 scaledTkin, 
							 cutEnergy, 
							 tmax);
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

  //G4cout << "G4PAIModel::SampleSecondaries: coupleIndex= "<<coupleIndex<<G4endl;
  if(0 > coupleIndex) { return; }

  SetParticle(dp->GetDefinition());
  G4double kineticEnergy = dp->GetKineticEnergy();

  G4double tmax = MaxSecondaryEnergy(fParticle, kineticEnergy);
  if(maxEnergy < tmax) { tmax = maxEnergy; }
  if(tmin >= tmax) { return; }

  G4ThreeVector direction= dp->GetMomentumDirection();
  G4double scaledTkin    = kineticEnergy*fRatio;
  G4double totalEnergy   = kineticEnergy + fMass;
  G4double totalMomentum = sqrt(kineticEnergy*(totalEnergy+fMass));

  G4double deltaTkin = 
    fModelData->SamplePostStepTransfer(coupleIndex, scaledTkin, tmin, tmax);

  //G4cout<<"G4PAIModel::SampleSecondaries; deltaKIn = "<<deltaTkin/keV
  //	<<" keV "<< " Escaled(MeV)= " << scaledTkin << G4endl;

  if( !(deltaTkin <= 0.) && !(deltaTkin > 0)) {
    G4cout<<"G4PAIModel::SampleSecondaries; deltaKIn = "<<deltaTkin/keV
	  <<" keV "<< " Escaled(MeV)= " << scaledTkin << G4endl;
    return;
  }
  if( deltaTkin <= 0.) { return; }

  if( deltaTkin > tmax) { deltaTkin = tmax; }

  const G4Element* anElement = SelectTargetAtom(matCC, fParticle, kineticEnergy,
                                                dp->GetLogKineticEnergy());

  G4int Z = G4lrint(anElement->GetZ());
 
  auto deltaRay = new G4DynamicParticle(fElectron,
	    GetAngularDistribution()->SampleDirection(dp, deltaTkin,
		 				      Z, matCC->GetMaterial()),
						      deltaTkin);

  // primary change
  kineticEnergy -= deltaTkin;
  G4ThreeVector dir = totalMomentum*direction - deltaRay->GetMomentum();
  direction = dir.unit();
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(direction);

  vdp->push_back(deltaRay);
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIModel::SampleFluctuations(const G4MaterialCutsCouple* matCC,
                                        const G4DynamicParticle* aParticle,
                                        const G4double tcut,
					const G4double, 
                                        const G4double step,
					const G4double eloss)
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
						      scaledTkin, tcut,
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
 				 const G4double tcut,
 				 const G4double tmax,
			         const G4double step )
{
  G4double particleMass  = aParticle->GetMass();
  G4double electronDensity = material->GetElectronDensity();
  G4double kineticEnergy = aParticle->GetKineticEnergy();
  G4double q = aParticle->GetCharge()/eplus;
  G4double etot = kineticEnergy + particleMass;
  G4double beta2 = kineticEnergy*(kineticEnergy + 2.0*particleMass)/(etot*etot);
  G4double siga  = (tmax/beta2 - 0.5*tcut) * twopi_mc2_rcl2 * step
                 * electronDensity * q * q;

  return siga;
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

///////////////////////////////////////////////////////////////

