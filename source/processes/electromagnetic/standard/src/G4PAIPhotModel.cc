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
// File name:     G4PAIPhotModel.cc
//
// Author: Vladimir.Grichine@cern.ch on base of G4PAIModel MT interface
//
// Creation date: 07.10.2013
//
// Modifications:
//
//

#include "G4PAIPhotModel.hh"

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
#include "G4RegionStore.hh"

#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4Poisson.hh"
#include "G4Step.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4PAIPhotData.hh"
#include "G4DeltaAngle.hh"

////////////////////////////////////////////////////////////////////////

using namespace std;

G4PAIPhotModel::G4PAIPhotModel(const G4ParticleDefinition* p, const G4String& nam)
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

G4PAIPhotModel::~G4PAIPhotModel()
{
  //G4cout << "G4PAIPhotModel::~G4PAIPhotModel() " << this << G4endl;
  if(IsMaster()) { delete fModelData; fModelData = nullptr; }
}

////////////////////////////////////////////////////////////////////////////

void G4PAIPhotModel::Initialise(const G4ParticleDefinition* p,
 			        const G4DataVector& cuts)
{
  if(fVerbose > 0) 
  {
    G4cout<<"G4PAIPhotModel::Initialise for "<<p->GetParticleName()<<G4endl;
  }
  SetParticle(p);
  fParticleChange = GetParticleChangeForLoss();

  if( IsMaster() ) 
  { 

    InitialiseElementSelectors(p, cuts);

    delete fModelData;
    fMaterialCutsCoupleVector.clear(); 

    G4double tmin = LowEnergyLimit()*fRatio;
    G4double tmax = HighEnergyLimit()*fRatio;
    fModelData = new G4PAIPhotData(tmin, tmax, fVerbose);
    
    // Prepare initialization
    const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    size_t numOfMat   = G4Material::GetNumberOfMaterials();
    size_t numRegions = fPAIRegionVector.size();

    // protect for unit tests
    if(0 == numRegions) {
      G4Exception("G4PAIModel::Initialise()","em0106",JustWarning,
                  "no G4Regions are registered for the PAI model - World is used");
      fPAIRegionVector.push_back(G4RegionStore::GetInstance()
				 ->GetRegion("DefaultRegionForTheWorld", false));
      numRegions = 1;
    }

    for( size_t iReg = 0; iReg < numRegions; ++iReg ) 
    {
      const G4Region* curReg = fPAIRegionVector[iReg];
      G4Region* reg = const_cast<G4Region*>(curReg);

      for(size_t jMat = 0; jMat < numOfMat; ++jMat) 
      {
	G4Material* mat = (*theMaterialTable)[jMat];
	const G4MaterialCutsCouple* cutCouple = reg->FindCouple(mat);
	//G4cout << "Couple <" << fCutCouple << G4endl;
	if(cutCouple) 
        {
	  if(fVerbose>0)
	  {
	    G4cout << "Reg <" <<curReg->GetName() << ">  mat <" 
	    << mat->GetName() << ">  fCouple= " 
	    << cutCouple << ", idx= " << cutCouple->GetIndex()
	    <<"  " << p->GetParticleName() 
	    <<", cuts.size() = " << cuts.size() << G4endl;
	  }
	  // check if this couple is not already initialized

	  size_t n = fMaterialCutsCoupleVector.size();

	  G4bool isnew = true;
	  if( 0 < n ) 
          {
	    for(size_t i=0; i<fMaterialCutsCoupleVector.size(); ++i) 
            {
	      if(cutCouple == fMaterialCutsCoupleVector[i]) {
		isnew = false;
		break;
	      }
	    }
	  }
	  // initialise data banks
	  if(isnew) {
	    fMaterialCutsCoupleVector.push_back(cutCouple);
	    G4double deltaCutInKinEnergy = cuts[cutCouple->GetIndex()];
	    fModelData->Initialise(cutCouple, deltaCutInKinEnergy, this);
	  }
	}
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////

void G4PAIPhotModel::InitialiseLocal(const G4ParticleDefinition*, 
				 G4VEmModel* masterModel)
{
  fModelData = static_cast<G4PAIPhotModel*>(masterModel)->GetPAIPhotData();
  fMaterialCutsCoupleVector = static_cast<G4PAIPhotModel*>(masterModel)->GetVectorOfCouples();
  SetElementSelectors( masterModel->GetElementSelectors() );
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotModel::MinEnergyCut(const G4ParticleDefinition*,
		 		      const G4MaterialCutsCouple*)
{
  return fLowestTcut;
}

//////////////////////////////////////////////////////////////////////////////

G4double G4PAIPhotModel::ComputeDEDXPerVolume(const G4Material*,
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

G4double G4PAIPhotModel::CrossSectionPerVolume( const G4Material*,
					    const G4ParticleDefinition* p,
					    G4double kineticEnergy,
					    G4double cutEnergy,
					    G4double maxEnergy  ) 
{
  G4int coupleIndex = FindCoupleIndex(CurrentCouple());

  if(0 > coupleIndex) return 0.0; 

  G4double tmax = std::min(MaxSecondaryEnergy(p, kineticEnergy), maxEnergy);

  if(tmax <= cutEnergy)  return 0.0; 

  G4double scaledTkin = kineticEnergy*fRatio;
  G4double xsc = fChargeSquare*fModelData->CrossSectionPerVolume(coupleIndex, 
							 scaledTkin, 
							 cutEnergy, tmax);
  return xsc;
}

///////////////////////////////////////////////////////////////////////////
//
// It is analog of PostStepDoIt in terms of secondary electron.
//
 
void G4PAIPhotModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
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

  if( maxEnergy <  tmax) tmax = maxEnergy; 
  if( tmin      >= tmax) return; 

  G4ThreeVector direction = dp->GetMomentumDirection();
  G4double scaledTkin     = kineticEnergy*fRatio;
  G4double totalEnergy    = kineticEnergy + fMass;
  G4double totalMomentum  = sqrt(kineticEnergy*(totalEnergy + fMass));
  G4double plRatio        = fModelData->GetPlasmonRatio(coupleIndex, scaledTkin);

  if( G4UniformRand() <= plRatio )
  {
    G4double deltaTkin = fModelData->SamplePostStepPlasmonTransfer(coupleIndex, scaledTkin);

    // G4cout<<"G4PAIPhotModel::SampleSecondaries; dp "<<dp->GetParticleDefinition()->GetParticleName()
    // <<"; Tkin = "<<kineticEnergy/keV<<" keV; transfer = "<<deltaTkin/keV<<" keV "<<G4endl;

    if( deltaTkin <= 0. && fVerbose > 0) 
    {
    G4cout<<"G4PAIPhotModel::SampleSecondary e- deltaTkin = "<<deltaTkin<<G4endl;
    }
    if( deltaTkin <= 0.) return; 

    if( deltaTkin > tmax) deltaTkin = tmax;

    const G4Element* anElement = SelectTargetAtom(matCC,fParticle,kineticEnergy,
                                                  dp->GetLogKineticEnergy());
    G4int Z = G4lrint(anElement->GetZ());
 
    auto deltaRay = new G4DynamicParticle(fElectron,
	    GetAngularDistribution()->SampleDirection(dp, deltaTkin,
		 				      Z, matCC->GetMaterial()),
						      deltaTkin);

    // primary change

    kineticEnergy -= deltaTkin;

    if( kineticEnergy <= 0. ) // kill primary as local? energy deposition
    {
      fParticleChange->SetProposedKineticEnergy(0.0);
      fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy+deltaTkin);
      // fParticleChange->ProposeTrackStatus(fStopAndKill);
      return; 
    }
    else
    {
      G4ThreeVector dir = totalMomentum*direction - deltaRay->GetMomentum();
      direction = dir.unit();
      fParticleChange->SetProposedKineticEnergy(kineticEnergy);
      fParticleChange->SetProposedMomentumDirection(direction);
      vdp->push_back(deltaRay);
    }
  }
  else // secondary X-ray CR photon
  {
    G4double deltaTkin     = fModelData->SamplePostStepPhotonTransfer(coupleIndex, scaledTkin);

    //  G4cout<<"PAIPhotonModel PhotonPostStepTransfer = "<<deltaTkin/keV<<" keV"<<G4endl; 

    if( deltaTkin <= 0. )
    {
      G4cout<<"G4PAIPhotonModel::SampleSecondary gamma deltaTkin = "<<deltaTkin<<G4endl;
    }
    if( deltaTkin <= 0.) return;

    if( deltaTkin >= kineticEnergy ) // stop primary
    {
      deltaTkin = kineticEnergy;
      kineticEnergy = 0.0;
    }
    G4double costheta = 0.; // G4UniformRand(); // VG: ??? for start only
    G4double sintheta = sqrt((1.+costheta)*(1.-costheta));

    //  direction of the 'Cherenkov' photon  
    G4double phi = twopi*G4UniformRand(); 
    G4double dirx = sintheta*cos(phi), diry = sintheta*sin(phi), dirz = costheta;

    G4ThreeVector deltaDirection(dirx,diry,dirz);
    deltaDirection.rotateUz(direction);

    if( kineticEnergy > 0.) // primary change
    {
      kineticEnergy -= deltaTkin;
      fParticleChange->SetProposedKineticEnergy(kineticEnergy);
    }
    else // stop primary, but pass X-ray CR
    {
      // fParticleChange->ProposeLocalEnergyDeposit(deltaTkin);
      fParticleChange->SetProposedKineticEnergy(0.0);
    }
    // create G4DynamicParticle object for photon ray
 
    auto photonRay = new G4DynamicParticle;
    photonRay->SetDefinition( G4Gamma::Gamma() );
    photonRay->SetKineticEnergy( deltaTkin );
    photonRay->SetMomentumDirection(deltaDirection); 

    vdp->push_back(photonRay);
  }
  return;
}

///////////////////////////////////////////////////////////////////////

G4double G4PAIPhotModel::SampleFluctuations(
                         const G4MaterialCutsCouple* matCC,
                         const G4DynamicParticle* aParticle,
                         const G4double, const G4double,
                         const G4double step, const G4double eloss)
{
  // return 0.;
  G4int coupleIndex = FindCoupleIndex(matCC);
  if(0 > coupleIndex) { return eloss; }

  SetParticle(aParticle->GetDefinition());

  
  // G4cout << "G4PAIPhotModel::SampleFluctuations step(mm)= "<< step/mm
  // << "  Eloss(keV)= " << eloss/keV  << " in " 
  // << matCC->GetMaterial()->GetName() << G4endl;
  

  G4double Tkin       = aParticle->GetKineticEnergy();
  G4double scaledTkin = Tkin*fRatio;

  G4double loss = fModelData->SampleAlongStepPhotonTransfer(coupleIndex, Tkin,
						      scaledTkin,
						      step*fChargeSquare);
           loss += fModelData->SampleAlongStepPlasmonTransfer(coupleIndex, Tkin,
						      scaledTkin,
							      step*fChargeSquare);

  
  // G4cout<<"  PAIPhotModel::SampleFluctuations loss = "<<loss/keV<<" keV, on step = "
  // <<step/mm<<" mm"<<G4endl; 
  return loss;

}

//////////////////////////////////////////////////////////////////////
//
// Returns the statistical estimation of the energy loss distribution variance
//


G4double G4PAIPhotModel::Dispersion(const G4Material* material, 
                                    const G4DynamicParticle* aParticle,
 				    const G4double tcut,
 				    const G4double tmax, 
			            const G4double step)
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

G4double G4PAIPhotModel::MaxSecondaryEnergy( const G4ParticleDefinition* p,
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

void G4PAIPhotModel::DefineForRegion(const G4Region* r) 
{
  fPAIRegionVector.push_back(r);
}

//
//
/////////////////////////////////////////////////
