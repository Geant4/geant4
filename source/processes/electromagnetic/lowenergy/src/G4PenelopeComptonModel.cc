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
// $Id: G4PenelopeComptonModel.cc 99415 2016-09-21 09:05:43Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// --------
// 15 Feb 2010   L Pandola  Implementation
// 18 Mar 2010   L Pandola  Removed GetAtomsPerMolecule(), now demanded
//                            to G4PenelopeOscillatorManager
// 01 Feb 2011   L Pandola  Suppress fake energy-violation warning when Auger is
//                            active.
//                          Make sure that fluorescence/Auger is generated only
//                            if above threshold
// 24 May 2011   L Pandola  Renamed (make v2008 as default Penelope)
// 10 Jun 2011   L Pandola  Migrate atomic deexcitation interface
// 09 Oct 2013   L Pandola  Migration to MT
//
#include "G4PenelopeComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4VEMDataSet.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeOscillator.hh"
#include "G4LossTableManager.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopeComptonModel::G4PenelopeComptonModel(const G4ParticleDefinition* part,
					       const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),fParticle(0),
   isInitialised(false),fAtomDeexcitation(0),
   oscManager(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  oscManager = G4PenelopeOscillatorManager::GetOscillatorManager();

  if (part)
    SetParticle(part);

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  fTransitionManager = G4AtomicTransitionManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeComptonModel::~G4PenelopeComptonModel()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::Initialise(const G4ParticleDefinition* part,
					  const G4DataVector&)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeComptonModel::Initialise()" << G4endl;

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  //Issue warning if the AtomicDeexcitation has not been declared
  if (!fAtomDeexcitation)
    {
      G4cout << G4endl;
      G4cout << "WARNING from G4PenelopeComptonModel " << G4endl;
      G4cout << "Atomic de-excitation module is not instantiated, so there will not be ";
      G4cout << "any fluorescence/Auger emission." << G4endl;
      G4cout << "Please make sure this is intended" << G4endl;
    }

  SetParticle(part);

  if (IsMaster() && part == fParticle)
    {

      if (verboseLevel > 0)
	{
	  G4cout << "Penelope Compton model v2008 is initialized " << G4endl
		 << "Energy range: "
		 << LowEnergyLimit() / keV << " keV - "
		 << HighEnergyLimit() / GeV << " GeV";
	}
      //Issue a warning, if the model is going to be used down to a
      //energy which is outside the validity of the model itself
      if (LowEnergyLimit() < fIntrinsicLowEnergyLimit)
	{
	  G4ExceptionDescription ed;
	  ed << "Using the Penelope Compton model outside its intrinsic validity range. "
	     << G4endl;
	  ed << "-> LowEnergyLimit() in process = " << LowEnergyLimit()/keV << "keV " << G4endl;
	  ed << "-> Instrinsic low-energy limit = " << fIntrinsicLowEnergyLimit/keV << "keV "
	     << G4endl;
	  ed << "Result of the simulation have to be taken with care" << G4endl;
	  G4Exception("G4PenelopeComptonModel::Initialise()",
		      "em2100",JustWarning,ed);
	}
    }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::InitialiseLocal(const G4ParticleDefinition* part,
                                                     G4VEmModel *masterModel)
{
  if (verboseLevel > 3)
    G4cout << "Calling  G4PenelopeComptonModel::InitialiseLocal()" << G4endl;

  //
  //Check that particle matches: one might have multiple master models (e.g.
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeComptonModel* theModel =
        static_cast<G4PenelopeComptonModel*> (masterModel);

      //Same verbosity for all workers, as the master
      verboseLevel = theModel->verboseLevel;
    }

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::CrossSectionPerVolume(const G4Material* material,
                                               const G4ParticleDefinition* p,
                                               G4double energy,
                                               G4double,
                                               G4double)
{
  // Penelope model v2008 to calculate the Compton scattering cross section:
  // D. Brusa et al., Nucl. Instrum. Meth. A 379 (1996) 167
  //
  // The cross section for Compton scattering is calculated according to the Klein-Nishina
  // formula for energy > 5 MeV.
  // For E < 5 MeV it is used a parametrization for the differential cross-section dSigma/dOmega,
  // which is integrated numerically in cos(theta), G4PenelopeComptonModel::DifferentialCrossSection().
  // The parametrization includes the J(p)
  // distribution profiles for the atomic shells, that are tabulated from Hartree-Fock calculations
  // from F. Biggs et al., At. Data Nucl. Data Tables 16 (1975) 201
  //
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4PenelopeComptonModel" << G4endl;
  SetupForMaterial(p, material, energy);


  G4double cs = 0;
  //Force null cross-section if below the low-energy edge of the table
  if (energy < LowEnergyLimit())
    return cs;

  //Retrieve the oscillator table for this material
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableCompton(material);

  if (energy < 5*MeV) //explicit calculation for E < 5 MeV
    {
      size_t numberOfOscillators = theTable->size();
      for (size_t i=0;i<numberOfOscillators;i++)
	{
	  G4PenelopeOscillator* theOsc = (*theTable)[i];
	  //sum contributions coming from each oscillator
	  cs += OscillatorTotalCrossSection(energy,theOsc);
	}
    }
  else //use Klein-Nishina for E>5 MeV
    cs = KleinNishinaCrossSection(energy,material);

  //cross sections are in units of pi*classic_electr_radius^2
  cs *= pi*classic_electr_radius*classic_electr_radius;

  //Now, cs is the cross section *per molecule*, let's calculate the
  //cross section per volume

  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double atPerMol =  oscManager->GetAtomsPerMolecule(material);

  if (verboseLevel > 3)
    G4cout << "Material " << material->GetName() << " has " << atPerMol <<
      "atoms per molecule" << G4endl;

  G4double moleculeDensity = 0.;

  if (atPerMol)
    moleculeDensity = atomDensity/atPerMol;

  G4double csvolume = cs*moleculeDensity;

  if (verboseLevel > 2)
    G4cout << "Compton mean free path at " << energy/keV << " keV for material " <<
            material->GetName() << " = " << (1./csvolume)/mm << " mm" << G4endl;
  return csvolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//This is a dummy method. Never inkoved by the tracking, it just issues
//a warning if one tries to get Cross Sections per Atom via the
//G4EmCalculator.
G4double G4PenelopeComptonModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                             G4double,
                                                             G4double,
                                                             G4double,
                                                             G4double,
                                                             G4double)
{
  G4cout << "*** G4PenelopeComptonModel -- WARNING ***" << G4endl;
  G4cout << "Penelope Compton model v2008 does not calculate cross section _per atom_ " << G4endl;
  G4cout << "so the result is always zero. For physics values, please invoke " << G4endl;
  G4cout << "GetCrossSectionPerVolume() or GetMeanFreePath() via the G4EmCalculator" << G4endl;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					       const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{

  // Penelope model v2008 to sample the Compton scattering final state.
  // D. Brusa et al., Nucl. Instrum. Meth. A 379 (1996) 167
  // The model determines also the original shell from which the electron is expelled,
  // in order to produce fluorescence de-excitation (from G4DeexcitationManager)
  //
  // The final state for Compton scattering is calculated according to the Klein-Nishina
  // formula for energy > 5 MeV. In this case, the Doppler broadening is negligible and
  // one can assume that the target electron is at rest.
  // For E < 5 MeV it is used the parametrization for the differential cross-section dSigma/dOmega,
  // to sample the scattering angle and the energy of the emerging electron, which is
  // G4PenelopeComptonModel::DifferentialCrossSection(). The rejection method is
  // used to sample cos(theta). The efficiency increases monotonically with photon energy and is
  // nearly independent on the Z; typical values are 35%, 80% and 95% for 1 keV, 1 MeV and 10 MeV,
  // respectively.
  // The parametrization includes the J(p) distribution profiles for the atomic shells, that are
  // tabulated
  // from Hartree-Fock calculations from F. Biggs et al., At. Data Nucl. Data Tables 16 (1975) 201.
  // Doppler broadening is included.
  //

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4PenelopeComptonModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  // do nothing below the threshold
  // should never get here because the XS is zero below the limit
  if(photonEnergy0 < LowEnergyLimit())
    return;

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();
  const G4Material* material = couple->GetMaterial();

  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableCompton(material);

  const G4int nmax = 64;
  G4double rn[nmax]={0.0};
  G4double pac[nmax]={0.0};

  G4double S=0.0;
  G4double epsilon = 0.0;
  G4double cosTheta = 1.0;
  G4double hartreeFunc = 0.0;
  G4double oscStren = 0.0;
  size_t numberOfOscillators = theTable->size();
  size_t targetOscillator = 0;
  G4double ionEnergy = 0.0*eV;

  G4double ek = photonEnergy0/electron_mass_c2;
  G4double ek2 = 2.*ek+1.0;
  G4double eks = ek*ek;
  G4double ek1 = eks-ek2-1.0;

  G4double taumin = 1.0/ek2;
  G4double a1 = std::log(ek2);
  G4double a2 = a1+2.0*ek*(1.0+ek)/(ek2*ek2);

  G4double TST = 0;
  G4double tau = 0.;

  //If the incoming photon is above 5 MeV, the quicker approach based on the
  //pure Klein-Nishina formula is used
  if (photonEnergy0 > 5*MeV)
    {
      do{
	do{
	  if ((a2*G4UniformRand()) < a1)
	    tau = std::pow(taumin,G4UniformRand());
	  else
	    tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	  //rejection function
	  TST = (1.0+tau*(ek1+tau*(ek2+tau*eks)))/(eks*tau*(1.0+tau*tau));
	}while (G4UniformRand()> TST);
	epsilon=tau;
	cosTheta = 1.0 - (1.0-tau)/(ek*tau);

	//Target shell electrons
	TST = oscManager->GetTotalZ(material)*G4UniformRand();
	targetOscillator = numberOfOscillators-1; //last level
	S=0.0;
	G4bool levelFound = false;
	for (size_t j=0;j<numberOfOscillators && !levelFound; j++)
	  {
	    S += (*theTable)[j]->GetOscillatorStrength();
	    if (S > TST)
	      {
		targetOscillator = j;
		levelFound = true;
	      }
	  }
	//check whether the level is valid
	ionEnergy = (*theTable)[targetOscillator]->GetIonisationEnergy();
      }while((epsilon*photonEnergy0-photonEnergy0+ionEnergy) >0);
    }
  else //photonEnergy0 < 5 MeV
    {
      //Incoherent scattering function for theta=PI
      G4double s0=0.0;
      G4double pzomc=0.0;
      G4double rni=0.0;
      G4double aux=0.0;
      for (size_t i=0;i<numberOfOscillators;i++)
	{
	  ionEnergy = (*theTable)[i]->GetIonisationEnergy();
	  if (photonEnergy0 > ionEnergy)
	    {
	      G4double aux2 = photonEnergy0*(photonEnergy0-ionEnergy)*2.0;
	      hartreeFunc = (*theTable)[i]->GetHartreeFactor();
	      oscStren = (*theTable)[i]->GetOscillatorStrength();
	      pzomc = hartreeFunc*(aux2-electron_mass_c2*ionEnergy)/
		(electron_mass_c2*std::sqrt(2.0*aux2+ionEnergy*ionEnergy));
	      if (pzomc > 0)
		rni = 1.0-0.5*G4Exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*
				       (std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
	      else
		rni = 0.5*G4Exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*
				   (std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
	      s0 += oscStren*rni;
	    }
	}
      //Sampling tau
      G4double cdt1 = 0.;
      do
	{
	  if ((G4UniformRand()*a2) < a1)
	    tau = std::pow(taumin,G4UniformRand());
	  else
	    tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	  cdt1 = (1.0-tau)/(ek*tau);
	  //Incoherent scattering function
	  S = 0.;
	  for (size_t i=0;i<numberOfOscillators;i++)
	    {
	      ionEnergy = (*theTable)[i]->GetIonisationEnergy();
	      if (photonEnergy0 > ionEnergy) //sum only on excitable levels
		{
		  aux = photonEnergy0*(photonEnergy0-ionEnergy)*cdt1;
		  hartreeFunc = (*theTable)[i]->GetHartreeFactor();
		  oscStren = (*theTable)[i]->GetOscillatorStrength();
		  pzomc = hartreeFunc*(aux-electron_mass_c2*ionEnergy)/
		    (electron_mass_c2*std::sqrt(2.0*aux+ionEnergy*ionEnergy));
		  if (pzomc > 0)
		    rn[i] = 1.0-0.5*G4Exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*
					     (std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
		  else
		    rn[i] = 0.5*G4Exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*
					 (std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
		  S += oscStren*rn[i];
		  pac[i] = S;
		}
	      else
		pac[i] = S-1e-6;
	    }
	  //Rejection function
	  TST = S*(1.0+tau*(ek1+tau*(ek2+tau*eks)))/(eks*tau*(1.0+tau*tau));
	}while ((G4UniformRand()*s0) > TST);

      cosTheta = 1.0 - cdt1;
      G4double fpzmax=0.0,fpz=0.0;
      G4double A=0.0;
      //Target electron shell
      do
	{
	  do
	    {
	      TST = S*G4UniformRand();
	      targetOscillator = numberOfOscillators-1; //last level
	      G4bool levelFound = false;
	      for (size_t i=0;i<numberOfOscillators && !levelFound;i++)
		{
		  if (pac[i]>TST)
		    {
		      targetOscillator = i;
		      levelFound = true;
		    }
		}
	      A = G4UniformRand()*rn[targetOscillator];
	      hartreeFunc = (*theTable)[targetOscillator]->GetHartreeFactor();
	      oscStren = (*theTable)[targetOscillator]->GetOscillatorStrength();
	      if (A < 0.5)
		pzomc = (std::sqrt(0.5)-std::sqrt(0.5-std::log(2.0*A)))/
		  (std::sqrt(2.0)*hartreeFunc);
	      else
		pzomc = (std::sqrt(0.5-std::log(2.0-2.0*A))-std::sqrt(0.5))/
		  (std::sqrt(2.0)*hartreeFunc);
	    } while (pzomc < -1);

	  // F(EP) rejection
	  G4double XQC = 1.0+tau*(tau-2.0*cosTheta);
	  G4double AF = std::sqrt(XQC)*(1.0+tau*(tau-cosTheta)/XQC);
	  if (AF > 0)
	    fpzmax = 1.0+AF*0.2;
	  else
	    fpzmax = 1.0-AF*0.2;
	  fpz = 1.0+AF*std::max(std::min(pzomc,0.2),-0.2);
	}while ((fpzmax*G4UniformRand())>fpz);

      //Energy of the scattered photon
      G4double T = pzomc*pzomc;
      G4double b1 = 1.0-T*tau*tau;
      G4double b2 = 1.0-T*tau*cosTheta;
      if (pzomc > 0.0)
	epsilon = (tau/b1)*(b2+std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
      else
	epsilon = (tau/b1)*(b2-std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
    } //energy < 5 MeV

  //Ok, the kinematics has been calculated.
  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  // Update G4VParticleChange for the scattered photon
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;

  G4double photonEnergy1 = epsilon * photonEnergy0;

  if (photonEnergy1 > 0.)
    fParticleChange->SetProposedKineticEnergy(photonEnergy1) ;
  else
  {
    fParticleChange->SetProposedKineticEnergy(0.) ;
    fParticleChange->ProposeTrackStatus(fStopAndKill);
  }

  //Create scattered electron
  G4double diffEnergy = photonEnergy0*(1-epsilon);
  ionEnergy = (*theTable)[targetOscillator]->GetIonisationEnergy();

  G4double Q2 =
    photonEnergy0*photonEnergy0+photonEnergy1*(photonEnergy1-2.0*photonEnergy0*cosTheta);
  G4double cosThetaE = 0.; //scattering angle for the electron

  if (Q2 > 1.0e-12)
    cosThetaE = (photonEnergy0-photonEnergy1*cosTheta)/std::sqrt(Q2);
  else
    cosThetaE = 1.0;
  G4double sinThetaE = std::sqrt(1-cosThetaE*cosThetaE);

  //Now, try to handle fluorescence
  //Notice: merged levels are indicated with Z=0 and flag=30
  G4int shFlag = (*theTable)[targetOscillator]->GetShellFlag();
  G4int Z = (G4int) (*theTable)[targetOscillator]->GetParentZ();

  //initialize here, then check photons created by Atomic-Deexcitation, and the final state e-
  G4double bindingEnergy = 0.*eV;
  const G4AtomicShell* shell = 0;

  //Real level
  if (Z > 0 && shFlag<30)
    {
      shell = fTransitionManager->Shell(Z,shFlag-1);
      bindingEnergy = shell->BindingEnergy();
    }

  G4double ionEnergyInPenelopeDatabase = ionEnergy;
  //protection against energy non-conservation
  ionEnergy = std::max(bindingEnergy,ionEnergyInPenelopeDatabase);

  //subtract the excitation energy. If not emitted by fluorescence
  //the ionization energy is deposited as local energy deposition
  G4double eKineticEnergy = diffEnergy - ionEnergy;
  G4double localEnergyDeposit = ionEnergy;
  G4double energyInFluorescence = 0.; //testing purposes only
  G4double energyInAuger = 0; //testing purposes

  if (eKineticEnergy < 0)
    {
      //It means that there was some problem/mismatch between the two databases.
      //Try to make it work
      //In this case available Energy (diffEnergy) < ionEnergy
      //Full residual energy is deposited locally
      localEnergyDeposit = diffEnergy;
      eKineticEnergy = 0.0;
    }

  //the local energy deposit is what remains: part of this may be spent for fluorescence.
  //Notice: shell might be nullptr (invalid!) if shFlag=30. Must be protected
  //Now, take care of fluorescence, if required
  if (fAtomDeexcitation && shell)
    {
      G4int index = couple->GetIndex();
      if (fAtomDeexcitation->CheckDeexcitationActiveRegion(index))
	{
	  size_t nBefore = fvect->size();
	  fAtomDeexcitation->GenerateParticles(fvect,shell,Z,index);
	  size_t nAfter = fvect->size();

	  if (nAfter > nBefore) //actual production of fluorescence
	    {
	      for (size_t j=nBefore;j<nAfter;j++) //loop on products
		{
		  G4double itsEnergy = ((*fvect)[j])->GetKineticEnergy();
		  localEnergyDeposit -= itsEnergy;
		  if (((*fvect)[j])->GetParticleDefinition() == G4Gamma::Definition())
		    energyInFluorescence += itsEnergy;
		  else if (((*fvect)[j])->GetParticleDefinition() == G4Electron::Definition())
		    energyInAuger += itsEnergy;

		}
	    }

	}
    }


  /*
  if(DeexcitationFlag() && Z > 5 && shellId>0) {

    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();

    size_t index = couple->GetIndex();
    G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[index];
    G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];

    // Generation of fluorescence
    // Data in EADL are available only for Z > 5
    // Protection to avoid generating photons in the unphysical case of
    // shell binding energy > photon energy
    if (localEnergyDeposit > cutg || localEnergyDeposit > cute)
      {
	G4DynamicParticle* aPhoton;
	deexcitationManager.SetCutForSecondaryPhotons(cutg);
	deexcitationManager.SetCutForAugerElectrons(cute);

	photonVector = deexcitationManager.GenerateParticles(Z,shellId);
	if(photonVector)
	  {
	    size_t nPhotons = photonVector->size();
	    for (size_t k=0; k<nPhotons; k++)
	      {
		aPhoton = (*photonVector)[k];
		if (aPhoton)
		  {
		    G4double itsEnergy = aPhoton->GetKineticEnergy();
		    G4bool keepIt = false;
		    if (itsEnergy <= localEnergyDeposit)
		      {
			//check if good!
			if(aPhoton->GetDefinition() == G4Gamma::Gamma()
			   && itsEnergy >= cutg)
			  {
			    keepIt = true;
			    energyInFluorescence += itsEnergy;
			  }
			if (aPhoton->GetDefinition() == G4Electron::Electron() &&
			    itsEnergy >= cute)
			  {
			    energyInAuger += itsEnergy;
			    keepIt = true;
			  }
		      }
		    //good secondary, register it
		    if (keepIt)
		      {
			localEnergyDeposit -= itsEnergy;
			fvect->push_back(aPhoton);
		      }
		    else
		      {
			delete aPhoton;
			(*photonVector)[k] = 0;
		      }
		  }
	      }
	    delete photonVector;
	  }
      }
  }
  */


  //Always produce explicitely the electron
  G4DynamicParticle* electron = 0;

  G4double xEl = sinThetaE * std::cos(phi+pi);
  G4double yEl = sinThetaE * std::sin(phi+pi);
  G4double zEl = cosThetaE;
  G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
  eDirection.rotateUz(photonDirection0);
  electron = new G4DynamicParticle (G4Electron::Electron(),
				    eDirection,eKineticEnergy) ;
  fvect->push_back(electron);


  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING-"
	     << "G4PenelopeComptonModel::SampleSecondaries - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit=0.;
    }
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  G4double electronEnergy = 0.;
  if (electron)
    electronEnergy = eKineticEnergy;
  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeCompton" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Scattered photon: " << photonEnergy1/keV << " keV" << G4endl;
      G4cout << "Scattered electron " << electronEnergy/keV << " keV" << G4endl;
      if (energyInFluorescence)
	G4cout << "Fluorescence x-rays: " << energyInFluorescence/keV << " keV" << G4endl;
      if (energyInAuger)
	G4cout << "Auger electrons: " << energyInAuger/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (photonEnergy1+electronEnergy+energyInFluorescence+
					  localEnergyDeposit+energyInAuger)/keV <<
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(photonEnergy1+
				      electronEnergy+energyInFluorescence+
				      localEnergyDeposit+energyInAuger-photonEnergy0);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeCompton: problem with energy conservation: " <<
	  (photonEnergy1+electronEnergy+energyInFluorescence+energyInAuger+localEnergyDeposit)/keV <<
	  " keV (final) vs. " <<
	  photonEnergy0/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::DifferentialCrossSection(G4double cosTheta,G4double energy,
							    G4PenelopeOscillator* osc)
{
  //
  // Penelope model v2008. Single differential cross section *per electron*
  // for photon Compton scattering by
  // electrons in the given atomic oscillator, differential in the direction of the
  // scattering photon. This is in units of pi*classic_electr_radius**2
  //
  // D. Brusa et al., Nucl. Instrum. Meth. A 379 (1996) 167
  // The parametrization includes the J(p) distribution profiles for the atomic shells,
  // that are tabulated from Hartree-Fock calculations
  // from F. Biggs et al., At. Data Nucl. Data Tables 16 (1975) 201
  //
  G4double ionEnergy = osc->GetIonisationEnergy();
  G4double harFunc = osc->GetHartreeFactor();

  static const G4double k2 = std::sqrt(2.);
  static const G4double k1 = 1./k2;

  if (energy < ionEnergy)
    return 0;

  //energy of the Compton line
  G4double cdt1 = 1.0-cosTheta;
  G4double EOEC = 1.0+(energy/electron_mass_c2)*cdt1;
  G4double ECOE = 1.0/EOEC;

  //Incoherent scattering function (analytical profile)
  G4double aux = energy*(energy-ionEnergy)*cdt1;
  G4double Pzimax =
    (aux - electron_mass_c2*ionEnergy)/(electron_mass_c2*std::sqrt(2*aux+ionEnergy*ionEnergy));
  G4double sia = 0.0;
  G4double x = harFunc*Pzimax;
  if (x > 0)
    sia = 1.0-0.5*G4Exp(0.5-(k1+k2*x)*(k1+k2*x));
  else
    sia = 0.5*G4Exp(0.5-(k1-k2*x)*(k1-k2*x));

  //1st order correction, integral of Pz times the Compton profile.
  //Calculated approximately using a free-electron gas profile
  G4double pf = 3.0/(4.0*harFunc);
  if (std::fabs(Pzimax) < pf)
    {
      G4double QCOE2 = 1.0+ECOE*ECOE-2.0*ECOE*cosTheta;
      G4double p2 = Pzimax*Pzimax;
      G4double dspz = std::sqrt(QCOE2)*
	(1.0+ECOE*(ECOE-cosTheta)/QCOE2)*harFunc
	*0.25*(2*p2-(p2*p2)/(pf*pf)-(pf*pf));
      sia += std::max(dspz,-1.0*sia);
    }

  G4double XKN = EOEC+ECOE-1.0+cosTheta*cosTheta;

  //Differential cross section (per electron, in units of pi*classic_electr_radius**2)
  G4double diffCS = ECOE*ECOE*XKN*sia;

  return diffCS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::OscillatorTotalCrossSection(G4double energy,G4PenelopeOscillator* osc)
{
  //Total cross section (integrated) for the given oscillator in units of
  //pi*classic_electr_radius^2

  //Integrate differential cross section for each oscillator
  G4double stre = osc->GetOscillatorStrength();

  // here one uses the  using the 20-point
  // Gauss quadrature method with an adaptive bipartition scheme
  const G4int npoints=10;
  const G4int ncallsmax=20000;
  const G4int nst=256;
  static const G4double Abscissas[10] = {7.652651133497334e-02,2.2778585114164508e-01,3.7370608871541956e-01,
				  5.1086700195082710e-01,6.3605368072651503e-01,7.4633190646015079e-01,
				  8.3911697182221882e-01,9.1223442825132591e-01,9.6397192727791379e-01,
				  9.9312859918509492e-01};
  static const G4double Weights[10] = {1.5275338713072585e-01,1.4917298647260375e-01,1.4209610931838205e-01,
				1.3168863844917663e-01,1.1819453196151842e-01,1.0193011981724044e-01,
				8.3276741576704749e-02,6.2672048334109064e-02,4.0601429800386941e-02,
				1.7614007139152118e-02};

  G4double MaxError = 1e-5;
  //Error control
  G4double Ctol = std::min(std::max(MaxError,1e-13),1e-02);
  G4double Ptol = 0.01*Ctol;
  G4double Err=1e35;

  //Gauss integration from -1 to 1
  G4double LowPoint = -1.0;
  G4double HighPoint = 1.0;

  G4double h=HighPoint-LowPoint;
  G4double sumga=0.0;
  G4double a=0.5*(HighPoint-LowPoint);
  G4double b=0.5*(HighPoint+LowPoint);
  G4double c=a*Abscissas[0];
  G4double d= Weights[0]*
    (DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));
  for (G4int i=2;i<=npoints;i++)
    {
      c=a*Abscissas[i-1];
      d += Weights[i-1]*
	(DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));
    }
  G4int icall = 2*npoints;
  G4int LH=1;
  G4double S[nst],x[nst],sn[nst],xrn[nst];
  S[0]=d*a;
  x[0]=LowPoint;

  G4bool loopAgain = true;

  //Adaptive bipartition scheme
  do{
    G4double h0=h;
    h=0.5*h; //bipartition
    G4double sumr=0;
    G4int LHN=0;
    G4double si,xa,xb,xc;
    for (G4int i=1;i<=LH;i++){
      si=S[i-1];
      xa=x[i-1];
      xb=xa+h;
      xc=xa+h0;
      a=0.5*(xb-xa);
      b=0.5*(xb+xa);
      c=a*Abscissas[0];
      G4double dLocal = Weights[0]*
	(DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));

      for (G4int j=1;j<npoints;j++)
	{
	  c=a*Abscissas[j];
	  dLocal += Weights[j]*
	    (DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));
	}
      G4double s1=dLocal*a;
      a=0.5*(xc-xb);
      b=0.5*(xc+xb);
      c=a*Abscissas[0];
      dLocal=Weights[0]*
	(DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));

      for (G4int j=1;j<npoints;j++)
	{
	  c=a*Abscissas[j];
	  dLocal += Weights[j]*
	    (DifferentialCrossSection(b+c,energy,osc)+DifferentialCrossSection(b-c,energy,osc));
	}
      G4double s2=dLocal*a;
      icall=icall+4*npoints;
      G4double s12=s1+s2;
      if (std::abs(s12-si)<std::max(Ptol*std::abs(s12),1e-35))
	sumga += s12;
      else
	{
	  sumr += s12;
	  LHN += 2;
	  sn[LHN-1]=s2;
	  xrn[LHN-1]=xb;
	  sn[LHN-2]=s1;
	  xrn[LHN-2]=xa;
	}

      if (icall>ncallsmax || LHN>nst)
	{
	  G4cout << "G4PenelopeComptonModel: " << G4endl;
	  G4cout << "LowPoint: " << LowPoint << ", High Point: " << HighPoint << G4endl;
	  G4cout << "Tolerance: " << MaxError << G4endl;
	  G4cout << "Calls: " << icall << ", Integral: " << sumga << ", Error: " << Err << G4endl;
	  G4cout << "Number of open subintervals: " << LHN << G4endl;
	  G4cout << "WARNING: the required accuracy has not been attained" << G4endl;
	  loopAgain = false;
	}
    }
    Err=std::abs(sumr)/std::max(std::abs(sumr+sumga),1e-35);
    if (Err < Ctol || LHN == 0)
      loopAgain = false; //end of cycle
    LH=LHN;
    for (G4int i=0;i<LH;i++)
      {
	S[i]=sn[i];
	x[i]=xrn[i];
      }
  }while(Ctol < 1.0 && loopAgain);


  G4double xs = stre*sumga;

  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::KleinNishinaCrossSection(G4double energy,
						   const G4Material* material)
{
  // use Klein-Nishina formula
  // total cross section in units of pi*classic_electr_radius^2

  G4double cs = 0;

  G4double ek =energy/electron_mass_c2;
  G4double eks = ek*ek;
  G4double ek2 = 1.0+ek+ek;
  G4double ek1 = eks-ek2-1.0;

  G4double t0 = 1.0/ek2;
  G4double csl = 0.5*eks*t0*t0+ek2*t0+ek1*std::log(t0)-(1.0/t0);

  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableCompton(material);

  for (size_t i=0;i<theTable->size();i++)
    {
      G4PenelopeOscillator* theOsc = (*theTable)[i];
      G4double ionEnergy = theOsc->GetIonisationEnergy();
      G4double tau=(energy-ionEnergy)/energy;
      if (tau > t0)
	{
	  G4double csu = 0.5*eks*tau*tau+ek2*tau+ek1*std::log(tau)-(1.0/tau);
	  G4double stre = theOsc->GetOscillatorStrength();

	  cs += stre*(csu-csl);
	}
    }

  cs /= (ek*eks);

  return cs;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeComptonModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;
  }
}
