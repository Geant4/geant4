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
// $Id: G4PenelopeComptonModel.cc,v 1.11 2010-12-01 15:20:26 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 02 Oct 2008   L Pandola    Migration from process to model 
// 28 Oct 2008   L Pandola    Treat the database data from Penelope according to the 
//                            original model, namely merging levels below 15 eV in 
//                            a single one. Still, it is not fully compliant with the 
//                            original Penelope model, because plasma excitation is not 
//                            considered.
// 22 Nov 2008   L Pandola    Make unit of measurements explicit for binding energies 
//			      that are read from the external files.
// 24 Nov 2008   L Pandola    Find a cleaner way to delete vectors.
// 16 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - do not apply production threshold on level of the model
// 21 Oct 2009   L Pandola    Remove un-necessary fUseAtomicDeexcitation flag - now managed by
//                            G4VEmModel::DeexcitationFlag()
//                            Add ActivateAuger() method
//

#include "G4PenelopeComptonModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4VEMDataSet.hh"
#include "G4PhysicsTable.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PenelopeIntegrator.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopeComptonModel::G4PenelopeComptonModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),ionizationEnergy(0),hartreeFunction(0),
   occupationNumber(0),isInitialised(false)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  energyForIntegration = 0.0;
  ZForIntegration = 1;

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //by default, the model will use atomic deexcitation
  SetDeexcitationFlag(true);
  ActivateAuger(false);

  //These vectors do not change when materials or cut change.
  //Therefore I can read it at the constructor
  ionizationEnergy = new std::map<G4int,G4DataVector*>;
  hartreeFunction  = new std::map<G4int,G4DataVector*>;
  occupationNumber = new std::map<G4int,G4DataVector*>;

  ReadData(); //Read data from file

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeComptonModel::~G4PenelopeComptonModel()
{  
  std::map <G4int,G4DataVector*>::iterator i;
  if (ionizationEnergy)
    {
      for (i=ionizationEnergy->begin();i != ionizationEnergy->end();i++)
	if (i->second) delete i->second;
      delete ionizationEnergy;
    }
  if (hartreeFunction)
    {
      for (i=hartreeFunction->begin();i != hartreeFunction->end();i++)
	if (i->second) delete i->second;
      delete hartreeFunction;
    }
  if (occupationNumber)
    {
      for (i=occupationNumber->begin();i != occupationNumber->end();i++)
	if (i->second) delete i->second;
      delete occupationNumber;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::Initialise(const G4ParticleDefinition* particle,
					const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeComptonModel::Initialise()" << G4endl;

  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 0) {
    G4cout << "Penelope Compton model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / keV << " keV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  // Penelope model to calculate the Compton scattering cross section:
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
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopeComptonModel" << G4endl;

  G4int iZ = (G4int) Z;
  G4double cs=0.0;
  energyForIntegration=energy; 
  ZForIntegration = iZ;
  if (energy< 5*MeV)
    {
      // numerically integrate differential cross section dSigma/dOmega
      G4PenelopeIntegrator<G4PenelopeComptonModel,G4double (G4PenelopeComptonModel::*)(G4double)> 
	theIntegrator;
      cs = theIntegrator.Calculate(this,&G4PenelopeComptonModel::DifferentialCrossSection,-1.0,1.0,1e-05);
    }
  else 
    {
      // use Klein-Nishina formula
      G4double ki=energy/electron_mass_c2;
      G4double ki3=ki*ki;
      G4double ki2=1.0+2*ki;
      G4double ki1=ki3-ki2-1.0;
      G4double t0=1.0/(ki2);
      G4double csl = 0.5*ki3*t0*t0+ki2*t0+ki1*std::log(t0)-(1.0/t0);
      G4int nosc = occupationNumber->find(iZ)->second->size();
      for (G4int i=0;i<nosc;i++)
	{
	  G4double ionEnergy = (*(ionizationEnergy->find(iZ)->second))[i];
	  G4double tau=(energy-ionEnergy)/energy;
	  if (tau > t0)
	    {
	      G4double csu = 0.5*ki3*tau*tau+ki2*tau+ki1*std::log(tau)-(1.0/tau);
	      G4int f = (G4int) (*(occupationNumber->find(iZ)->second))[i];
	      cs = cs + f*(csu-csl);
	    }
	}
      cs=pi*classic_electr_radius*classic_electr_radius*cs/(ki*ki3);
    }
  
  if (verboseLevel > 2)
    G4cout << "Compton cross Section at " << energy/keV << " keV for Z=" << Z << 
      " = " << cs/barn << " barn" << G4endl;
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  
  // Penelope model to sample the Compton scattering final state.
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

  if (photonEnergy0 <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  if (verboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;
  // atom can be selected effitiantly if element selectors are initialised 
  const G4Element* anElement = 
    SelectRandomAtom(couple,G4Gamma::GammaDefinition(),photonEnergy0);
  G4int Z = (G4int) anElement->GetZ();
  if (verboseLevel > 2)
    G4cout << "Selected " << anElement->GetName() << G4endl;

  const G4int nmax = 64;
  G4double rn[nmax]={0.0};
  G4double pac[nmax]={0.0};
  
  G4double ki,ki1,ki2,ki3,taumin,a1,a2;
  G4double tau,TST;
  G4double S=0.0;
  G4double epsilon,cosTheta;
  G4double harFunc = 0.0;
  G4int occupNb= 0;
  G4double ionEnergy=0.0;
  G4int nosc = occupationNumber->find(Z)->second->size();
  G4int iosc = nosc;
  ki = photonEnergy0/electron_mass_c2;
  ki2 = 2*ki+1.0;
  ki3 = ki*ki;
  ki1 = ki3-ki2-1.0;
  taumin = 1.0/ki2;
  a1 = std::log(ki2);
  a2 = a1+2.0*ki*(1.0+ki)/(ki2*ki2);
  //If the incoming photon is above 5 MeV, the quicker approach based on the 
  //pure Klein-Nishina formula is used
  if (photonEnergy0 > 5*MeV)
    {
      do{
	do{
	  if ((a2*G4UniformRand()) < a1)
	    {
	      tau = std::pow(taumin,G4UniformRand());
	    }
	  else
	    {
	      tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  //rejection function
	  TST = (1+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));
	}while (G4UniformRand()> TST);
	epsilon=tau;
	cosTheta = 1.0 - (1.0-tau)/(ki*tau);
	//Target shell electrons
	TST = Z*G4UniformRand();
	iosc = nosc;
	S=0.0;
	for (G4int j=0;j<nosc;j++)
	  {
	    occupNb = (G4int) (*(occupationNumber->find(Z)->second))[j];
	    S = S + occupNb;
	    if (S > TST) iosc = j;
	    if (S > TST) break; 
	  }
	ionEnergy = (*(ionizationEnergy->find(Z)->second))[iosc];
      }while((epsilon*photonEnergy0-photonEnergy0+ionEnergy) >0);
    }
  else //photonEnergy0 < 5 MeV
    {
      //Incoherent scattering function for theta=PI
      G4double s0=0.0;
      G4double pzomc=0.0,rni=0.0;
      G4double aux=0.0;
      for (G4int i=0;i<nosc;i++){
	ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
	if (photonEnergy0 > ionEnergy)
	  {
	    G4double aux = photonEnergy0*(photonEnergy0-ionEnergy)*2.0;
	    harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
	    occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
	    pzomc = harFunc*(aux-electron_mass_c2*ionEnergy)/
	       (electron_mass_c2*std::sqrt(2.0*aux+ionEnergy*ionEnergy));
	    if (pzomc > 0) 
	      {
		rni = 1.0-0.5*std::exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*
				       (std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
	      }
	    else
	      {
		rni = 0.5*std::exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*
				   (std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
	      }
	    s0 = s0 + occupNb*rni;
	  }
      }
      
      //Sampling tau
      G4double cdt1;
      do
	{
	  if ((G4UniformRand()*a2) < a1)
	    {
	      tau = std::pow(taumin,G4UniformRand());
	    }
	  else
	    {
	      tau = std::sqrt(1.0+G4UniformRand()*(taumin*taumin-1.0));
	    }
	  cdt1 = (1.0-tau)/(ki*tau);
	  S=0.0;
	  //Incoherent scattering function
	  for (G4int i=0;i<nosc;i++){
	    ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
	    if (photonEnergy0 > ionEnergy) //sum only on excitable levels
	      {
		aux = photonEnergy0*(photonEnergy0-ionEnergy)*cdt1;
		harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
		occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
		pzomc = harFunc*(aux-electron_mass_c2*ionEnergy)/
		  (electron_mass_c2*std::sqrt(2.0*aux+ionEnergy*ionEnergy));
		if (pzomc > 0) 
		  {
		    rn[i] = 1.0-0.5*std::exp(0.5-(std::sqrt(0.5)+std::sqrt(2.0)*pzomc)*
					     (std::sqrt(0.5)+std::sqrt(2.0)*pzomc));
		  }
		else
		  {
		    rn[i] = 0.5*std::exp(0.5-(std::sqrt(0.5)-std::sqrt(2.0)*pzomc)*
					 (std::sqrt(0.5)-std::sqrt(2.0)*pzomc));
		  }
		S = S + occupNb*rn[i];
		pac[i] = S;
	      }
	    else
	      {
		pac[i] = S-(1e-06);
	      }
	  }
	  //Rejection function
	  TST = S*(1.0+tau*(ki1+tau*(ki2+tau*ki3)))/(ki3*tau*(1.0+tau*tau));  
	}while ((G4UniformRand()*s0) > TST);
      //Target electron shell
      cosTheta = 1.0 - cdt1;
      G4double fpzmax=0.0,fpz=0.0;
      G4double A=0.0;
      do
	{
	  do
	    {
	      TST =S*G4UniformRand();
	      iosc=nosc;
	      for (G4int i=0;i<nosc;i++){
		if (pac[i]>TST) iosc = i;
		if (pac[i]>TST) break; 
	      }
	      A = G4UniformRand()*rn[iosc];
	      harFunc = (*(hartreeFunction->find(Z)->second))[iosc]/fine_structure_const;
	      occupNb = (G4int) (*(occupationNumber->find(Z)->second))[iosc];
	      if (A < 0.5) {
		pzomc = (std::sqrt(0.5)-std::sqrt(0.5-std::log(2.0*A)))/
		  (std::sqrt(2.0)*harFunc);
	      }
	      else
		{
		  pzomc = (std::sqrt(0.5-std::log(2.0-2.0*A))-std::sqrt(0.5))/
		    (std::sqrt(2.0)*harFunc);
		}
	    } while (pzomc < -1);
	  // F(EP) rejection
	  G4double XQC = 1.0+tau*(tau-2.0*cosTheta);
	  G4double AF = std::sqrt(XQC)*(1.0+tau*(tau-cosTheta)/XQC);
	  if (AF > 0) {
	    fpzmax = 1.0+AF*0.2;
	  }
	  else
	    {
	      fpzmax = 1.0-AF*0.2;
	    }
	  fpz = 1.0+AF*std::max(std::min(pzomc,0.2),-0.2);
	}while ((fpzmax*G4UniformRand())>fpz);
  
      //Energy of the scattered photon
      G4double T = pzomc*pzomc;
      G4double b1 = 1.0-T*tau*tau;
      G4double b2 = 1.0-T*tau*cosTheta;
      if (pzomc > 0.0)
	{
	  epsilon = (tau/b1)*(b2+std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
	}
      else
	{
	  epsilon = (tau/b1)*(b2-std::sqrt(std::abs(b2*b2-b1*(1.0-T))));
	}
    }
  

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
  {
    fParticleChange->SetProposedKineticEnergy(photonEnergy1) ;
  }
  else
  {
    fParticleChange->SetProposedKineticEnergy(0.) ;
    fParticleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  //Create scattered electron
  G4double diffEnergy = photonEnergy0*(1-epsilon);
  ionEnergy = (*(ionizationEnergy->find(Z)->second))[iosc];
  G4double Q2 = 
    photonEnergy0*photonEnergy0+photonEnergy1*(photonEnergy1-2.0*photonEnergy0*cosTheta);
  G4double cosThetaE; //scattering angle for the electron
  if (Q2 > 1.0e-12)
    {
      cosThetaE = (photonEnergy0-photonEnergy1*cosTheta)/std::sqrt(Q2);
    }
  else
    {
      cosThetaE = 1.0;
    }
  G4double sinThetaE = std::sqrt(1-cosThetaE*cosThetaE);

  //initialize here, then check photons created by Atomic-Deexcitation, and the final state e-
  std::vector<G4DynamicParticle*>* photonVector=0;

  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,iosc);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();
  G4double ionEnergyInPenelopeDatabase = ionEnergy;
  //protection against energy non-conservation
  ionEnergy = std::max(bindingEnergy,ionEnergyInPenelopeDatabase);  

  //subtract the excitation energy. If not emitted by fluorescence
  //the ionization energy is deposited as local energy deposition
  G4double eKineticEnergy = diffEnergy - ionEnergy; 
  G4double localEnergyDeposit = ionEnergy; 
  G4double energyInFluorescence = 0.; //testing purposes only

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
  if(DeexcitationFlag() && Z > 5) {

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
		    if (itsEnergy <= localEnergyDeposit)
		      {
			localEnergyDeposit -= itsEnergy;
			if (aPhoton->GetDefinition() == G4Gamma::Gamma()) 
			  energyInFluorescence += itsEnergy;;
			fvect->push_back(aPhoton);		    
		      }
		    else
		      {
			delete aPhoton;
			(*photonVector)[k]=0;
		      }
		  }
	      }
	    delete photonVector;
	  }
      }
  }

  //Produce explicitely the electron only if its proposed kinetic energy is 
  //above the cut, otherwise add local energy deposition
  G4DynamicParticle* electron = 0;
  //  if (eKineticEnergy > cutE) // VI: may be consistency problem if cut is applied here
  if (eKineticEnergy > 0.0)
    {
      G4double xEl = sinThetaE * std::cos(phi+pi); 
      G4double yEl = sinThetaE * std::sin(phi+pi);
      G4double zEl = cosThetaE;
      G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
      eDirection.rotateUz(photonDirection0);
      electron = new G4DynamicParticle (G4Electron::Electron(),
					eDirection,eKineticEnergy) ;
      fvect->push_back(electron);
    }
  else
    {
      localEnergyDeposit += eKineticEnergy;
    }

  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING-" 
	     << "G4PenelopeComptonModel::SampleSecondaries - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit=0.;
    }
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);
  
  G4double electronEnergy = 0.;
  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeCompton" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Scattered photon: " << photonEnergy1/keV << " keV" << G4endl;
      if (electron)
	electronEnergy = eKineticEnergy;
      G4cout << "Scattered electron " << electronEnergy/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (photonEnergy1+electronEnergy+energyInFluorescence+
					  localEnergyDeposit)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(photonEnergy1+
				      electronEnergy+energyInFluorescence+
				      localEnergyDeposit-photonEnergy0);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeCompton: problem with energy conservation: " << 
	  (photonEnergy1+electronEnergy+energyInFluorescence+localEnergyDeposit)/keV << 
	  " keV (final) vs. " << 
	  photonEnergy0/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::ReadData()
{
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeComptonModel - G4LEDATA environment variable not set!";
      G4Exception(excep);
      return;
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/compton-pen.dat";
  std::ifstream file(pathFile);
  
  if (!file.is_open())
    {
      G4String excep = "G4PenelopeComptonModel - data file " + pathFile + " not found!";
      G4Exception(excep);
    }

  G4int k1,test,test1;
  G4double a1,a2;
  G4int Z=1,nLevels=0;

  if (!ionizationEnergy || !hartreeFunction || !occupationNumber)
    {
      G4String excep = "G4PenelopeComptonModel: problem with reading data from file";
      G4Exception(excep);
      return;
    }

  do{
    G4double harOfElectronsBelowThreshold = 0;
    G4int nbOfElectronsBelowThreshold = 0; 
    file >> Z >> nLevels;
    //Check for nLevels validity, before using it in a loop
    if (nLevels<0 || nLevels>64)
      {
	G4String excep = "G4PenelopeComptonModel: corrupted data file?";
	G4Exception(excep);
	return;
      }
    G4DataVector* occVector = new G4DataVector;
    G4DataVector* harVector = new G4DataVector;
    G4DataVector* bindingEVector = new G4DataVector;
    for (G4int h=0;h<nLevels;h++)
      {
	file >> k1 >> a1 >> a2;
	//Make explicit unit of measurements for ionisation energy, which is MeV
        a1 *= MeV; 
	if (a1 > 15*eV)
	  {
	    occVector->push_back((G4double) k1);
	    bindingEVector->push_back(a1);
	    harVector->push_back(a2);
	  }
	else
	  {
	    nbOfElectronsBelowThreshold += k1;
	    harOfElectronsBelowThreshold += k1*a2;
	  }
      }
    //Add the "final" level
    if (nbOfElectronsBelowThreshold)
      {
	occVector->push_back(nbOfElectronsBelowThreshold);
	bindingEVector->push_back(0*eV);
	G4double averageHartree = 
	  harOfElectronsBelowThreshold/((G4double) nbOfElectronsBelowThreshold);
	harVector->push_back(averageHartree);
      }
    //Ok, done for element Z
    occupationNumber->insert(std::make_pair(Z,occVector));
    ionizationEnergy->insert(std::make_pair(Z,bindingEVector));
    hartreeFunction->insert(std::make_pair(Z,harVector));
    file >> test >> test1; //-1 -1 close the data for each Z
    if (test > 0) {
      G4String excep = "G4PenelopeComptonModel - data file corrupted!";
      G4Exception(excep);
    }
  }while (test != -2); //the very last Z is closed with -2 instead of -1
  file.close();
  if (verboseLevel > 2)
    {
      G4cout << "Data from G4PenelopeComptonModel read " << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeComptonModel::DifferentialCrossSection(G4double cosTheta)
{
  //
  // Penelope model for the Compton scattering differential cross section 
  // dSigma/dOmega.
  // D. Brusa et al., Nucl. Instrum. Meth. A 379 (1996) 167
  // The parametrization includes the J(p) distribution profiles for the atomic shells, 
  // that are tabulated from Hartree-Fock calculations 
  // from F. Biggs et al., At. Data Nucl. Data Tables 16 (1975) 201
  //
  const G4double k2 = std::sqrt(2.0);
  const G4double k1 = std::sqrt(0.5);
  const G4double k12 = 0.5;
  G4double cdt1 = 1.0-cosTheta;
  G4double energy = energyForIntegration;
  G4int Z = ZForIntegration;
  //energy of Compton line;
  G4double EOEC = 1.0+(energy/electron_mass_c2)*cdt1; 
  G4double ECOE = 1.0/EOEC;
  //Incoherent scattering function (analytical profile)
  G4double sia = 0.0;
  G4int nosc = occupationNumber->find(Z)->second->size();
  for (G4int i=0;i<nosc;i++){
    G4double ionEnergy = (*(ionizationEnergy->find(Z)->second))[i];
    //Sum only of those shells for which E>Eion
    if (energy > ionEnergy)
      {
        G4double aux = energy * (energy-ionEnergy)*cdt1;
	G4double Pzimax = 
	  (aux - electron_mass_c2*ionEnergy)/(electron_mass_c2*std::sqrt(2*aux+ionEnergy*ionEnergy));
	G4double harFunc = (*(hartreeFunction->find(Z)->second))[i]/fine_structure_const;
	G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
	G4double x = harFunc*Pzimax;
	G4double siap = 0;
	if (x > 0) 
	  {
	    siap = 1.0-0.5*std::exp(k12-(k1+k2*x)*(k1+k2*x));
	  }
	else
	  {
	    siap = 0.5*std::exp(k12-(k1-k2*x)*(k1-k2*x));
	  }
	sia = sia + occupNb*siap; //sum of all contributions;
      }
  }
  G4double XKN = EOEC+ECOE-1+cosTheta*cosTheta;
  G4double diffCS = pi*classic_electr_radius*classic_electr_radius*ECOE*ECOE*XKN*sia;
  return diffCS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeComptonModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4PenelopeComptonModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
