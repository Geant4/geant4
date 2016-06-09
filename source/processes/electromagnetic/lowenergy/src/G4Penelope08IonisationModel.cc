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
// $Id: G4Penelope08IonisationModel.cc,v 1.3 2010-12-15 10:26:41 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 27 Jul 2010   L Pandola    First complete implementation
//

#include "G4Penelope08IonisationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeOscillator.hh"
#include "G4PenelopeCrossSection.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
 
G4Penelope08IonisationModel::G4Penelope08IonisationModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),
   kineticEnergy1(0.*eV),
   cosThetaPrimary(1.0),
   energySecondary(0.*eV),
   cosThetaSecondary(0.0),
   targetOscillator(-1),XSTableElectron(0),XSTablePositron(0),
   theDeltaTable(0),energyGrid(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  nBins = 200;
  //
  oscManager = G4PenelopeOscillatorManager::GetOscillatorManager();
  //
  verboseLevel= 0;
   
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  // Atomic deexcitation model activated by default
  SetDeexcitationFlag(true);
  ActivateAuger(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4Penelope08IonisationModel::~G4Penelope08IonisationModel()
{
  ClearTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4Penelope08IonisationModel::Initialise(const G4ParticleDefinition*,
					     const G4DataVector&)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4Penelope08IonisationModel::Initialise()" << G4endl;
 
  //Clear and re-build the tables
  ClearTables();
  XSTableElectron = new 
    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;
  XSTablePositron = new 
    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;

  theDeltaTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;

  //Set the number of bins for the tables. 20 points per decade
  nBins = (size_t) (20*std::log10(HighEnergyLimit()/LowEnergyLimit()));
  nBins = std::max(nBins,(size_t)100);

  energyGrid = new G4PhysicsLogVector(LowEnergyLimit(),
				      HighEnergyLimit(), 
				      nBins-1); //one hidden bin is added

  if (verboseLevel > 2) {
    G4cout << "Penelope Ionisation model is initialized " << G4endl
	   << "Energy range: "
           << LowEnergyLimit() / keV << " keV - "
           << HighEnergyLimit() / GeV << " GeV. Using " 
	   << nBins << " bins." 
           << G4endl;
  }
 
  if(isInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4Penelope08IonisationModel::CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* theParticle,
                                           G4double energy,
                                           G4double cutEnergy,
                                           G4double)
{
  // Penelope model to calculate the cross section for inelastic collisions above the
  // threshold. It makes use of the Generalised Oscillator Strength (GOS) model from
  //  D. Liljequist, J. Phys. D: Appl. Phys. 16 (1983) 1567
  //
  // The total cross section is calculated analytically by taking
  // into account the atomic oscillators coming into the play for a given threshold.
  //
  // For incident e- the maximum energy allowed for the delta-rays is energy/2.
  // because particles are undistinghishable.
  //
  // The contribution is splitted in three parts: distant longitudinal collisions,
  // distant transverse collisions and close collisions. Each term is described by
  // its own analytical function.
  // Fermi density correction is calculated analytically according to
  //  U. Fano, Ann. Rev. Nucl. Sci. 13 (1963),1
  //
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4Penelope08IonisationModel" << G4endl;
 
  SetupForMaterial(theParticle, material, energy);
   
  G4double totalCross = 0.0;
  G4double crossPerMolecule = 0.;

  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(theParticle,material,
								cutEnergy);

  if (theXS)
    crossPerMolecule = theXS->GetHardCrossSection(energy);

  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double atPerMol =  oscManager->GetAtomsPerMolecule(material);
 
  if (verboseLevel > 3)
    G4cout << "Material " << material->GetName() << " has " << atPerMol <<
      "atoms per molecule" << G4endl;

  G4double moleculeDensity = 0.; 
  if (atPerMol)
    moleculeDensity = atomDensity/atPerMol;
 
  G4double crossPerVolume = crossPerMolecule*moleculeDensity;

  if (verboseLevel > 2)
  {
    G4cout << "G4Penelope08IonisationModel " << G4endl;
    G4cout << "Mean free path for delta emission > " << cutEnergy/keV << " keV at " <<
      energy/keV << " keV = " << (1./crossPerVolume)/mm << " mm" << G4endl;
    if (theXS)
      totalCross = (theXS->GetTotalCrossSection(energy))*moleculeDensity;
    G4cout << "Total free path for ionisation (no threshold) at " <<
      energy/keV << " keV = " << (1./totalCross)/mm << " mm" << G4endl;
  }
  return crossPerVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
                                                                                                     
//This is a dummy method. Never inkoved by the tracking, it just issues
//a warning if one tries to get Cross Sections per Atom via the
//G4EmCalculator.
G4double G4Penelope08IonisationModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                             G4double,
                                                             G4double,
                                                             G4double,
                                                             G4double,
                                                             G4double)
{
  G4cout << "*** G4Penelope08IonisationModel -- WARNING ***" << G4endl;
  G4cout << "Penelope08 Ionisation model does not calculate cross section _per atom_ " << G4endl;
  G4cout << "so the result is always zero. For physics values, please invoke " << G4endl;
  G4cout << "GetCrossSectionPerVolume() or GetMeanFreePath() via the G4EmCalculator" << G4endl;
  return 0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08IonisationModel::ComputeDEDXPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* theParticle,
                                           G4double kineticEnergy,
                                           G4double cutEnergy)
{
  // Penelope model to calculate the stopping power for soft inelastic collisions
  // below the threshold. It makes use of the Generalised Oscillator Strength (GOS)
  // model from
  //  D. Liljequist, J. Phys. D: Appl. Phys. 16 (1983) 1567
  //
  // The stopping power is calculated analytically using the dsigma/dW cross
  // section from the GOS models, which includes separate contributions from
  // distant longitudinal collisions, distant transverse collisions and
  // close collisions. Only the atomic oscillators that come in the play
  // (according to the threshold) are considered for the calculation.
  // Differential cross sections have a different form for e+ and e-.
  //
  // Fermi density correction is calculated analytically according to
  //  U. Fano, Ann. Rev. Nucl. Sci. 13 (1963),1
 
  if (verboseLevel > 3)
    G4cout << "Calling ComputeDEDX() of G4Penelope08IonisationModel" << G4endl;
 

  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(theParticle,material,
								cutEnergy);
  G4double sPowerPerMolecule = 0.0;
  if (theXS)
    sPowerPerMolecule = theXS->GetSoftStoppingPower(kineticEnergy);

  
  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double atPerMol =  oscManager->GetAtomsPerMolecule(material);
                                                                                        
  G4double moleculeDensity = 0.; 
  if (atPerMol)
    moleculeDensity = atomDensity/atPerMol;
 
  G4double sPowerPerVolume = sPowerPerMolecule*moleculeDensity;

  if (verboseLevel > 2)
    {
      G4cout << "G4Penelope08IonisationModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " <<
        kineticEnergy/keV << " keV = " << 
	sPowerPerVolume/(keV/mm) << " keV/mm" << G4endl;
    }
  return sPowerPerVolume;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08IonisationModel::MinEnergyCut(const G4ParticleDefinition*,
						 const G4MaterialCutsCouple*)
{
  return fIntrinsicLowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08IonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                              const G4MaterialCutsCouple* couple,
                                              const G4DynamicParticle* aDynamicParticle,
                                              G4double cutE, G4double)
{
  // Penelope model to sample the final state following an hard inelastic interaction.
  // It makes use of the Generalised Oscillator Strength (GOS) model from
  //  D. Liljequist, J. Phys. D: Appl. Phys. 16 (1983) 1567
  //
  // The GOS model is used to calculate the individual cross sections for all
  // the atomic oscillators coming in the play, taking into account the three
  // contributions (distant longitudinal collisions, distant transverse collisions and
  // close collisions). Then the target shell and the interaction channel are
  // sampled. Final state of the delta-ray (energy, angle) are generated according
  // to the analytical distributions (dSigma/dW) for the selected interaction
  // channels.
  // For e-, the maximum energy for the delta-ray is initialEnergy/2. (because
  // particles are indistinghusbable), while it is the full initialEnergy for
  // e+.
  // The efficiency of the random sampling algorithm (e.g. for close collisions)
  // decreases when initial and cutoff energy increase (e.g. 87% for 10-keV primary
  // and 1 keV threshold, 99% for 10-MeV primary and 10-keV threshold).
  // Differential cross sections have a different form for e+ and e-.
  //
  // WARNING: The model provides an _average_ description of inelastic collisions.
  // Anyway, the energy spectrum associated to distant excitations of a given
  // atomic shell is approximated as a single resonance. The simulated energy spectra
  // show _unphysical_ narrow peaks at energies that are multiple of the shell
  // resonance energies. The spurious speaks are automatically smoothed out after
  // multiple inelastic collisions.
  //
  // The model determines also the original shell from which the delta-ray is expelled,
  // in order to produce fluorescence de-excitation (from G4DeexcitationManager)
  //
  // Fermi density correction is calculated analytically according to
  //  U. Fano, Ann. Rev. Nucl. Sci. 13 (1963),1

  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4Penelope08IonisationModel" << G4endl;
 
  G4double kineticEnergy0 = aDynamicParticle->GetKineticEnergy();
  const G4ParticleDefinition* theParticle = aDynamicParticle->GetDefinition();
 
  if (kineticEnergy0 <= fIntrinsicLowEnergyLimit)
  {
    fParticleChange->SetProposedKineticEnergy(0.);
    fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy0);
    return ;
  }
  const G4Material* material = couple->GetMaterial();
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(material); 

  G4ParticleMomentum particleDirection0 = aDynamicParticle->GetMomentumDirection();
 
  //Initialise final-state variables. The proper values will be set by the methods
  // SampleFinalStateElectron() and SampleFinalStatePositron()
  kineticEnergy1=kineticEnergy0;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;
  targetOscillator = -1;
     
  if (theParticle == G4Electron::Electron())
    SampleFinalStateElectron(material,cutE,kineticEnergy0);
  else if (theParticle == G4Positron::Positron())
    SampleFinalStatePositron(material,cutE,kineticEnergy0);
  else
    {
      G4cout << "G4Penelope08IonisationModel::SamplingSecondaries() - " ;
      G4cout << "Invalid particle " << theParticle->GetParticleName() << G4endl;
      G4Exception();	
    }
 
  if (verboseLevel > 3)
  {
     G4cout << "G4Penelope08IonisationModel::SamplingSecondaries() for " << 
	theParticle->GetParticleName() << G4endl;
      G4cout << "Final eKin = " << kineticEnergy1 << " keV" << G4endl;
      G4cout << "Final cosTheta = " << cosThetaPrimary << G4endl;
      G4cout << "Delta-ray eKin = " << energySecondary << " keV" << G4endl;
      G4cout << "Delta-ray cosTheta = " << cosThetaSecondary << G4endl;
      G4cout << "Oscillator: " << targetOscillator << G4endl;
   }
   
  //Update the primary particle
  G4double sint = std::sqrt(1. - cosThetaPrimary*cosThetaPrimary);
  G4double phiPrimary  = twopi * G4UniformRand();
  G4double dirx = sint * std::cos(phiPrimary);
  G4double diry = sint * std::sin(phiPrimary);
  G4double dirz = cosThetaPrimary;
 
  G4ThreeVector electronDirection1(dirx,diry,dirz);
  electronDirection1.rotateUz(particleDirection0);
   
  if (kineticEnergy1 > 0)
    {
      fParticleChange->ProposeMomentumDirection(electronDirection1);
      fParticleChange->SetProposedKineticEnergy(kineticEnergy1);
    }
  else
    {
      fParticleChange->SetProposedKineticEnergy(0.);
    }
  
  //Generate the delta ray
  G4double ionEnergyInPenelopeDatabase = 
    (*theTable)[targetOscillator]->GetIonisationEnergy();
  
  //Now, try to handle fluorescence
  //Notice: merged levels are indicated with Z=0 and flag=30
  G4int shFlag = (*theTable)[targetOscillator]->GetShellFlag(); 
  G4int Z = (G4int) (*theTable)[targetOscillator]->GetParentZ();

  //initialize here, then check photons created by Atomic-Deexcitation, and the final state e-
  std::vector<G4DynamicParticle*>* photonVector=0;  
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4double bindingEnergy = 0.*eV;
  G4int shellId = 0;

  //Real level
  if (Z > 0 && shFlag<30)
    {
      const G4AtomicShell* shell = transitionManager->Shell(Z,shFlag-1);
      bindingEnergy = shell->BindingEnergy();
      shellId = shell->ShellId();
    }

  //correct the energySecondary to account for the fact that the Penelope 
  //database of ionisation energies is in general (slightly) different 
  //from the fluorescence database used in Geant4.
  energySecondary += ionEnergyInPenelopeDatabase-bindingEnergy;
  
  G4double localEnergyDeposit = bindingEnergy;
  //testing purposes only
  G4double energyInFluorescence = 0;

  if (energySecondary < 0)
    {
      //It means that there was some problem/mismatch between the two databases. 
      //In this case, the available energy is ok to excite the level according 
      //to the Penelope database, but not according to the Geant4 database
      //Full residual energy is deposited locally
      localEnergyDeposit += energySecondary;
      energySecondary = 0.0;
    }
  
  //the local energy deposit is what remains: part of this may be spent for fluorescence.
  if(DeexcitationFlag() && Z > 5) 
    {
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

  //Always produce explicitely the delta electron 
  G4DynamicParticle* electron = 0;
  G4double sinThetaE = std::sqrt(1.-cosThetaSecondary*cosThetaSecondary);
  G4double phiEl = phiPrimary+pi; //pi with respect to the primary electron/positron
  G4double xEl = sinThetaE * std::cos(phiEl);
  G4double yEl = sinThetaE * std::sin(phiEl);
  G4double zEl = cosThetaSecondary;
  G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
  eDirection.rotateUz(particleDirection0);
  electron = new G4DynamicParticle (G4Electron::Electron(),
				    eDirection,energySecondary) ;
  fvect->push_back(electron);


  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING-" 
	     << "G4Penelope08IonisationModel::SampleSecondaries - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit=0.;
    }
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4Penelope08Ionisation" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << kineticEnergy1/keV << " keV" << G4endl;
      G4cout << "Delta ray " << energySecondary/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (energySecondary+energyInFluorescence+kineticEnergy1+
					      localEnergyDeposit)/keV <<
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
      
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(energySecondary+energyInFluorescence+kineticEnergy1+
				      localEnergyDeposit-kineticEnergy0);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeIonisation: problem with energy conservation: " <<  
	  (energySecondary+energyInFluorescence+kineticEnergy1+localEnergyDeposit)/keV <<
	  " keV (final) vs. " <<
	  kineticEnergy0/keV << " keV (initial)" << G4endl;
    }
    
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4Penelope08IonisationModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4Penelope08IonisationModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08IonisationModel::ClearTables()
{  
  std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>::iterator i;
  if (XSTableElectron)
    {
      for (i=XSTableElectron->begin(); i != XSTableElectron->end(); i++)
	{
	  G4PenelopeCrossSection* tab = i->second;
	  delete tab;
	}
      delete XSTableElectron;
      XSTableElectron = 0;
    }

  if (XSTablePositron)
    {
      for (i=XSTablePositron->begin(); i != XSTablePositron->end(); i++)
	{
	  G4PenelopeCrossSection* tab = i->second;
	  delete tab;
	}
      delete XSTablePositron;
      XSTablePositron = 0;
    }

  std::map<const G4Material*,G4PhysicsFreeVector*>::iterator k;
  if (theDeltaTable)
    {
      for (k=theDeltaTable->begin();k!=theDeltaTable->end();k++)	
	delete k->second;
      delete theDeltaTable;
      theDeltaTable = 0;
    }
    
  if (energyGrid)
    delete energyGrid;

  if (verboseLevel > 2)
    G4cout << "G4Penelope08IonisationModel: cleared tables" << G4endl;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeCrossSection* 
G4Penelope08IonisationModel::GetCrossSectionTableForCouple(const G4ParticleDefinition* part,
							   const G4Material* mat,
							   G4double cut)
{
  if (part != G4Electron::Electron() && part != G4Positron::Positron())
    {
      G4cout << "G4Penelope08IonisationModel::GetCrossSectionTableForCouple" << G4endl;
      G4cout << "Invalid particle: " << part->GetParticleName() << G4endl;
      G4Exception();
      return NULL;
    }

  if (part == G4Electron::Electron())
    {
      if (!XSTableElectron)
	{
	  G4cout << "G4Penelope08IonisationModel::GetCrossSectionTableForCouple()" << G4endl;
	  G4cout << "The Cross Section Table for e- was not initialized correctly!" << G4endl;
	  G4Exception();	  
	  return NULL;
	}
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (XSTableElectron->count(theKey)) //table already built	
	return XSTableElectron->find(theKey)->second;
      else
	{
	  BuildXSTable(mat,cut,part);
	  if (XSTableElectron->count(theKey)) //now it should be ok!
	    return XSTableElectron->find(theKey)->second;
	  else
	    {
	      G4cout << "G4Penelope08IonisationModel::GetCrossSectionTableForCouple()" << G4endl;
	      G4cout << "Unable to build e- table for " << mat->GetName() << G4endl;
	      G4Exception();	   
	    }
	}
    }

  if (part == G4Positron::Positron())
    {
      if (!XSTablePositron)
	{
	  G4cout << "G4Penelope08IonisationModel::GetCrossSectionTableForCouple()" << G4endl;
	  G4cout << "The Cross Section Table for e+ was not initialized correctly!" << G4endl;
	  G4Exception();
	  return NULL;
	}
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (XSTablePositron->count(theKey)) //table already built	
	return XSTablePositron->find(theKey)->second;
      else
	{
	  BuildXSTable(mat,cut,part);
	  if (XSTablePositron->count(theKey)) //now it should be ok!
	    return XSTablePositron->find(theKey)->second;
	  else
	    {
	      G4cout << "G4Penelope08IonisationModel::GetCrossSectionTableForCouple()" << G4endl;
	      G4cout << "Unable to build e+ table for " << mat->GetName() << G4endl;
	      G4Exception();	   
	    }
	}
    }
  return NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08IonisationModel::BuildXSTable(const G4Material* mat,G4double cut,
					       const G4ParticleDefinition* part)
{
  //
  //This method fills the G4PenelopeCrossSection containers for electrons or positrons
  //and for the given material/cut couple. The calculation is done as sum over the 
  //individual shells.
  //Equivalent of subroutines EINaT and PINaT of Penelope
  //
  if (verboseLevel > 2)
    {
      G4cout << "G4Penelope08IonisationModel: going to build cross section table " << G4endl;
      G4cout << "for " << part->GetParticleName() << " in " << mat->GetName() << G4endl;
    }

  //Tables have been already created (checked by GetCrossSectionTableForCouple)
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  size_t numberOfOscillators = theTable->size();

  if (energyGrid->GetVectorLength() != nBins) 
    {
      G4cout << "G4Penelope08IonisationModel::BuildXSTable" << G4endl;
      G4cout << "Energy Grid looks not initialized" << G4endl;
      G4cout << nBins << " " << energyGrid->GetVectorLength() << G4endl;
      G4Exception();
    }

  G4PenelopeCrossSection* XSEntry = new G4PenelopeCrossSection(nBins,numberOfOscillators);
 
  //loop on the energy grid
  for (size_t bin=0;bin<nBins;bin++)
    {
       G4double energy = energyGrid->GetLowEdgeEnergy(bin);
       G4double XH0=0, XH1=0, XH2=0;
       G4double XS0=0, XS1=0, XS2=0;
   
       //oscillator loop
       for (size_t iosc=0;iosc<numberOfOscillators;iosc++)
	 {
	   G4DataVector* tempStorage = 0;

	   G4PenelopeOscillator* theOsc = (*theTable)[iosc];
	   G4double delta = GetDensityCorrection(mat,energy);
	   if (part == G4Electron::Electron())	     
	     tempStorage = ComputeShellCrossSectionsElectron(theOsc,energy,cut,delta);
	   else if (part == G4Positron::Positron())
	     tempStorage = ComputeShellCrossSectionsPositron(theOsc,energy,cut,delta);
	   //check results are all right
	   if (!tempStorage)
	     {
	       G4cout << "G4Penelope08IonisationModel::BuildXSTable" << G4endl;
	       G4cout << "Problem in calculating the shell XS " << G4endl;
	       G4Exception();
	       return;
	     }
	   if (tempStorage->size() != 6)
	     {
	       G4cout << "G4Penelope08IonisationModel::BuildXSTable" << G4endl;
	       G4cout << "Problem in calculating the shell XS " << G4endl;
	       G4cout << "Result has dimension " << tempStorage->size() << " instead of 6" << G4endl;
	       G4Exception();
	     }
	   G4double stre = theOsc->GetOscillatorStrength();

	   XH0 += stre*(*tempStorage)[0];
	   XH1 += stre*(*tempStorage)[1];
	   XH2 += stre*(*tempStorage)[2];
	   XS0 += stre*(*tempStorage)[3];
	   XS1 += stre*(*tempStorage)[4];
	   XS2 += stre*(*tempStorage)[5];
	   XSEntry->AddShellCrossSectionPoint(bin,iosc,energy,stre*(*tempStorage)[0]);
	   if (tempStorage)
	     {
	       delete tempStorage;
	       tempStorage = 0;
	     }
	 }       
       XSEntry->AddCrossSectionPoint(bin,energy,XH0,XH1,XH2,XS0,XS1,XS2);
    }


  //Normalize shell cross sections before insertion in the table
  XSEntry->NormalizeShellCrossSections();


  //Insert in the appropriate table
  std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
  if (part == G4Electron::Electron())      
    XSTableElectron->insert(std::make_pair(theKey,XSEntry));
  else if (part == G4Positron::Positron())
    XSTablePositron->insert(std::make_pair(theKey,XSEntry));
  else
    delete XSEntry;
  
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08IonisationModel::GetDensityCorrection(const G4Material* mat,
							   G4double energy)
{
  G4double result = 0;
  if (!theDeltaTable)
    {
      G4cout << "G4Penelope08IonisationModel::GetDensityCorrection()" << G4endl;
      G4cout << "Delta Table not initialized. Was Initialise() run?" << G4endl;
      G4Exception();
      return 0;
    }
  if (energy <= 0*eV)
    {
      G4cout << "G4Penelope08IonisationModel::GetDensityCorrection()" << G4endl;
      G4cout << "Invalid energy " << energy/eV << " eV " << G4endl;
      return 0;
    }
  G4double logene = std::log(energy);

  //check if the material has been built
  if (!(theDeltaTable->count(mat)))
    BuildDeltaTable(mat);

  if (theDeltaTable->count(mat))
    {
      G4PhysicsFreeVector* vec = theDeltaTable->find(mat)->second;
      result = vec->Value(logene); //the table has delta vs. ln(E)      
    }
  else
    {
      G4cout << "G4Penelope08IonisationModel::GetDensityCorrection()" << G4endl;
      G4cout << "Unable to build table for " << mat->GetName() << G4endl;
      G4Exception();
    }

  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08IonisationModel::BuildDeltaTable(const G4Material* mat)
{
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  G4double plasmaSq = oscManager->GetPlasmaEnergySquared(mat);
  G4double totalZ = oscManager->GetTotalZ(mat);
  size_t numberOfOscillators = theTable->size();

  if (energyGrid->GetVectorLength() != nBins) 
    {
      G4cout << "G4Penelope08IonisationModel::BuildDeltaTable" << G4endl;
      G4cout << "Energy Grid for Delta looks not initialized" << G4endl;
      G4cout << nBins << " " << energyGrid->GetVectorLength() << G4endl;
      G4Exception();
    }

  G4PhysicsFreeVector* theVector = new G4PhysicsFreeVector(nBins);

  //loop on the energy grid
  for (size_t bin=0;bin<nBins;bin++)
    {
      G4double delta = 0.;
      G4double energy = energyGrid->GetLowEdgeEnergy(bin);

      //Here calculate delta
      G4double gam = 1.0+(energy/electron_mass_c2);
      G4double gamSq = gam*gam;

      G4double TST = totalZ/(gamSq*plasmaSq);
      G4double wl2 = 0;
      G4double fdel = 0;

      //loop on oscillators
      for (size_t i=0;i<numberOfOscillators;i++)
	{
	  G4PenelopeOscillator* theOsc = (*theTable)[i];
	  G4double wri = theOsc->GetResonanceEnergy();
	  fdel += theOsc->GetOscillatorStrength()/(wri*wri+wl2);
	}      
      if (fdel >= TST) //if fdel < TST, delta = 0
	{
	  //get last oscillator
	  G4PenelopeOscillator* theOsc = (*theTable)[numberOfOscillators-1]; 
	  wl2 = theOsc->GetResonanceEnergy()*theOsc->GetResonanceEnergy();
	  
	  //First iteration
	  G4bool loopAgain = false;
	  do
	    {
	      loopAgain = false;
	      wl2 += wl2;
	      fdel = 0.;
	      for (size_t i=0;i<numberOfOscillators;i++)
		{
		  G4PenelopeOscillator* theOsc = (*theTable)[i];
		  G4double wri = theOsc->GetResonanceEnergy();
		  fdel += theOsc->GetOscillatorStrength()/(wri*wri+wl2);
		}
	      if (fdel > TST)
		loopAgain = true;
	    }while(loopAgain);

	  G4double wl2l = 0;
	  G4double wl2u = wl2;
	  //second iteration
	  do
	    {	     
	      loopAgain = false;
	      wl2 = 0.5*(wl2l+wl2u);
	      fdel = 0;
	      for (size_t i=0;i<numberOfOscillators;i++)
		{
		  G4PenelopeOscillator* theOsc = (*theTable)[i];
		  G4double wri = theOsc->GetResonanceEnergy();
		  fdel += theOsc->GetOscillatorStrength()/(wri*wri+wl2);
		}
	      if (fdel > TST)
		wl2l = wl2;
	      else
		wl2u = wl2;
	      if ((wl2u-wl2l)>1e-12*wl2)
		loopAgain = true;
	    }while(loopAgain);
	  
	  //Eventually get density correction
	  delta = 0.;
	  for (size_t i=0;i<numberOfOscillators;i++)
	    {
	      G4PenelopeOscillator* theOsc = (*theTable)[i];
	      G4double wri = theOsc->GetResonanceEnergy();
	      delta += theOsc->GetOscillatorStrength()*
		std::log(1.0+(wl2/(wri*wri)));	  	     
	    }
	  delta = (delta/totalZ)-wl2/(gamSq*plasmaSq);
	}
      energy = std::max(1e-9*eV,energy); //prevents log(0)
      theVector->PutValue(bin,std::log(energy),delta);
    }
  theDeltaTable->insert(std::make_pair(mat,theVector));
  return;
}
							
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4DataVector* G4Penelope08IonisationModel::ComputeShellCrossSectionsElectron(G4PenelopeOscillator* theOsc,
									     G4double energy,
									     G4double cut,
									     G4double delta)
{
  //
  //This method calculates the hard and soft cross sections (H0-H1-H2-S0-S1-S2) for 
  //the given oscillator/cut and at the given energy.
  //It returns a G4DataVector* with 6 entries (H0-H1-H2-S0-S1-S2)
  //Equivalent of subroutines EINaT1 of Penelope
  //
  // Results are _per target electron_
  //
  G4DataVector* result = new G4DataVector();
  for (size_t i=0;i<6;i++)
    result->push_back(0.);
  G4double ionEnergy = theOsc->GetIonisationEnergy();
  
  //return a set of zero's if the energy it too low to excite the current oscillator
  if (energy < ionEnergy)
    return result;

  G4double H0=0.,H1=0.,H2=0.;
  G4double S0=0.,S1=0.,S2=0.;

  //Define useful constants to be used in the calculation
  G4double gamma = 1.0+energy/electron_mass_c2;
  G4double gammaSq = gamma*gamma;
  G4double beta = (gammaSq-1.0)/gammaSq;
  G4double pielr2 = pi*classic_electr_radius*classic_electr_radius; //pi*re^2
  G4double constant = pielr2*2.0*electron_mass_c2/beta;
  G4double XHDT0 = std::log(gammaSq)-beta;

  G4double cpSq = energy*(energy+2.0*electron_mass_c2);
  G4double cp = std::sqrt(cpSq);
  G4double amol = (energy/(energy+electron_mass_c2))*(energy/(energy+electron_mass_c2));

  //
  // Distant interactions
  //
  G4double resEne = theOsc->GetResonanceEnergy();
  G4double cutoffEne = theOsc->GetCutoffRecoilResonantEnergy();
  if (energy > resEne)
    {
      G4double cp1Sq = (energy-resEne)*(energy-resEne+2.0*electron_mass_c2);
      G4double cp1 = std::sqrt(cp1Sq);
      
      //Distant longitudinal interactions
      G4double QM = 0;
      if (resEne > 1e-6*energy)
	QM = std::sqrt((cp-cp1)*(cp-cp1)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
      else
	{
	  QM = resEne*resEne/(beta*2.0*electron_mass_c2);
	  QM = QM*(1.0-0.5*QM/electron_mass_c2);
	}
      G4double SDL1 = 0;
      if (QM < cutoffEne)
	SDL1 = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)));
      
      //Distant transverse interactions
      if (SDL1)
	{
	  G4double SDT1 = std::max(XHDT0-delta,0.0);
	  G4double SD1 = SDL1+SDT1;
	  if (cut > resEne)
	    {
	      S1 = SD1; //XS1
	      S0 = SD1/resEne; //XS0
	      S2 = SD1*resEne; //XS2
	    }
	  else
	    {
	      H1 = SD1; //XH1
	      H0 = SD1/resEne; //XH0
	      H2 = SD1*resEne; //XH2
	    }
	}
    }
  //
  // Close collisions (Moller's cross section)
  //
  G4double wl = std::max(cut,cutoffEne);
  G4double ee = energy + ionEnergy;
  G4double wu = 0.5*ee;
  if (wl < wu-(1e-5*eV))
    {
      H0 += (1.0/(ee-wu)) - (1.0/(ee-wl)) - (1.0/wu) + (1.0/wl) + 
	(1.0-amol)*std::log(((ee-wu)*wl)/((ee-wl)*wu))/ee + 
	amol*(wu-wl)/(ee*ee);
      H1 += std::log(wu/wl)+(ee/(ee-wu))-(ee/(ee-wl)) + 
	(2.0-amol)*std::log((ee-wu)/(ee-wl)) + 
	amol*(wu*wu-wl*wl)/(2.0*ee*ee);
      H2 += (2.0-amol)*(wu-wl)+(wu*(2.0*ee-wu)/(ee-wu)) - 
	(wl*(2.0*ee-wl)/(ee-wl)) + 
	(3.0-amol)*ee*std::log((ee-wu)/(ee-wl)) + 	
	amol*(wu*wu*wu-wl*wl*wl)/(3.0*ee*ee);
      wu = wl;
    }
  wl = cutoffEne;
  
  if (wl > wu-(1e-5*eV))
    {
      (*result)[0] = constant*H0;
      (*result)[1] = constant*H1;
      (*result)[2] = constant*H2;
      (*result)[3] = constant*S0;
      (*result)[4] = constant*S1;
      (*result)[5] = constant*S2;
      return result;
    }

  S0 += (1.0/(ee-wu))-(1.0/(ee-wl)) - (1.0/wu) + (1.0/wl) + 
    (1.0-amol)*std::log(((ee-wu)*wl)/((ee-wl)*wu))/ee +
    amol*(wu-wl)/(ee*ee);
  S1 += std::log(wu/wl)+(ee/(ee-wu))-(ee/(ee-wl)) + 
    (2.0-amol)*std::log((ee-wu)/(ee-wl)) + 
    amol*(wu*wu-wl*wl)/(2.0*ee*ee);
  S2 += (2.0-amol)*(wu-wl)+(wu*(2.0*ee-wu)/(ee-wu)) - 
    (wl*(2.0*ee-wl)/(ee-wl)) + 
    (3.0-amol)*ee*std::log((ee-wu)/(ee-wl)) + 
    amol*(wu*wu*wu-wl*wl*wl)/(3.0*ee*ee);

  (*result)[0] = constant*H0;
  (*result)[1] = constant*H1;
  (*result)[2] = constant*H2;
  (*result)[3] = constant*S0;
  (*result)[4] = constant*S1;
  (*result)[5] = constant*S2;
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4DataVector* G4Penelope08IonisationModel::ComputeShellCrossSectionsPositron(G4PenelopeOscillator* theOsc,
									     G4double energy,
									     G4double cut,
									     G4double delta)
{
  //
  //This method calculates the hard and soft cross sections (H0-H1-H2-S0-S1-S2) for 
  //the given oscillator/cut and at the given energy.
  //It returns a G4DataVector* with 6 entries (H0-H1-H2-S0-S1-S2)
  //Equivalent of subroutines PINaT1 of Penelope
  //
  // Results are _per target electron_
  //
  G4DataVector* result = new G4DataVector();
  for (size_t i=0;i<6;i++)
    result->push_back(0.);
  G4double ionEnergy = theOsc->GetIonisationEnergy();
  
  //return a set of zero's if the energy it too low to excite the current oscillator
  if (energy < ionEnergy)
    return result;

  G4double H0=0.,H1=0.,H2=0.;
  G4double S0=0.,S1=0.,S2=0.;

  //Define useful constants to be used in the calculation
  G4double gamma = 1.0+energy/electron_mass_c2;
  G4double gammaSq = gamma*gamma;
  G4double beta = (gammaSq-1.0)/gammaSq;
  G4double pielr2 = pi*classic_electr_radius*classic_electr_radius; //pi*re^2
  G4double constant = pielr2*2.0*electron_mass_c2/beta;
  G4double XHDT0 = std::log(gammaSq)-beta;

  G4double cpSq = energy*(energy+2.0*electron_mass_c2);
  G4double cp = std::sqrt(cpSq);
  G4double amol = (energy/(energy+electron_mass_c2))*(energy/(energy+electron_mass_c2));
  G4double g12 = (gamma+1.0)*(gamma+1.0);
  //Bhabha coefficients
  G4double bha1 = amol*(2.0*g12-1.0)/(gammaSq-1.0);
  G4double bha2 = amol*(3.0+1.0/g12);
  G4double bha3 = amol*2.0*gamma*(gamma-1.0)/g12;
  G4double bha4 = amol*(gamma-1.0)*(gamma-1.0)/g12;

  //
  // Distant interactions
  //
  G4double resEne = theOsc->GetResonanceEnergy();
  G4double cutoffEne = theOsc->GetCutoffRecoilResonantEnergy();
  if (energy > resEne)
    {
      G4double cp1Sq = (energy-resEne)*(energy-resEne+2.0*electron_mass_c2);
      G4double cp1 = std::sqrt(cp1Sq);
      
      //Distant longitudinal interactions
      G4double QM = 0;
      if (resEne > 1e-6*energy)
	QM = std::sqrt((cp-cp1)*(cp-cp1)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
      else
	{
	  QM = resEne*resEne/(beta*2.0*electron_mass_c2);
	  QM = QM*(1.0-0.5*QM/electron_mass_c2);
	}
      G4double SDL1 = 0;
      if (QM < cutoffEne)
	SDL1 = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)));
      
      //Distant transverse interactions
      if (SDL1)
	{
	  G4double SDT1 = std::max(XHDT0-delta,0.0);
	  G4double SD1 = SDL1+SDT1;
	  if (cut > resEne)
	    {
	      S1 = SD1; //XS1
	      S0 = SD1/resEne; //XS0
	      S2 = SD1*resEne; //XS2
	    }
	  else
	    {
	      H1 = SD1; //XH1
	      H0 = SD1/resEne; //XH0
	      H2 = SD1*resEne; //XH2
	    }
	}
    }

  //
  // Close collisions (Bhabha's cross section)
  //
  G4double wl = std::max(cut,cutoffEne);
  G4double wu = energy; 
  G4double energySq = energy*energy;
  if (wl < wu-(1e-5*eV))
    {
      G4double wlSq = wl*wl;
      G4double wuSq = wu*wu;
      H0 += (1.0/wl) - (1.0/wu)- bha1*std::log(wu/wl)/energy  
	+ bha2*(wu-wl)/energySq  
	- bha3*(wuSq-wlSq)/(2.0*energySq*energy)
	+ bha4*(wuSq*wu-wlSq*wl)/(3.0*energySq*energySq);
      H1 += std::log(wu/wl) - bha1*(wu-wl)/energy
	+ bha2*(wuSq-wlSq)/(2.0*energySq)
	- bha3*(wuSq*wu-wlSq*wl)/(3.0*energySq*energy)
	+ bha4*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energySq);
      H2 += wu - wl - bha1*(wuSq-wlSq)/(2.0*energy)
	+ bha2*(wuSq*wu-wlSq*wl)/(3.0*energySq)
	- bha3*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energy)
	+ bha4*(wuSq*wuSq*wu-wlSq*wlSq*wl)/(5.0*energySq*energySq);
      wu = wl;
    }
  wl = cutoffEne;
  
  if (wl > wu-(1e-5*eV))
    {
      (*result)[0] = constant*H0;
      (*result)[1] = constant*H1;
      (*result)[2] = constant*H2;
      (*result)[3] = constant*S0;
      (*result)[4] = constant*S1;
      (*result)[5] = constant*S2;
      return result;
    }

  G4double wlSq = wl*wl;
  G4double wuSq = wu*wu;

  S0 += (1.0/wl) - (1.0/wu) - bha1*std::log(wu/wl)/energy 
    + bha2*(wu-wl)/energySq  
    - bha3*(wuSq-wlSq)/(2.0*energySq*energy)
    + bha4*(wuSq*wu-wlSq*wl)/(3.0*energySq*energySq);

  S1 += std::log(wu/wl) - bha1*(wu-wl)/energy
    + bha2*(wuSq-wlSq)/(2.0*energySq)
    - bha3*(wuSq*wu-wlSq*wl)/(3.0*energySq*energy)
    + bha4*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energySq);

  S2 += wu - wl - bha1*(wuSq-wlSq)/(2.0*energy)
    + bha2*(wuSq*wu-wlSq*wl)/(3.0*energySq)
    - bha3*(wuSq*wuSq-wlSq*wlSq)/(4.0*energySq*energy)
    + bha4*(wuSq*wuSq*wu-wlSq*wlSq*wl)/(5.0*energySq*energySq);

 (*result)[0] = constant*H0;
 (*result)[1] = constant*H1;
 (*result)[2] = constant*H2;
 (*result)[3] = constant*S0;
 (*result)[4] = constant*S1;
 (*result)[5] = constant*S2;

 return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4Penelope08IonisationModel::SampleFinalStateElectron(const G4Material* mat,
							   G4double cutEnergy,
							   G4double kineticEnergy)
{
  // This method sets the final ionisation parameters
  // kineticEnergy1, cosThetaPrimary (= updates of the primary e-)
  // energySecondary, cosThetaSecondary (= info of delta-ray)
  // targetOscillator (= ionised oscillator)
  //
  // The method implements SUBROUTINE EINa of Penelope
  //

  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  size_t numberOfOscillators = theTable->size();
  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(G4Electron::Electron(),mat,
								cutEnergy);
  G4double delta = GetDensityCorrection(mat,kineticEnergy);
 
  // Selection of the active oscillator
  G4double TST = G4UniformRand();
  targetOscillator = numberOfOscillators-1; //initialization, last oscillator
  G4double XSsum = 0.;

  for (size_t i=0;i<numberOfOscillators-1;i++)
    {
      XSsum += theXS->GetShellCrossSection(i,kineticEnergy); 
     	
      if (XSsum > TST)
	{
	  targetOscillator = (G4int) i;
	  break;
	}
    }

  
  if (verboseLevel > 3)
    {
      G4cout << "SampleFinalStateElectron: sampled oscillator #" << targetOscillator << "." << G4endl;
      G4cout << "Ionisation energy: " << (*theTable)[targetOscillator]->GetIonisationEnergy()/eV << 
	" eV " << G4endl;
      G4cout << "Resonance energy: : " << (*theTable)[targetOscillator]->GetResonanceEnergy()/eV << " eV "
	     << G4endl;
    }
  //Constants
  G4double rb = kineticEnergy + 2.0*electron_mass_c2;
  G4double gam = 1.0+kineticEnergy/electron_mass_c2;
  G4double gam2 = gam*gam;
  G4double beta2 = (gam2-1.0)/gam2;
  G4double amol = ((gam-1.0)/gam)*((gam-1.0)/gam);

  //Partial cross section of the active oscillator
  G4double resEne = (*theTable)[targetOscillator]->GetResonanceEnergy();
  G4double invResEne = 1.0/resEne;
  G4double ionEne = (*theTable)[targetOscillator]->GetIonisationEnergy();
  G4double cutoffEne = (*theTable)[targetOscillator]->GetCutoffRecoilResonantEnergy();
  G4double XHDL = 0.;
  G4double XHDT = 0.;
  G4double QM = 0.;
  G4double cps = 0.;
  G4double cp = 0.;

  //Distant excitations
  if (resEne > cutEnergy && resEne < kineticEnergy)
    {
      cps = kineticEnergy*rb;
      cp = std::sqrt(cps);
      G4double XHDT0 = std::max(std::log(gam2)-beta2-delta,0.);
      if (resEne > 1.0e-6*kineticEnergy)
	{
	  G4double cpp = std::sqrt((kineticEnergy-resEne)*(kineticEnergy-resEne+2.0*electron_mass_c2));
	  QM = std::sqrt((cp-cpp)*(cp-cpp)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
	}
      else
	{
	  QM = resEne*resEne/(beta2*2.0*electron_mass_c2);
	  QM *= (1.0-QM*0.5/electron_mass_c2);
	}
      if (QM < cutoffEne)
	{
	  XHDL = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)))
	    *invResEne;
	  XHDT = XHDT0*invResEne;	  
	}
      else
	{
	  QM = cutoffEne;
	  XHDL = 0.;
	  XHDT = 0.;
	}
    }
  else
    {
      QM = cutoffEne;
      cps = 0.;
      cp = 0.;
      XHDL = 0.;
      XHDT = 0.;
    }

  //Close collisions
  G4double EE = kineticEnergy + ionEne;
  G4double wmaxc = 0.5*EE;
  G4double wcl = std::max(cutEnergy,cutoffEne);
  G4double rcl = wcl/EE;
  G4double XHC = 0.;
  if (wcl < wmaxc)
    {
      G4double rl1 = 1.0-rcl;
      G4double rrl1 = 1.0/rl1;
      XHC = (amol*(0.5-rcl)+1.0/rcl-rrl1+
	     (1.0-amol)*std::log(rcl*rrl1))/EE;
    }

  //Total cross section per molecule for the active shell, in cm2
  G4double XHTOT = XHC + XHDL + XHDT;

  //very small cross section, do nothing
  if (XHTOT < 1.e-14*barn)
    {
      kineticEnergy1=kineticEnergy;
      cosThetaPrimary=1.0;
      energySecondary=0.0;
      cosThetaSecondary=1.0;
      targetOscillator = numberOfOscillators-1;
      return;
    }
  
  //decide which kind of interaction we'll have
  TST = XHTOT*G4UniformRand();
  

  // Hard close collision
  G4double TS1 = XHC;
  
  if (TST < TS1)
    {
      G4double A = 5.0*amol;
      G4double ARCL = A*0.5*rcl;
      G4double rk=0.;
      G4bool loopAgain = false;
      do
	{
	  loopAgain = false;
	  G4double fb = (1.0+ARCL)*G4UniformRand();
	  if (fb < 1)
	    rk = rcl/(1.0-fb*(1.0-(rcl+rcl)));	     
	  else
	    rk = rcl + (fb-1.0)*(0.5-rcl)/ARCL;
	  G4double rk2 = rk*rk;
	  G4double rkf = rk/(1.0-rk);
	  G4double phi = 1.0+rkf*rkf-rkf+amol*(rk2+rkf);
	  if (G4UniformRand()*(1.0+A*rk2) > phi)
	    loopAgain = true;
	}while(loopAgain);
      //energy and scattering angle (primary electron)
      G4double deltaE = rk*EE;
      
      kineticEnergy1 = kineticEnergy - deltaE;
      cosThetaPrimary = std::sqrt(kineticEnergy1*rb/(kineticEnergy*(rb-deltaE)));
      //energy and scattering angle of the delta ray
      energySecondary = deltaE - ionEne; //subtract ionisation energy
      cosThetaSecondary= std::sqrt(deltaE*rb/(kineticEnergy*(deltaE+2.0*electron_mass_c2)));
      if (verboseLevel > 3)
	G4cout << "SampleFinalStateElectron: sampled close collision " << G4endl;
      return;     			     
    }

  //Hard distant longitudinal collisions
  TS1 += XHDL;
  G4double deltaE = resEne;
  kineticEnergy1 = kineticEnergy - deltaE;

  if (TST < TS1)
    {
      G4double QS = QM/(1.0+QM*0.5/electron_mass_c2);
      G4double Q = QS/(std::pow((QS/cutoffEne)*(1.0+cutoffEne*0.5/electron_mass_c2),G4UniformRand())
		       - (QS*0.5/electron_mass_c2));
      G4double QTREV = Q*(Q+2.0*electron_mass_c2);
      G4double cpps = kineticEnergy1*(kineticEnergy1+2.0*electron_mass_c2);
      cosThetaPrimary = (cpps+cps-QTREV)/(2.0*cp*std::sqrt(cpps));
      if (cosThetaPrimary > 1.) 
	cosThetaPrimary = 1.0;
      //energy and emission angle of the delta ray
      energySecondary = deltaE - ionEne;
      cosThetaSecondary = 0.5*(deltaE*(kineticEnergy+rb-deltaE)+QTREV)/std::sqrt(cps*QTREV);
      if (cosThetaSecondary > 1.0)
	cosThetaSecondary = 1.0;
      if (verboseLevel > 3)
	G4cout << "SampleFinalStateElectron: sampled distant longitudinal collision " << G4endl;
      return;      
    }

  //Hard distant transverse collisions
  cosThetaPrimary = 1.0;
  //energy and emission angle of the delta ray
  energySecondary = deltaE - ionEne;
  cosThetaSecondary = 0.5;
  if (verboseLevel > 3)
    G4cout << "SampleFinalStateElectron: sampled distant transverse collision " << G4endl;

  return;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4Penelope08IonisationModel::SampleFinalStatePositron(const G4Material* mat,
							   G4double cutEnergy,
							   G4double kineticEnergy)
{
  // This method sets the final ionisation parameters
  // kineticEnergy1, cosThetaPrimary (= updates of the primary e-)
  // energySecondary, cosThetaSecondary (= info of delta-ray)
  // targetOscillator (= ionised oscillator)
  //
  // The method implements SUBROUTINE PINa of Penelope
  // 
 
  G4PenelopeOscillatorTable* theTable = oscManager->GetOscillatorTableIonisation(mat);
  size_t numberOfOscillators = theTable->size();
  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(G4Positron::Positron(),mat,
								cutEnergy);
  G4double delta = GetDensityCorrection(mat,kineticEnergy);

  // Selection of the active oscillator
  G4double TST = G4UniformRand();
  targetOscillator = numberOfOscillators-1; //initialization, last oscillator
  G4double XSsum = 0.;
  for (size_t i=0;i<numberOfOscillators-1;i++)
    {
      XSsum += theXS->GetShellCrossSection(i,kineticEnergy);   	
      if (XSsum > TST)
	{
	  targetOscillator = (G4int) i;
	  break;
	}
    }

  if (verboseLevel > 3)
    {
      G4cout << "SampleFinalStatePositron: sampled oscillator #" << targetOscillator << "." << G4endl;
      G4cout << "Ionisation energy: " << (*theTable)[targetOscillator]->GetIonisationEnergy()/eV 
	     << " eV " << G4endl; 
      G4cout << "Resonance energy: : " << (*theTable)[targetOscillator]->GetResonanceEnergy()/eV 
	     << " eV " << G4endl;
    }

  //Constants
  G4double rb = kineticEnergy + 2.0*electron_mass_c2;
  G4double gam = 1.0+kineticEnergy/electron_mass_c2;
  G4double gam2 = gam*gam;
  G4double beta2 = (gam2-1.0)/gam2;
  G4double g12 = (gam+1.0)*(gam+1.0);
  G4double amol = ((gam-1.0)/gam)*((gam-1.0)/gam);
  //Bhabha coefficients
  G4double bha1 = amol*(2.0*g12-1.0)/(gam2-1.0);
  G4double bha2 = amol*(3.0+1.0/g12);
  G4double bha3 = amol*2.0*gam*(gam-1.0)/g12;
  G4double bha4 = amol*(gam-1.0)*(gam-1.0)/g12;

  //
  //Partial cross section of the active oscillator
  //
  G4double resEne = (*theTable)[targetOscillator]->GetResonanceEnergy();
  G4double invResEne = 1.0/resEne;
  G4double ionEne = (*theTable)[targetOscillator]->GetIonisationEnergy();
  G4double cutoffEne = (*theTable)[targetOscillator]->GetCutoffRecoilResonantEnergy();

  G4double XHDL = 0.;
  G4double XHDT = 0.;
  G4double QM = 0.;
  G4double cps = 0.;
  G4double cp = 0.;

  //Distant excitations XS (same as for electrons)
  if (resEne > cutEnergy && resEne < kineticEnergy)
    {
      cps = kineticEnergy*rb;
      cp = std::sqrt(cps);
      G4double XHDT0 = std::max(std::log(gam2)-beta2-delta,0.);
      if (resEne > 1.0e-6*kineticEnergy)
	{
	  G4double cpp = std::sqrt((kineticEnergy-resEne)*(kineticEnergy-resEne+2.0*electron_mass_c2));
	  QM = std::sqrt((cp-cpp)*(cp-cpp)+electron_mass_c2*electron_mass_c2)-electron_mass_c2;
	}
      else
	{
	  QM = resEne*resEne/(beta2*2.0*electron_mass_c2);
	  QM *= (1.0-QM*0.5/electron_mass_c2);
	}
      if (QM < cutoffEne)
	{
	  XHDL = std::log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)))
	    *invResEne;
	  XHDT = XHDT0*invResEne;	  
	}
      else
	{
	  QM = cutoffEne;
	  XHDL = 0.;
	  XHDT = 0.;
	}
    }
  else
    {
      QM = cutoffEne;
      cps = 0.;
      cp = 0.;
      XHDL = 0.;
      XHDT = 0.;
    }
  //Close collisions (Bhabha)
  G4double wmaxc = kineticEnergy;
  G4double wcl = std::max(cutEnergy,cutoffEne);
  G4double rcl = wcl/kineticEnergy;
  G4double XHC = 0.;
  if (wcl < wmaxc)
    {
      G4double rl1 = 1.0-rcl;
      XHC = ((1.0/rcl-1.0)+bha1*std::log(rcl)+bha2*rl1
	     + (bha3/2.0)*(rcl*rcl-1.0) 
	     + (bha4/3.0)*(1.0-rcl*rcl*rcl))/kineticEnergy;
    }

  //Total cross section per molecule for the active shell, in cm2
  G4double XHTOT = XHC + XHDL + XHDT;

  //very small cross section, do nothing
  if (XHTOT < 1.e-14*barn)
    {
      kineticEnergy1=kineticEnergy;
      cosThetaPrimary=1.0;
      energySecondary=0.0;
      cosThetaSecondary=1.0;
      targetOscillator = numberOfOscillators-1;
      return;
    }

  //decide which kind of interaction we'll have
  TST = XHTOT*G4UniformRand();
  
  // Hard close collision
  G4double TS1 = XHC;
  if (TST < TS1)
    {
      G4double rl1 = 1.0-rcl;
      G4double rk=0.;
      G4bool loopAgain = false;
      do
	{
	  loopAgain = false;
	  rk = rcl/(1.0-G4UniformRand()*rl1);
	  G4double phi = 1.0-rk*(bha1-rk*(bha2-rk*(bha3-bha4*rk)));
	  if (G4UniformRand() > phi)
	    loopAgain = true;
	}while(loopAgain);
      //energy and scattering angle (primary electron)
      G4double deltaE = rk*kineticEnergy;      
      kineticEnergy1 = kineticEnergy - deltaE;
      cosThetaPrimary = std::sqrt(kineticEnergy1*rb/(kineticEnergy*(rb-deltaE)));
      //energy and scattering angle of the delta ray
      energySecondary = deltaE - ionEne; //subtract ionisation energy
      cosThetaSecondary= std::sqrt(deltaE*rb/(kineticEnergy*(deltaE+2.0*electron_mass_c2)));
      if (verboseLevel > 3)
	G4cout << "SampleFinalStatePositron: sampled close collision " << G4endl;
      return;     			     
    }

  //Hard distant longitudinal collisions
  TS1 += XHDL;
  G4double deltaE = resEne;
  kineticEnergy1 = kineticEnergy - deltaE;
  if (TST < TS1)
    {
      G4double QS = QM/(1.0+QM*0.5/electron_mass_c2);
      G4double Q = QS/(std::pow((QS/cutoffEne)*(1.0+cutoffEne*0.5/electron_mass_c2),G4UniformRand())
		       - (QS*0.5/electron_mass_c2));
      G4double QTREV = Q*(Q+2.0*electron_mass_c2);
      G4double cpps = kineticEnergy1*(kineticEnergy1+2.0*electron_mass_c2);
      cosThetaPrimary = (cpps+cps-QTREV)/(2.0*cp*std::sqrt(cpps));
      if (cosThetaPrimary > 1.) 
	cosThetaPrimary = 1.0;
      //energy and emission angle of the delta ray
      energySecondary = deltaE - ionEne;
      cosThetaSecondary = 0.5*(deltaE*(kineticEnergy+rb-deltaE)+QTREV)/std::sqrt(cps*QTREV);
      if (cosThetaSecondary > 1.0)
	cosThetaSecondary = 1.0;
      if (verboseLevel > 3)
	G4cout << "SampleFinalStatePositron: sampled distant longitudinal collision " << G4endl;
      return;      
    }

  //Hard distant transverse collisions
  cosThetaPrimary = 1.0;
  //energy and emission angle of the delta ray
  energySecondary = deltaE - ionEne;
  cosThetaSecondary = 0.5;

  if (verboseLevel > 3)    
    G4cout << "SampleFinalStatePositron: sampled distant transverse collision " << G4endl;
    
  return;
}

