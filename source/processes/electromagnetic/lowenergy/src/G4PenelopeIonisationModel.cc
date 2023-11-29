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
// Author: Luciano Pandola
//
// History:
// --------
// 27 Jul 2010   L Pandola    First complete implementation
// 18 Jan 2011   L.Pandola    Stricter check on production of sub-treshold delta-rays. 
//                            Should never happen now
// 01 Feb 2011   L Pandola  Suppress fake energy-violation warning when Auger is active.
//                          Make sure that fluorescence/Auger is generated only if 
//                          above threshold
// 25 May 2011   L Pandola  Renamed (make v2008 as default Penelope)
// 26 Jan 2012   L Pandola  Migration of AtomicDeexcitation to the new interface 
// 09 Mar 2012   L Pandola  Moved the management and calculation of
//                          cross sections to a separate class. Use a different method to 
//                          get normalized shell cross sections
// 07 Oct 2013   L. Pandola Migration to MT      
// 23 Jun 2015   L. Pandola Keep track of the PIXE flag, to avoid double-production of 
//                           atomic de-excitation (bug #1761)                  
// 29 Aug 2018   L. Pandola Fix bug causing energy non-conservation
//

#include "G4PenelopeIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeOscillator.hh"
#include "G4PenelopeCrossSection.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh" 
#include "G4LossTableManager.hh"
#include "G4PenelopeIonisationXSHandler.hh"
#include "G4EmParameters.hh"
#include "G4AutoLock.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
namespace { G4Mutex  PenelopeIonisationModelMutex = G4MUTEX_INITIALIZER; }
 
G4PenelopeIonisationModel::G4PenelopeIonisationModel(const G4ParticleDefinition* part,
						     const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),fParticle(nullptr),
   fCrossSectionHandler(nullptr), 
   fAtomDeexcitation(nullptr), fKineticEnergy1(0.*eV),
   fCosThetaPrimary(1.0),fEnergySecondary(0.*eV),
   fCosThetaSecondary(0.0),fTargetOscillator(-1),
   fIsInitialised(false),fPIXEflag(false),fLocalTable(false)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  fNBins = 200;

  if (part)
    SetParticle(part);

  //
  fOscManager = G4PenelopeOscillatorManager::GetOscillatorManager();
  //
  fVerboseLevel= 0;   
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  // Atomic deexcitation model activated by default
  SetDeexcitationFlag(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4PenelopeIonisationModel::~G4PenelopeIonisationModel()
{
  if (IsMaster() || fLocalTable)
    {
      if (fCrossSectionHandler)
	delete fCrossSectionHandler;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::Initialise(const G4ParticleDefinition* particle,
					   const G4DataVector& theCuts)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling G4PenelopeIonisationModel::Initialise()" << G4endl;

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  //Issue warning if the AtomicDeexcitation has not been declared
  if (!fAtomDeexcitation)
    {
      G4cout << G4endl;
      G4cout << "WARNING from G4PenelopeIonisationModel " << G4endl;
      G4cout << "Atomic de-excitation module is not instantiated, so there will not be ";
      G4cout << "any fluorescence/Auger emission." << G4endl;
      G4cout << "Please make sure this is intended" << G4endl;
    }

  if (fAtomDeexcitation)
    fPIXEflag = fAtomDeexcitation->IsPIXEActive();

  //If the PIXE flag is active, the PIXE interface will take care of the 
  //atomic de-excitation. The model does not need to do that.
  //Issue warnings here
  if (fPIXEflag && IsMaster() && particle==G4Electron::Electron())
    {
      G4String theModel = G4EmParameters::Instance()->PIXEElectronCrossSectionModel();
      G4cout << "======================================================================" << G4endl;
      G4cout << "The G4PenelopeIonisationModel is being used with the PIXE flag ON." << G4endl;
      G4cout << "Atomic de-excitation will be produced statistically by the PIXE " << G4endl;
      G4cout << "interface by using the shell cross section --> " << theModel << G4endl;
      G4cout << "The built-in model procedure for atomic de-excitation is disabled. " << G4endl;
      G4cout << "*Please be sure this is intended*, or disable PIXE by" << G4endl;
      G4cout << "/process/em/pixe false" << G4endl;    
      G4cout << "======================================================================" << G4endl;
    }

  SetParticle(particle);

  //Only the master model creates/manages the tables. All workers get the 
  //pointer to the table, and use it as readonly
  if (IsMaster() && particle == fParticle)
    {
      //Set the number of bins for the tables. 20 points per decade
      fNBins = (std::size_t) (20*std::log10(HighEnergyLimit()/LowEnergyLimit()));
      fNBins = std::max(fNBins,(std::size_t)100);
      
      //Clear and re-build the tables
      if (fCrossSectionHandler)
	{
	  delete fCrossSectionHandler;
	  fCrossSectionHandler = 0;
	}
      fCrossSectionHandler = new G4PenelopeIonisationXSHandler(fNBins);
      fCrossSectionHandler->SetVerboseLevel(fVerboseLevel);

      //Build tables for all materials
      G4ProductionCutsTable* theCoupleTable = 
	G4ProductionCutsTable::GetProductionCutsTable();      
      for (G4int i=0;i<(G4int)theCoupleTable->GetTableSize();++i)
	{
	  const G4Material* theMat = 
	    theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
	  //Forces the building of the cross section tables
	  fCrossSectionHandler->BuildXSTable(theMat,theCuts.at(i),particle,
					       IsMaster());	 
	}

      if (fVerboseLevel > 2) {
	G4cout << "Penelope Ionisation model v2008 is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / keV << " keV - "
	       << HighEnergyLimit() / GeV << " GeV. Using " 
	       << fNBins << " bins." 
	       << G4endl;
      }
    }
 
  if(fIsInitialised) 
    return;
  fParticleChange = GetParticleChangeForLoss();
  fIsInitialised = true;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4PenelopeIonisationModel::InitialiseLocal(const G4ParticleDefinition* part,
						 G4VEmModel *masterModel)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopeIonisationModel::InitialiseLocal()" << G4endl;
  //
  //Check that particle matches: one might have multiple master models (e.g. 
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeIonisationModel* theModel = 
	static_cast<G4PenelopeIonisationModel*> (masterModel);
      
      //Copy pointers to the data tables
      fCrossSectionHandler = theModel->fCrossSectionHandler;
      
      //copy data
      fNBins = theModel->fNBins;
      
      //Same verbosity for all workers, as the master
      fVerboseLevel = theModel->fVerboseLevel;
    }
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::CrossSectionPerVolume(const G4Material* material,
							  const G4ParticleDefinition* 
							  theParticle,
							  G4double energy,
							  G4double cutEnergy,
							  G4double)
{
  // Penelope model v2008 to calculate the cross section for inelastic collisions above the
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
  if (fVerboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4PenelopeIonisationModel" << G4endl;
 
  SetupForMaterial(theParticle, material, energy);
   
  G4double totalCross = 0.0;
  G4double crossPerMolecule = 0.;

  //Either Initialize() was not called, or we are in a slave and InitializeLocal() was 
  //not invoked
  if (!fCrossSectionHandler)
    {
      //create a **thread-local** version of the table. Used only for G4EmCalculator and 
      //Unit Tests
      fLocalTable = true;
      fCrossSectionHandler = new G4PenelopeIonisationXSHandler(fNBins);	
    }
  
  const G4PenelopeCrossSection* theXS = 
    fCrossSectionHandler->GetCrossSectionTableForCouple(theParticle,
							  material,
							  cutEnergy);
  if (!theXS)
    {
      //If we are here, it means that Initialize() was inkoved, but the MaterialTable was 
      //not filled up. This can happen in a UnitTest or via G4EmCalculator
      if (fVerboseLevel > 0)
	{
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to retrieve the cross section table for " << 
	    theParticle->GetParticleName() << 
	    " in " << material->GetName() << ", cut = " << cutEnergy/keV << " keV " << G4endl;
	  ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	  G4Exception("G4PenelopeIonisationModel::CrossSectionPerVolume()",
		      "em2038",JustWarning,ed);
	}
      //protect file reading via autolock
      G4AutoLock lock(&PenelopeIonisationModelMutex);
      fCrossSectionHandler->BuildXSTable(material,cutEnergy,theParticle);
      lock.unlock();
      //now it should be ok
      theXS = 
	fCrossSectionHandler->GetCrossSectionTableForCouple(theParticle,
							  material,
							  cutEnergy);
    }

  if (theXS)
    crossPerMolecule = theXS->GetHardCrossSection(energy);

  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double atPerMol =  fOscManager->GetAtomsPerMolecule(material);
 
  if (fVerboseLevel > 3)
    G4cout << "Material " << material->GetName() << " has " << atPerMol <<
      "atoms per molecule" << G4endl;

  G4double moleculeDensity = 0.; 
  if (atPerMol)
    moleculeDensity = atomDensity/atPerMol;
  G4double crossPerVolume = crossPerMolecule*moleculeDensity;

  if (fVerboseLevel > 2)
  {
    G4cout << "G4PenelopeIonisationModel " << G4endl;
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
G4double G4PenelopeIonisationModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
							       G4double,
							       G4double,
							       G4double,
							       G4double,
							       G4double)
{
  G4cout << "*** G4PenelopeIonisationModel -- WARNING ***" << G4endl;
  G4cout << "Penelope Ionisation model v2008 does not calculate cross section _per atom_ " << G4endl;
  G4cout << "so the result is always zero. For physics values, please invoke " << G4endl;
  G4cout << "GetCrossSectionPerVolume() or GetMeanFreePath() via the G4EmCalculator" << G4endl;
  return 0;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::ComputeDEDXPerVolume(const G4Material* material,
							 const G4ParticleDefinition* theParticle,
							 G4double kineticEnergy,
							 G4double cutEnergy)
{
  // Penelope model v2008 to calculate the stopping power for soft inelastic collisions
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
 
  if (fVerboseLevel > 3)
    G4cout << "Calling ComputeDEDX() of G4PenelopeIonisationModel" << G4endl;
 
  //Either Initialize() was not called, or we are in a slave and InitializeLocal() was 
  //not invoked
  if (!fCrossSectionHandler)
    {
      //create a **thread-local** version of the table. Used only for G4EmCalculator and 
      //Unit Tests
      fLocalTable = true;
      fCrossSectionHandler = new G4PenelopeIonisationXSHandler(fNBins);	
    }

  const G4PenelopeCrossSection* theXS = 
    fCrossSectionHandler->GetCrossSectionTableForCouple(theParticle,material,
							  cutEnergy);
  if (!theXS)
    {
      //If we are here, it means that Initialize() was inkoved, but the MaterialTable was 
      //not filled up. This can happen in a UnitTest or via G4EmCalculator
      if (fVerboseLevel > 0)
	{   
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to retrieve the cross section table for " << 
	    theParticle->GetParticleName() <<
	    " in " << material->GetName() << ", cut = " << cutEnergy/keV << " keV " << G4endl;
	  ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	  G4Exception("G4PenelopeIonisationModel::ComputeDEDXPerVolume()",
		      "em2038",JustWarning,ed);
	}
      //protect file reading via autolock
      G4AutoLock lock(&PenelopeIonisationModelMutex);
      fCrossSectionHandler->BuildXSTable(material,cutEnergy,theParticle);
      lock.unlock();
      //now it should be ok
      theXS = 
	fCrossSectionHandler->GetCrossSectionTableForCouple(theParticle,
							  material,
							  cutEnergy);
    }

  G4double sPowerPerMolecule = 0.0;
  if (theXS)
    sPowerPerMolecule = theXS->GetSoftStoppingPower(kineticEnergy);

  G4double atomDensity = material->GetTotNbOfAtomsPerVolume();
  G4double atPerMol =  fOscManager->GetAtomsPerMolecule(material);
                                                                                        
  G4double moleculeDensity = 0.; 
  if (atPerMol)
    moleculeDensity = atomDensity/atPerMol;
  G4double sPowerPerVolume = sPowerPerMolecule*moleculeDensity;

  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopeIonisationModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " <<
        kineticEnergy/keV << " keV = " << 
	sPowerPerVolume/(keV/mm) << " keV/mm" << G4endl;
    }
  return sPowerPerVolume;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::MinEnergyCut(const G4ParticleDefinition*,
						 const G4MaterialCutsCouple*)
{
  return fIntrinsicLowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						  const G4MaterialCutsCouple* couple,
						  const G4DynamicParticle* aDynamicParticle,
						  G4double cutE, G4double)
{
  // Penelope model v2008 to sample the final state following an hard inelastic interaction.
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

  if (fVerboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeIonisationModel" << G4endl;
 
  G4double kineticEnergy0 = aDynamicParticle->GetKineticEnergy();
  const G4ParticleDefinition* theParticle = aDynamicParticle->GetDefinition();
 
  if (kineticEnergy0 <= fIntrinsicLowEnergyLimit)
  {
    fParticleChange->SetProposedKineticEnergy(0.);
    fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy0);
    return ;
  }

  const G4Material* material = couple->GetMaterial();
  const G4PenelopeOscillatorTable* theTable = fOscManager->GetOscillatorTableIonisation(material); 

  G4ParticleMomentum particleDirection0 = aDynamicParticle->GetMomentumDirection();
 
  //Initialise final-state variables. The proper values will be set by the methods
  // SampleFinalStateElectron() and SampleFinalStatePositron()
  fKineticEnergy1=kineticEnergy0;
  fCosThetaPrimary=1.0;
  fEnergySecondary=0.0;
  fCosThetaSecondary=1.0;
  fTargetOscillator = -1;
     
  if (theParticle == G4Electron::Electron())
    SampleFinalStateElectron(material,cutE,kineticEnergy0);
  else if (theParticle == G4Positron::Positron())
    SampleFinalStatePositron(material,cutE,kineticEnergy0);
  else
    {
      G4ExceptionDescription ed;
      ed << "Invalid particle " << theParticle->GetParticleName() << G4endl;
      G4Exception("G4PenelopeIonisationModel::SamplingSecondaries()",
		  "em0001",FatalException,ed);
      
    }
  if (fEnergySecondary == 0) return;

  if (fVerboseLevel > 3)
  {
     G4cout << "G4PenelopeIonisationModel::SamplingSecondaries() for " << 
	theParticle->GetParticleName() << G4endl;
      G4cout << "Final eKin = " << fKineticEnergy1 << " keV" << G4endl;
      G4cout << "Final cosTheta = " << fCosThetaPrimary << G4endl;
      G4cout << "Delta-ray eKin = " << fEnergySecondary << " keV" << G4endl;
      G4cout << "Delta-ray cosTheta = " << fCosThetaSecondary << G4endl;
      G4cout << "Oscillator: " << fTargetOscillator << G4endl;
   }
   
  //Update the primary particle
  G4double sint = std::sqrt(1. - fCosThetaPrimary*fCosThetaPrimary);
  G4double phiPrimary  = twopi * G4UniformRand();
  G4double dirx = sint * std::cos(phiPrimary);
  G4double diry = sint * std::sin(phiPrimary);
  G4double dirz = fCosThetaPrimary;
 
  G4ThreeVector electronDirection1(dirx,diry,dirz);
  electronDirection1.rotateUz(particleDirection0);
   
  if (fKineticEnergy1 > 0)
    {
      fParticleChange->ProposeMomentumDirection(electronDirection1);
      fParticleChange->SetProposedKineticEnergy(fKineticEnergy1);
    }
  else    
    fParticleChange->SetProposedKineticEnergy(0.);
    
  //Generate the delta ray
  G4double ionEnergyInPenelopeDatabase = 
    (*theTable)[fTargetOscillator]->GetIonisationEnergy();
  
  //Now, try to handle fluorescence
  //Notice: merged levels are indicated with Z=0 and flag=30
  G4int shFlag = (*theTable)[fTargetOscillator]->GetShellFlag(); 
  G4int Z = (G4int) (*theTable)[fTargetOscillator]->GetParentZ();

  //initialize here, then check photons created by Atomic-Deexcitation, and the final state e-
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4double bindingEnergy = 0.*eV;  

  const G4AtomicShell* shell = nullptr;
  //Real level
  if (Z > 0 && shFlag<30)
    {
      shell = transitionManager->Shell(Z,shFlag-1);
      bindingEnergy = shell->BindingEnergy();
      //shellId = shell->ShellId();
    }

  //correct the fEnergySecondary to account for the fact that the Penelope 
  //database of ionisation energies is in general (slightly) different 
  //from the fluorescence database used in Geant4.
  fEnergySecondary += ionEnergyInPenelopeDatabase-bindingEnergy;
  
  G4double localEnergyDeposit = bindingEnergy;
  //testing purposes only
  G4double energyInFluorescence = 0;
  G4double energyInAuger = 0; 

  if (fEnergySecondary < 0)
    {
      //It means that there was some problem/mismatch between the two databases. 
      //In this case, the available energy is ok to excite the level according 
      //to the Penelope database, but not according to the Geant4 database
      //Full residual energy is deposited locally
      localEnergyDeposit += fEnergySecondary;
      fEnergySecondary = 0.0;
    }

  //Notice: shell might be nullptr (invalid!) if shFlag=30. Must be protected
  //Disable the built-in de-excitation of the PIXE flag is active. In this 
  //case, the PIXE interface takes care (statistically) of producing the 
  //de-excitation.
  //Now, take care of fluorescence, if required
  if (fAtomDeexcitation && !fPIXEflag && shell)
    {
      G4int index = couple->GetIndex();
      if (fAtomDeexcitation->CheckDeexcitationActiveRegion(index))
	{
	  std::size_t nBefore = fvect->size();
	  fAtomDeexcitation->GenerateParticles(fvect,shell,Z,index);
	  std::size_t nAfter = fvect->size(); 
      
	  if (nAfter>nBefore) //actual production of fluorescence
	    {
	      for (std::size_t j=nBefore;j<nAfter;++j) //loop on products
		{
		  G4double itsEnergy = ((*fvect)[j])->GetKineticEnergy();
                  if (itsEnergy < localEnergyDeposit) // valid secondary, generate it
                    {
		      localEnergyDeposit -= itsEnergy;
		      if (((*fvect)[j])->GetParticleDefinition() == G4Gamma::Definition())
		        energyInFluorescence += itsEnergy;
		      else if (((*fvect)[j])->GetParticleDefinition() == G4Electron::Definition())
		        energyInAuger += itsEnergy;
                    }
                  else //invalid secondary: takes more than the available energy: delete it
		    {
		      delete (*fvect)[j];
		      (*fvect)[j] = nullptr;
		    }		    
		}
	    }
	}
    }

  // Generate the delta ray --> to be done only if above cut
  if (fEnergySecondary > cutE)
    {
      G4DynamicParticle* electron = nullptr;
      G4double sinThetaE = std::sqrt(1.-fCosThetaSecondary*fCosThetaSecondary);
      G4double phiEl = phiPrimary+pi; //pi with respect to the primary electron/positron
      G4double xEl = sinThetaE * std::cos(phiEl);
      G4double yEl = sinThetaE * std::sin(phiEl);
      G4double zEl = fCosThetaSecondary;
      G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
      eDirection.rotateUz(particleDirection0);
      electron = new G4DynamicParticle (G4Electron::Electron(),
					eDirection,fEnergySecondary) ;
      fvect->push_back(electron);
    }
  else
    {
      localEnergyDeposit += fEnergySecondary;
      fEnergySecondary = 0;
    }

  if (localEnergyDeposit < 0) //Should not be: issue a G4Exception (warning)
    {
      G4Exception("G4PenelopeIonisationModel::SampleSecondaries()",
		  "em2099",JustWarning,"WARNING: Negative local energy deposit");
      localEnergyDeposit=0.;
    }
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (fVerboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeIonisation" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << fKineticEnergy1/keV << " keV" << G4endl;
      G4cout << "Delta ray " << fEnergySecondary/keV << " keV" << G4endl;
      if (energyInFluorescence)
	G4cout << "Fluorescence x-rays: " << energyInFluorescence/keV << " keV" << G4endl;
      if (energyInAuger)
	G4cout << "Auger electrons: " << energyInAuger/keV << " keV" << G4endl;     
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (fEnergySecondary+energyInFluorescence+fKineticEnergy1+
					  localEnergyDeposit+energyInAuger)/keV <<
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
      
  if (fVerboseLevel > 0)
    {
      G4double energyDiff = std::fabs(fEnergySecondary+energyInFluorescence+fKineticEnergy1+
				      localEnergyDeposit+energyInAuger-kineticEnergy0);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeIonisation: problem with energy conservation: " <<  
	  (fEnergySecondary+energyInFluorescence+fKineticEnergy1+localEnergyDeposit+energyInAuger)/keV <<
	  " keV (final) vs. " <<
	  kineticEnergy0/keV << " keV (initial)" << G4endl;      
    }
    
}
  						       
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PenelopeIonisationModel::SampleFinalStateElectron(const G4Material* mat,
							 G4double cutEnergy,
							 G4double kineticEnergy)
{
  // This method sets the final ionisation parameters
  // fKineticEnergy1, fCosThetaPrimary (= updates of the primary e-)
  // fEnergySecondary, fCosThetaSecondary (= info of delta-ray)
  // fTargetOscillator (= ionised oscillator)
  //
  // The method implements SUBROUTINE EINa of Penelope
  //

  G4PenelopeOscillatorTable* theTable = fOscManager->GetOscillatorTableIonisation(mat);
  std::size_t numberOfOscillators = theTable->size();
  const G4PenelopeCrossSection* theXS = 
    fCrossSectionHandler->GetCrossSectionTableForCouple(G4Electron::Electron(),mat,
							  cutEnergy);
  G4double delta = fCrossSectionHandler->GetDensityCorrection(mat,kineticEnergy);
 
  // Selection of the active oscillator
  G4double TST = G4UniformRand();
  fTargetOscillator = G4int(numberOfOscillators-1); //initialization, last oscillator
  G4double XSsum = 0.;

  for (std::size_t i=0;i<numberOfOscillators-1;++i)
    {     
      XSsum += theXS->GetNormalizedShellCrossSection(i,kineticEnergy); 
     	
      if (XSsum > TST)
	{
	  fTargetOscillator = (G4int) i;
	  break;
	}
    }

  if (fVerboseLevel > 3)
    {
      G4cout << "SampleFinalStateElectron: sampled oscillator #" << 
	fTargetOscillator << "." << G4endl;
      G4cout << "Ionisation energy: " << 
	(*theTable)[fTargetOscillator]->GetIonisationEnergy()/eV << 
	" eV " << G4endl;
      G4cout << "Resonance energy: : " << 
	(*theTable)[fTargetOscillator]->GetResonanceEnergy()/eV << " eV "
	     << G4endl;
    }
  //Constants
  G4double rb = kineticEnergy + 2.0*electron_mass_c2;
  G4double gam = 1.0+kineticEnergy/electron_mass_c2;
  G4double gam2 = gam*gam;
  G4double beta2 = (gam2-1.0)/gam2;
  G4double amol = ((gam-1.0)/gam)*((gam-1.0)/gam);

  //Partial cross section of the active oscillator
  G4double resEne = (*theTable)[fTargetOscillator]->GetResonanceEnergy();
  G4double invResEne = 1.0/resEne;
  G4double ionEne = (*theTable)[fTargetOscillator]->GetIonisationEnergy();
  G4double cutoffEne = (*theTable)[fTargetOscillator]->GetCutoffRecoilResonantEnergy();
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
      G4double XHDT0 = std::max(G4Log(gam2)-beta2-delta,0.);
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
	  XHDL = G4Log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)))
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
	     (1.0-amol)*G4Log(rcl*rrl1))/EE;
    }

  //Total cross section per molecule for the active shell, in cm2
  G4double XHTOT = XHC + XHDL + XHDT;

  //very small cross section, do nothing
  if (XHTOT < 1.e-14*barn)
    {
      fKineticEnergy1=kineticEnergy;
      fCosThetaPrimary=1.0;
      fEnergySecondary=0.0;
      fCosThetaSecondary=1.0;
      fTargetOscillator = G4int(numberOfOscillators-1);
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
      
      fKineticEnergy1 = kineticEnergy - deltaE;
      fCosThetaPrimary = std::sqrt(fKineticEnergy1*rb/(kineticEnergy*(rb-deltaE)));
      //energy and scattering angle of the delta ray
      fEnergySecondary = deltaE - ionEne; //subtract ionisation energy
      fCosThetaSecondary= std::sqrt(deltaE*rb/(kineticEnergy*(deltaE+2.0*electron_mass_c2)));
      if (fVerboseLevel > 3)
	G4cout << "SampleFinalStateElectron: sampled close collision " << G4endl;
      return;     			     
    }

  //Hard distant longitudinal collisions
  TS1 += XHDL;
  G4double deltaE = resEne;
  fKineticEnergy1 = kineticEnergy - deltaE;

  if (TST < TS1)
    {
      G4double QS = QM/(1.0+QM*0.5/electron_mass_c2);
      G4double Q = QS/(std::pow((QS/cutoffEne)*(1.0+cutoffEne*0.5/electron_mass_c2),G4UniformRand())
		       - (QS*0.5/electron_mass_c2));
      G4double QTREV = Q*(Q+2.0*electron_mass_c2);
      G4double cpps = fKineticEnergy1*(fKineticEnergy1+2.0*electron_mass_c2);
      fCosThetaPrimary = (cpps+cps-QTREV)/(2.0*cp*std::sqrt(cpps));
      if (fCosThetaPrimary > 1.) 
	fCosThetaPrimary = 1.0;
      //energy and emission angle of the delta ray
      fEnergySecondary = deltaE - ionEne;
      fCosThetaSecondary = 0.5*(deltaE*(kineticEnergy+rb-deltaE)+QTREV)/std::sqrt(cps*QTREV);
      if (fCosThetaSecondary > 1.0)
	fCosThetaSecondary = 1.0;
      if (fVerboseLevel > 3)
	G4cout << "SampleFinalStateElectron: sampled distant longitudinal collision " << G4endl;
      return;      
    }

  //Hard distant transverse collisions
  fCosThetaPrimary = 1.0;
  //energy and emission angle of the delta ray
  fEnergySecondary = deltaE - ionEne;
  fCosThetaSecondary = 0.5;
  if (fVerboseLevel > 3)
    G4cout << "SampleFinalStateElectron: sampled distant transverse collision " << G4endl;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::SampleFinalStatePositron(const G4Material* mat,
							 G4double cutEnergy,
							 G4double kineticEnergy)
{
  // This method sets the final ionisation parameters
  // fKineticEnergy1, fCosThetaPrimary (= updates of the primary e-)
  // fEnergySecondary, fCosThetaSecondary (= info of delta-ray)
  // fTargetOscillator (= ionised oscillator)
  //
  // The method implements SUBROUTINE PINa of Penelope
  // 
 
  G4PenelopeOscillatorTable* theTable = fOscManager->GetOscillatorTableIonisation(mat);
  std::size_t numberOfOscillators = theTable->size();
  const G4PenelopeCrossSection* theXS = 
    fCrossSectionHandler->GetCrossSectionTableForCouple(G4Positron::Positron(),mat,
							cutEnergy);
  G4double delta = fCrossSectionHandler->GetDensityCorrection(mat,kineticEnergy);

  // Selection of the active oscillator
  G4double TST = G4UniformRand();
  fTargetOscillator = G4int(numberOfOscillators-1); //initialization, last oscillator
  G4double XSsum = 0.;
  for (std::size_t i=0;i<numberOfOscillators-1;++i)
    {
      XSsum += theXS->GetNormalizedShellCrossSection(i,kineticEnergy);   	
      if (XSsum > TST)
	{
	  fTargetOscillator = (G4int) i;
	  break;
	}
    }

  if (fVerboseLevel > 3)
    {
      G4cout << "SampleFinalStatePositron: sampled oscillator #" << 
	fTargetOscillator << "." << G4endl;
      G4cout << "Ionisation energy: " << (*theTable)[fTargetOscillator]->GetIonisationEnergy()/eV 
	     << " eV " << G4endl; 
      G4cout << "Resonance energy: : " << (*theTable)[fTargetOscillator]->GetResonanceEnergy()/eV 
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
  G4double resEne = (*theTable)[fTargetOscillator]->GetResonanceEnergy();
  G4double invResEne = 1.0/resEne;
  G4double ionEne = (*theTable)[fTargetOscillator]->GetIonisationEnergy();
  G4double cutoffEne = (*theTable)[fTargetOscillator]->GetCutoffRecoilResonantEnergy();

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
      G4double XHDT0 = std::max(G4Log(gam2)-beta2-delta,0.);
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
	  XHDL = G4Log(cutoffEne*(QM+2.0*electron_mass_c2)/(QM*(cutoffEne+2.0*electron_mass_c2)))
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
      XHC = ((1.0/rcl-1.0)+bha1*G4Log(rcl)+bha2*rl1
	     + (bha3/2.0)*(rcl*rcl-1.0) 
	     + (bha4/3.0)*(1.0-rcl*rcl*rcl))/kineticEnergy;
    }

  //Total cross section per molecule for the active shell, in cm2
  G4double XHTOT = XHC + XHDL + XHDT;

  //very small cross section, do nothing
  if (XHTOT < 1.e-14*barn)
    {
      fKineticEnergy1=kineticEnergy;
      fCosThetaPrimary=1.0;
      fEnergySecondary=0.0;
      fCosThetaSecondary=1.0;
      fTargetOscillator = G4int(numberOfOscillators-1);
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
      fKineticEnergy1 = kineticEnergy - deltaE;
      fCosThetaPrimary = std::sqrt(fKineticEnergy1*rb/(kineticEnergy*(rb-deltaE)));
      //energy and scattering angle of the delta ray
      fEnergySecondary = deltaE - ionEne; //subtract ionisation energy
      fCosThetaSecondary= std::sqrt(deltaE*rb/(kineticEnergy*(deltaE+2.0*electron_mass_c2)));
      if (fVerboseLevel > 3)
	G4cout << "SampleFinalStatePositron: sampled close collision " << G4endl;
      return;     			     
    }

  //Hard distant longitudinal collisions
  TS1 += XHDL;
  G4double deltaE = resEne;
  fKineticEnergy1 = kineticEnergy - deltaE;
  if (TST < TS1)
    {
      G4double QS = QM/(1.0+QM*0.5/electron_mass_c2);
      G4double Q = QS/(std::pow((QS/cutoffEne)*(1.0+cutoffEne*0.5/electron_mass_c2),G4UniformRand())
		       - (QS*0.5/electron_mass_c2));
      G4double QTREV = Q*(Q+2.0*electron_mass_c2);
      G4double cpps = fKineticEnergy1*(fKineticEnergy1+2.0*electron_mass_c2);
      fCosThetaPrimary = (cpps+cps-QTREV)/(2.0*cp*std::sqrt(cpps));
      if (fCosThetaPrimary > 1.) 
	fCosThetaPrimary = 1.0;
      //energy and emission angle of the delta ray
      fEnergySecondary = deltaE - ionEne;
      fCosThetaSecondary = 0.5*(deltaE*(kineticEnergy+rb-deltaE)+QTREV)/std::sqrt(cps*QTREV);
      if (fCosThetaSecondary > 1.0)
	fCosThetaSecondary = 1.0;
      if (fVerboseLevel > 3)
	G4cout << "SampleFinalStatePositron: sampled distant longitudinal collision " << G4endl;
      return;      
    }

  //Hard distant transverse collisions
  fCosThetaPrimary = 1.0;
  //energy and emission angle of the delta ray
  fEnergySecondary = deltaE - ionEne;
  fCosThetaSecondary = 0.5;

  if (fVerboseLevel > 3)    
    G4cout << "SampleFinalStatePositron: sampled distant transverse collision " << G4endl;
    
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeIonisationModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;  
  }
}
