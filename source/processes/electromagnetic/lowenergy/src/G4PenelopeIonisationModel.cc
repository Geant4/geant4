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
// $Id: G4PenelopeIonisationModel.cc,v 1.18 2010-12-01 15:20:35 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 26 Nov 2008   L Pandola    Migration from process to model 
// 17 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - added MinEnergyCut method
//                  - do not change track status
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they might be checked later on
// 21 Oct 2009   L Pandola    Remove un-necessary fUseAtomicDeexcitation flag - now managed by
//                            G4VEmModel::DeexcitationFlag()
//			      Add ActivateAuger() method 
// 15 Mar 2010   L Pandola    Explicitely initialize Auger to false
// 29 Mar 2010   L Pandola    Added a dummy ComputeCrossSectionPerAtom() method issueing a 
//                            warning if users try to access atomic cross sections via 
//                            G4EmCalculator
// 15 Apr 2010   L. Pandola   Implemented model's own version of MinEnergyCut()
// 23 Apr 2010   L. Pandola   Removed InitialiseElementSelectors() call. Useless here and 
//                            triggers fake warning messages
//

#include "G4PenelopeIonisationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4CrossSectionHandler.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopeIonisationModel::G4PenelopeIonisationModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),
   kineticEnergy1(0.*eV),
   cosThetaPrimary(1.0),
   energySecondary(0.*eV),
   cosThetaSecondary(0.0),
   iOsc(-1),
   crossSectionHandler(0),ionizationEnergy(0),
   resonanceEnergy(0),occupationNumber(0),shellFlag(0),
   theXSTable(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
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
  
  //These vectors do not change when materials or cut change.
  //Therefore I can read it at the constructor
  ionizationEnergy = new std::map<G4int,G4DataVector*>;
  resonanceEnergy  = new std::map<G4int,G4DataVector*>;
  occupationNumber = new std::map<G4int,G4DataVector*>;
  shellFlag = new std::map<G4int,G4DataVector*>;

  ReadData(); //Read data from file

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeIonisationModel::~G4PenelopeIonisationModel()
{
  
  if (crossSectionHandler) delete crossSectionHandler;

  if (theXSTable)
    {
      for (size_t i=0; i<theXSTable->size(); i++)
	delete (*theXSTable)[i];
      delete theXSTable;
     }

  
  std::map <G4int,G4DataVector*>::iterator i;
  if (ionizationEnergy)
    {
      for (i=ionizationEnergy->begin();i != ionizationEnergy->end();i++)
	if (i->second) delete i->second;
      delete ionizationEnergy;
    }
  if (resonanceEnergy)
    {
      for (i=resonanceEnergy->begin();i != resonanceEnergy->end();i++)
	if (i->second) delete i->second;
      delete resonanceEnergy;
    }
  if (occupationNumber)
    {
      for (i=occupationNumber->begin();i != occupationNumber->end();i++)
	if (i->second) delete i->second;
      delete occupationNumber;
    }
  if (shellFlag)
    {
      for (i=shellFlag->begin();i != shellFlag->end();i++)
	if (i->second) delete i->second;
      delete shellFlag;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& )
{
  if (verboseLevel > 3)
    G4cout << "Calling G4PenelopeIonisationModel::Initialise()" << G4endl;

  //Delete and re-initialize the cross section handler
  if (crossSectionHandler)
    {
      crossSectionHandler->Clear();
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }

   if (theXSTable)
    {
      for (size_t i=0; i<theXSTable->size(); i++)
	{
	  delete (*theXSTable)[i];
	  (*theXSTable)[i] = 0;
	}
      delete theXSTable;
      theXSTable = 0;
     }

  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Clear();
  G4String crossSectionFile = "NULL";
  
  if (particle == G4Electron::Electron())
     crossSectionFile = "penelope/ion-cs-el-";
  else if (particle == G4Positron::Positron())
     crossSectionFile = "penelope/ion-cs-po-"; 
  crossSectionHandler->LoadData(crossSectionFile);
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata =  
    crossSectionHandler->BuildMeanFreePathForMaterials();
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;

  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for PenelopeIonisationModel" << G4endl;

  if (verboseLevel > 2) {
    G4cout << "Penelope Ionisation model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / keV << " keV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* theParticle,
                                           G4double energy,
                                           G4double cutEnergy,
                                           G4double)		
{  
  // Penelope model to calculate the cross section for inelastic collisions above the 
  // threshold. It makes use of the Generalised Oscillator Strength (GOS) model from 
  //  D. Liljequist, J. Phys. D: Appl. Phys. 16 (1983) 1567 
  //
  // The total cross section (hard+soft) is read from a database file (element per 
  // element), while the ratio hard-to-total is calculated analytically by taking 
  // into account the atomic oscillators coming into the play for a given threshold. 
  // This is done by the method CalculateCrossSectionsRatio().
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
    G4cout << "Calling CrossSectionPerVolume() of G4PenelopeIonisationModel" << G4endl;

  SetupForMaterial(theParticle, material, energy);

  // VI - should be check at initialisation not in run time
  /*
  if (!crossSectionHandler)
    {
      G4cout << "G4PenelopeIonisationModel::CrossSectionPerVolume" << G4endl;
      G4cout << "The cross section handler is not correctly initialized" << G4endl;
      G4Exception();
    }
  */  
  if (!theXSTable)
    {
      if (verboseLevel > 2)
	{
	  G4cout << "G4PenelopeIonisationModel::CrossSectionPerVolume" << G4endl;
	  G4cout << "Going to build Cross Section table " << G4endl;
	}
      theXSTable = new std::vector<G4VEMDataSet*>;
      theXSTable = BuildCrossSectionTable(theParticle); 
    }
  
  G4double totalCross = 0.0; 
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4double electronVolumeDensity =
    material->GetTotNbOfElectPerVolume();  //electron density

  if (verboseLevel > 4)
    G4cout << "Electron volume density of " << material->GetName() << ": " << 
      electronVolumeDensity*cm3 << " electrons/cm3" << G4endl;

  G4int nelm = material->GetNumberOfElements();
  for (G4int i=0; i<nelm; i++) 
    {
      G4int iZ = (G4int) (*theElementVector)[i]->GetZ();
      G4double ratio = CalculateCrossSectionsRatio(energy,cutEnergy,iZ,electronVolumeDensity,
						   theParticle);
      cross += theAtomNumDensityVector[i]*
	crossSectionHandler->FindValue(iZ,energy)*ratio;
      totalCross += theAtomNumDensityVector[i]*
	crossSectionHandler->FindValue(iZ,energy);
    }

  if (verboseLevel > 2)
  {
    G4cout << "G4PenelopeIonisationModel " << G4endl;
    G4cout << "Mean free path for delta emission > " << cutEnergy/keV << " keV at " << 
      energy/keV << " keV = " << (1./cross)/mm << " mm" << G4endl;
    G4cout << "Total free path for ionisation (no threshold) at " << 
      energy/keV << " keV = " << (1./totalCross)/mm << " mm" << G4endl;
  }
  return cross;
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
  G4cout << "Penelope Ionisation model does not calculate cross section _per atom_ " << G4endl;
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
    G4cout << "Calling ComputeDEDX() of G4PenelopeIonisationModel" << G4endl;

  G4double sPower = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  G4double electronVolumeDensity =
    material->GetTotNbOfElectPerVolume();  //electron density

  //Loop on the elements in the material
  G4int nelm = material->GetNumberOfElements();
  for (G4int i=0; i<nelm; i++) 
    {
      G4int iZ = (G4int) (*theElementVector)[i]->GetZ();
    
      //Calculate stopping power contribution from each element
      //Constants
      G4double gamma = 1.0+kineticEnergy/electron_mass_c2;
      G4double gamma2 = gamma*gamma;
      G4double beta2 = (gamma2-1.0)/gamma2;
      G4double constant = pi*classic_electr_radius*classic_electr_radius
	*2.0*electron_mass_c2/beta2;
 
      G4double delta = CalculateDeltaFermi(kineticEnergy,iZ,electronVolumeDensity);
      G4int nbOsc = (G4int) resonanceEnergy->find(iZ)->second->size();
      G4double stoppingPowerForElement = 0.0;
      //Loop on oscillators of element Z
      for (G4int iosc=0;iosc<nbOsc;iosc++)
	{
	  G4double S1 = 0.0;
	  G4double resEnergy = (*(resonanceEnergy->find(iZ)->second))[iosc];
	  if (theParticle == G4Electron::Electron())
	    S1 = ComputeStoppingPowerForElectrons(kineticEnergy,cutEnergy,delta,resEnergy);
	  else if (theParticle == G4Positron::Positron())
	    S1 = ComputeStoppingPowerForPositrons(kineticEnergy,cutEnergy,delta,resEnergy);
	  G4double occupNb = (*(occupationNumber->find(iZ)->second))[iosc];
	  stoppingPowerForElement += occupNb*constant*S1;
	}

      sPower += stoppingPowerForElement*theAtomNumDensityVector[i];
    }
  if (verboseLevel > 2)
    {
      G4cout << "G4PenelopeIonisationModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " << 
	kineticEnergy/keV << " keV = " << sPower/(keV/mm) << " keV/mm" << G4endl;
    }  
  return sPower;
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
    G4cout << "Calling SamplingSecondaries() of G4PenelopeIonisationModel" << G4endl;

  G4double kineticEnergy0 = aDynamicParticle->GetKineticEnergy();
  const G4ParticleDefinition* theParticle = aDynamicParticle->GetDefinition();

  if (kineticEnergy0 <= fIntrinsicLowEnergyLimit)
  {
    fParticleChange->SetProposedKineticEnergy(0.);
    fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy0);
    return ;
  }
  const G4double electronVolumeDensity = 
    couple->GetMaterial()->GetTotNbOfElectPerVolume();  
  G4ParticleMomentum particleDirection0 = aDynamicParticle->GetMomentumDirection();

  //Initialise final-state variables. The proper values will be set by the methods 
  // CalculateDiscreteForElectrons() and CalculateDiscreteForPositrons()
  kineticEnergy1=kineticEnergy0;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;

  // Select randomly one element in the current material
  if (verboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;

  //Use sampler of G4CrossSectionHandler for now
  // atom can be selected effitiantly if element selectors are initialised   
  //const G4Element* anElement = SelectRandomAtom(couple,theParticle,kineticEnergy0);

  G4int iZ = SampleRandomAtom(couple,kineticEnergy0);
  
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t indx = couple->GetIndex();
  G4double cutG = (*(theCoupleTable->GetEnergyCutsVector(0)))[indx];
  
  if (verboseLevel > 2)
    G4cout << "Selected Z = " << iZ << G4endl;
 
  // The method CalculateDiscreteForXXXX() set the private variables:
  // kineticEnergy1 = energy of the primary electron/positron after the interaction
  // cosThetaPrimary = std::cos(theta) of the primary after the interaction
  // energySecondary = energy of the secondary electron
  // cosThetaSecondary = std::cos(theta) of the secondary
  if (theParticle == G4Electron::Electron())
    CalculateDiscreteForElectrons(kineticEnergy0,cutE,iZ,electronVolumeDensity);
  else if (theParticle == G4Positron::Positron())
    CalculateDiscreteForPositrons(kineticEnergy0,cutE,iZ,electronVolumeDensity);

  // if (energySecondary == 0) return;

  //Update the primary particle
  G4double sint = std::sqrt(1. - cosThetaPrimary*cosThetaPrimary);
  G4double phi  = twopi * G4UniformRand();
  G4double dirx = sint * std::cos(phi);
  G4double diry = sint * std::sin(phi);
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
  G4int iosc2 = 0;
  G4double ioniEnergy = 0.0;
  if (iOsc > 0) {
    ioniEnergy=(*(ionizationEnergy->find(iZ)->second))[iOsc];
    iosc2 = (ionizationEnergy->find(iZ)->second->size()) - iOsc; //they are in reversed order
  }

  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  G4double bindingEnergy = 0.0;
  G4int shellId = 0;
  if (iOsc > 0)
    {
      const G4AtomicShell* shell = transitionManager->Shell(iZ,iosc2-1); // Modified by Alf
      bindingEnergy = shell->BindingEnergy();
      shellId = shell->ShellId();
    }

  G4double ionEnergy = bindingEnergy; //energy spent to ionise the atom according to G4dabatase
  G4double eKineticEnergy = energySecondary;

  //This is an awful thing: Penelope generates the fluorescence only for L and K shells 
  //(i.e. Osc = 1 --> 4). For high-Z, the other shells can be quite relevant. In this case 
  //one MUST ensure ''by hand'' the energy conservation. Then there is the other problem that 
  //the fluorescence database of Penelope doesn not match that of Geant4.

  G4double energyBalance = kineticEnergy0 - kineticEnergy1 - energySecondary; //Penelope Balance

  if (std::abs(energyBalance) < 1*eV)
    //in this case Penelope didn't subtract the fluorescence energy: do here by hand
    eKineticEnergy = energySecondary - bindingEnergy;    
  else 
    //Penelope subtracted the fluorescence, but one has to match the databases
    eKineticEnergy = energySecondary+ioniEnergy-bindingEnergy;
  
  G4double localEnergyDeposit = ionEnergy; 
  G4double energyInFluorescence = 0.0*eV;

  if(DeexcitationFlag() && iZ > 5) 
    {
      if (ionEnergy > cutG || ionEnergy > cutE)
	{
	  deexcitationManager.SetCutForSecondaryPhotons(cutG);
	  deexcitationManager.SetCutForAugerElectrons(cutE);
	  std::vector<G4DynamicParticle*> *photonVector =
	    deexcitationManager.GenerateParticles(iZ,shellId);
	  //Check for secondaries
          if(photonVector) 
	    {
	      for (size_t k=0;k<photonVector->size();k++)
		{
		  G4DynamicParticle* aPhoton = (*photonVector)[k];
		  if (aPhoton)
		    {
		      G4double itsEnergy = aPhoton->GetKineticEnergy();
		      if (itsEnergy <= localEnergyDeposit)
			{
			  if(aPhoton->GetDefinition() == G4Gamma::Gamma())
			    energyInFluorescence += itsEnergy;
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

  // Generate the delta ray 
  G4double sin2 = std::sqrt(1. - cosThetaSecondary*cosThetaSecondary);
  G4double phi2  = twopi * G4UniformRand();
  
  G4double xEl = sin2 * std::cos(phi2); 
  G4double yEl = sin2 * std::sin(phi2);
  G4double zEl = cosThetaSecondary;
  G4ThreeVector eDirection(xEl,yEl,zEl); //electron direction
  eDirection.rotateUz(particleDirection0);

  G4DynamicParticle* deltaElectron = new G4DynamicParticle (G4Electron::Electron(),
							    eDirection,eKineticEnergy) ;
  fvect->push_back(deltaElectron);
  
  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING-" 
	     << "G4PenelopeIonisationModel::SampleSecondaries - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit=0.;
    }
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);
  
  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeIonisation" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy0/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << kineticEnergy1/keV << " keV" << G4endl;
      G4cout << "Delta ray " << eKineticEnergy/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (eKineticEnergy+energyInFluorescence+kineticEnergy1+
					  localEnergyDeposit)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = std::fabs(eKineticEnergy+energyInFluorescence+kineticEnergy1+
				      localEnergyDeposit-kineticEnergy0);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeIonisation: problem with energy conservation: " << 
	  (eKineticEnergy+energyInFluorescence+kineticEnergy1+localEnergyDeposit)/keV << 
	  " keV (final) vs. " << 
	  kineticEnergy0/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::ReadData()
{
  if (verboseLevel > 2)
    {
      G4cout << "Data from G4PenelopeIonisationModel read " << G4endl;
    }
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeIonisationModel - G4LEDATA environment variable not set!";
      G4Exception(excep);
      return;
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/ion-pen.dat";
  std::ifstream file(pathFile);
  
  if (!file.is_open())
    {
      G4String excep = "G4PenelopeIonisationModel - data file " + pathFile + " not found!";
      G4Exception(excep);
    }

 if (!ionizationEnergy || !resonanceEnergy || !occupationNumber || !shellFlag)
    {
      G4String excep = "G4PenelopeIonisationModel: problem with reading data from file";
      G4Exception(excep);
      return;
    }

 G4int Z=1,nLevels=0;
 G4int test,test1;

 do{
   file >> Z >> nLevels;
   //Check for nLevels validity, before using it in a loop
   if (nLevels<0 || nLevels>64)
     {
       G4String excep = "G4PenelopeIonisationModel: corrupted data file ?";
       G4Exception(excep);
       return;
     }
   //Allocate space for storage
   G4DataVector* occVector = new G4DataVector;
   G4DataVector* ionEVector = new G4DataVector;
   G4DataVector* resEVector = new G4DataVector;
   G4DataVector* shellIndVector = new G4DataVector;
   //
   G4double a1,a2,a3,a4;
   G4int k1,k2,k3;
   for (G4int h=0;h<nLevels;h++)
     {
       //index,occup number,ion energy,res energy,fj0,kz,shell flag
       file >> k1 >> a1 >> a2 >> a3 >> a4 >> k2 >> k3;
       //Make explicit unit of measurements for ionisation and resonance energy, 
       // which is MeV
       a2 *= MeV;
       a3 *= MeV;
       //
       occVector->push_back(a1);
       ionEVector->push_back(a2);
       resEVector->push_back(a3);
       shellIndVector->push_back((G4double) k3);
     }
   //Ok, done for element Z
   occupationNumber->insert(std::make_pair(Z,occVector));
   ionizationEnergy->insert(std::make_pair(Z,ionEVector));
   resonanceEnergy->insert(std::make_pair(Z,resEVector));
   shellFlag->insert(std::make_pair(Z,shellIndVector));
   file >> test >> test1; //-1 -1 close the data for each Z
   if (test > 0) {
     G4String excep = "G4PenelopeIonisationModel - data file corrupted!";
     G4Exception(excep);
   }
 }while (test != -2); //the very last Z is closed with -2 instead of -1
}
				   
					   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::CalculateDeltaFermi(G4double kinEnergy ,G4int Z,
							G4double electronVolumeDensity)
{
  G4double plasmaEnergyCoefficient = 1.377e-39*(MeV*MeV*m3); //(e*hbar)^2/(epsilon0*electron_mass)
  G4double plasmaEnergySquared = plasmaEnergyCoefficient*electronVolumeDensity;
  // std::sqrt(plasmaEnergySquared) is the plasma energy of the solid (MeV)
  G4double gam = 1.0+kinEnergy/electron_mass_c2;
  G4double gam2=gam*gam;
  G4double delta = 0.0;

  //Density effect
  G4double TST = ((G4double) Z)/(gam2*plasmaEnergySquared); 
 
  G4double wl2 = 0.0;
  G4double fdel=0.0;
  G4double wr=0;
  size_t nbOsc = resonanceEnergy->find(Z)->second->size();
  for(size_t i=0;i<nbOsc;i++)
    {
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
      wr = (*(resonanceEnergy->find(Z)->second))[i];
      fdel += occupNb/(wr*wr+wl2);
    }
  if (fdel < TST) return delta;
  G4double help1 = (*(resonanceEnergy->find(Z)->second))[nbOsc-1];
  wl2 = help1*help1;
  do{
    wl2=wl2*2.0;
    fdel = 0.0;
    for (size_t ii=0;ii<nbOsc;ii++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[ii];
      wr = (*(resonanceEnergy->find(Z)->second))[ii];
      fdel += occupNb/(wr*wr+wl2);
    }
  }while (fdel > TST);
  G4double wl2l=0.0;
  G4double wl2u = wl2;
  G4double control = 0.0;
  do{
    wl2=0.5*(wl2l+wl2u);
    fdel = 0.0;
    for (size_t jj=0;jj<nbOsc;jj++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[jj];
      wr = (*(resonanceEnergy->find(Z)->second))[jj];
      fdel += occupNb/(wr*wr+wl2);
    }
    if (fdel > TST)
      wl2l = wl2;
    else
      wl2u = wl2;
    control = wl2u-wl2l-wl2*1e-12; 
  }while(control>0);

  //Density correction effect
   for (size_t kk=0;kk<nbOsc;kk++){
      G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[kk];
      wr = (*(resonanceEnergy->find(Z)->second))[kk];
      delta += occupNb*std::log(1.0+wl2/(wr*wr));
    }
   delta = (delta/((G4double) Z))-wl2/(gam2*plasmaEnergySquared); 
   return delta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4PenelopeIonisationModel::CalculateDiscreteForElectrons(G4double kinEnergy,
							 G4double cutoffEnergy,
							 G4int Z,
							 G4double electronVolumeDensity)
{
  if (verboseLevel > 2)
    G4cout << "Entering in CalculateDiscreteForElectrons() for energy " << 
      kinEnergy/keV << " keV " << G4endl;

  //Initialization of variables to be calculated in this routine
  kineticEnergy1=kinEnergy;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;
  iOsc=-1;

  //constants
  G4double rb=kinEnergy+2.0*electron_mass_c2;
  G4double gamma = 1.0+kinEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double amol = (gamma-1.0)*(gamma-1.0)/gamma2;
  G4double cps = kinEnergy*rb;
  G4double cp = std::sqrt(cps);
  
  G4double delta = CalculateDeltaFermi(kinEnergy,Z,electronVolumeDensity);
  G4double distantTransvCS0 = std::max(std::log(gamma2)-beta2-delta,0.0);

  G4double rl,rl1;

  if (cutoffEnergy > kinEnergy) return; //delta rays are not generated

  G4DataVector* qm = new G4DataVector();
  G4DataVector* cumulHardCS = new G4DataVector();
  std::vector<G4int> *typeOfInteraction = new std::vector<G4int>;
  G4DataVector* nbOfLevel = new G4DataVector();
 
  //Hard close collisions with outer shells
  G4double wmaxc = 0.5*kinEnergy;
  G4double closeCS0 = 0.0;
  G4double closeCS = 0.0;
  if (cutoffEnergy>0.1*eV) 
    {
      rl=cutoffEnergy/kinEnergy;
      rl1=1.0-rl;
      if (rl < 0.5)
	closeCS0 = (amol*(0.5-rl)+(1.0/rl)-(1.0/rl1)+(1.0-amol)*std::log(rl/rl1))/kinEnergy;
    }

  // Cross sections for the different oscillators

  // totalHardCS contains the cumulative hard interaction cross section for the different
  // excitable levels and the different interaction channels (close, distant, etc.),
  // i.e.
  // cumulHardCS[0] = 0.0
  // cumulHardCS[1] = 1st excitable level (distant longitudinal only)
  // cumulHardCS[2] = 1st excitable level (distant longitudinal + transverse)
  // cumulHardCS[3] = 1st excitable level (distant longitudinal + transverse + close)
  // cumulHardCS[4] = 1st excitable level (all channels) + 2nd excitable level (distant long only)
  // etc.
  // This is used for sampling the atomic level which is ionised and the channel of the
  // interaction.
  //
  // For each index iFill of the cumulHardCS vector,
  // nbOfLevel[iFill] contains the current excitable atomic level and
  // typeOfInteraction[iFill] contains the current interaction channel, with the legenda:
  //   1 = distant longitudinal interaction
  //   2 = distant transverse interaction
  //   3 = close collision
  //   4 = close collision with outer shells (in this case nbOfLevel < 0 --> no binding energy)


  G4int nOscil = ionizationEnergy->find(Z)->second->size();
  G4double totalHardCS = 0.0;
  G4double involvedElectrons = 0.0;
  for (G4int i=0;i<nOscil;i++){
    G4double wi = (*(resonanceEnergy->find(Z)->second))[i];
    G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
 
    //Distant excitations
    if (wi>cutoffEnergy && wi<kinEnergy)
      {
	if (wi>(1e-6*kinEnergy))
	  {
	    G4double cpp=std::sqrt((kinEnergy-wi)*(kinEnergy-wi+2.0*electron_mass_c2));
	    qm->push_back(std::sqrt((cp-cpp)*(cp-cpp)+electron_mass_c2*electron_mass_c2)-electron_mass_c2);
	  }
	else
	  {
	    qm->push_back((wi*wi)/(beta2*2.0*electron_mass_c2));
	  }
	if ((*qm)[i] < wi)
	  {
	    
	    G4double distantLongitCS =  occupNb*std::log(wi*((*qm)[i]+2.0*electron_mass_c2)/
					 ((*qm)[i]*(wi+2.0*electron_mass_c2)))/wi;
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(1); //distant longitudinal
	    nbOfLevel->push_back((G4double) i); //only excitable level are counted 
	    totalHardCS += distantLongitCS;
	    
	    G4double distantTransvCS = occupNb*distantTransvCS0/wi;
	    
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(2); //distant tranverse
	    nbOfLevel->push_back((G4double) i);
	    totalHardCS += distantTransvCS;
	  }
      }
    else 
      qm->push_back(wi);
      
    //close collisions
    if(wi < wmaxc)
      {
	if (wi < cutoffEnergy) 
	  involvedElectrons += occupNb;
	else
	  {
	    rl=wi/kinEnergy;
	    rl1=1.0-rl;
	    closeCS = occupNb*(amol*(0.5-rl)+(1.0/rl)-(1.0/rl1)+(1.0-amol)*std::log(rl/rl1))/kinEnergy;
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(3); //close
	    nbOfLevel->push_back((G4double) i);
	    totalHardCS += closeCS;
	  }
      }
  } // loop on the levels
  
  cumulHardCS->push_back(totalHardCS);
  typeOfInteraction->push_back(4); //close interaction with outer shells
  nbOfLevel->push_back(-1.0);
  totalHardCS += involvedElectrons*closeCS0;
  cumulHardCS->push_back(totalHardCS); //this is the final value of the totalHardCS

  if (totalHardCS < 1e-30) {
    kineticEnergy1=kinEnergy;
    cosThetaPrimary=1.0;
    energySecondary=0.0;
    cosThetaSecondary=0.0;
    iOsc=-1;
    delete qm;
    delete cumulHardCS;
    delete typeOfInteraction;
    delete nbOfLevel;
    return;
  }

  //Testing purposes
  if (verboseLevel > 6)
    {
      for (size_t t=0;t<cumulHardCS->size();t++)
	G4cout << (*cumulHardCS)[t] << " " << (*typeOfInteraction)[t] << 
	  " " << (*nbOfLevel)[t] << G4endl;
    }

  //Selection of the active oscillator on the basis of the cumulative cross sections
  G4double TST = totalHardCS*G4UniformRand();
  G4int is=0;
  G4int js= nbOfLevel->size();
  do{
    G4int it=(is+js)/2;
    if (TST > (*cumulHardCS)[it]) is=it;
    if (TST <= (*cumulHardCS)[it]) js=it;
  }while((js-is) > 1);

  G4double UII=0.0;
  G4double rkc=cutoffEnergy/kinEnergy;
  G4double dde;
  G4int kks;

  G4int sampledInteraction = (*typeOfInteraction)[is];
  iOsc = (G4int) (*nbOfLevel)[is];

  if (verboseLevel > 4) 
    G4cout << "Chosen interaction #:" << sampledInteraction << " on level " << iOsc << G4endl;

  //Generates the final state according to the sampled level and 
  //interaction channel
  
  if (sampledInteraction == 1)  //Hard distant longitudinal collisions
    {
      dde= (*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=kinEnergy-dde;
      G4double qs=(*qm)[iOsc]/(1.0+((*qm)[iOsc]/(2.0*electron_mass_c2)));
      G4double q=qs/(std::pow((qs/dde)*(1.0+(0.5*dde/electron_mass_c2)),G4UniformRand())-(0.5*qs/electron_mass_c2));
      G4double qtrev = q*(q+2.0*electron_mass_c2);
      G4double cpps = kineticEnergy1*(kineticEnergy1+2.0*electron_mass_c2);
      cosThetaPrimary = (cpps+cps-qtrev)/(2.0*cp*std::sqrt(cpps));
      if (cosThetaPrimary>1.0) cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      //kks > 4 means that we are in an outer shell
      if (kks>4) 
	energySecondary=dde;
      else
	energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
	
      cosThetaSecondary = 0.5*(dde*(kinEnergy+rb-dde)+qtrev)/std::sqrt(cps*qtrev);
      if (cosThetaSecondary>1.0) cosThetaSecondary=1.0;
    }
  else if (sampledInteraction == 2)  //Hard distant transverse collisions
    {
      dde=(*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=kinEnergy-dde;
      cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      if (kks>4)
	{
	  energySecondary=dde;
	}
      else
	{
	  energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
	}
      cosThetaSecondary = 1.0;
    }
  else if (sampledInteraction == 3 || sampledInteraction == 4) //Close interaction
    {
      if (sampledInteraction == 4) //interaction with inner shells
	{
	  UII=0.0;
	  rkc = cutoffEnergy/kinEnergy;
	  iOsc = -1;
	}
      else
	{
	  kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
	  if (kks > 4) {
	    UII=0.0;
	  }
	  else
	    {
	      UII = (*(ionizationEnergy->find(Z)->second))[iOsc];
	    }
	  rkc = (*(resonanceEnergy->find(Z)->second))[iOsc]/kinEnergy;
	}
    G4double A = 0.5*amol;
    G4double arkc = A*0.5*rkc;
    G4double phi,rk2,rk,rkf;
    do{
      G4double fb = (1.0+arkc)*G4UniformRand();
      if (fb<1.0)
	{
	  rk=rkc/(1.0-fb*(1.0-(rkc*2.0)));
	}
      else{
	rk = rkc+(fb-1.0)*(0.5-rkc)/arkc;
      }
      rk2 = rk*rk;
      rkf = rk/(1.0-rk);
      phi = 1.0+(rkf*rkf)-rkf+amol*(rk2+rkf);
    }while ((G4UniformRand()*(1.0+A*rk2)) > phi);
    //Energy and scattering angle (primary electron);
    kineticEnergy1 = kinEnergy*(1.0-rk);
    cosThetaPrimary = std::sqrt(kineticEnergy1*rb/(kinEnergy*(rb-(rk*kinEnergy))));
    //Energy and scattering angle of the delta ray
    energySecondary = kinEnergy-kineticEnergy1-UII;
    cosThetaSecondary = std::sqrt(rk*kinEnergy*rb/(kinEnergy*(rk*kinEnergy+2.0*electron_mass_c2)));
    }
  else
    {
      G4String excep = "G4PenelopeIonisationModel - Error in the calculation of the final state";
      G4Exception(excep);
    }

  delete qm;
  delete cumulHardCS;
  
  delete typeOfInteraction;
  delete nbOfLevel;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 void 
G4PenelopeIonisationModel::CalculateDiscreteForPositrons(G4double kinEnergy,
							 G4double cutoffEnergy,
							 G4int Z,
							 G4double electronVolumeDensity)
{
  kineticEnergy1=kinEnergy;
  cosThetaPrimary=1.0;
  energySecondary=0.0;
  cosThetaSecondary=1.0;
  iOsc=-1;
  //constants
  G4double rb=kinEnergy+2.0*electron_mass_c2;
  G4double gamma = 1.0+kinEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double amol = (gamma-1.0)*(gamma-1.0)/gamma2;
  G4double cps = kinEnergy*rb;
  G4double cp = std::sqrt(cps);
  G4double help = (gamma+1.0)*(gamma+1.0);
  G4double bha1 = amol*(2.0*help-1.0)/(gamma2-1.0);
  G4double bha2 = amol*(3.0+1.0/help);
  G4double bha3 = amol*2.0*gamma*(gamma-1.0)/help;
  G4double bha4 = amol*(gamma-1.0)*(gamma-1.0)/help;
  
  G4double delta = CalculateDeltaFermi(kinEnergy,Z,electronVolumeDensity);
  G4double distantTransvCS0 = std::max(std::log(gamma2)-beta2-delta,0.0);

  G4double rl,rl1;

  if (cutoffEnergy > kinEnergy) return; //delta rays are not generated

  G4DataVector* qm = new G4DataVector();
  G4DataVector* cumulHardCS = new G4DataVector();
  std::vector<G4int> *typeOfInteraction = new std::vector<G4int>;
  G4DataVector* nbOfLevel = new G4DataVector();
 

  //Hard close collisions with outer shells
  G4double wmaxc = kinEnergy;
  G4double closeCS0 = 0.0;
  G4double closeCS = 0.0;
  if (cutoffEnergy>0.1*eV) 
    {
      rl=cutoffEnergy/kinEnergy;
      rl1=1.0-rl;
      if (rl < 1.0)
	closeCS0 = (((1.0/rl)-1.0) + bha1*std::log(rl) + bha2*rl1
		    + (bha3/2.0)*((rl*rl)-1.0)
		    + (bha4/3.0)*(1.0-(rl*rl*rl)))/kinEnergy;
    }

  // Cross sections for the different oscillators

  // totalHardCS contains the cumulative hard interaction cross section for the different
  // excitable levels and the different interaction channels (close, distant, etc.),
  // i.e.
  // cumulHardCS[0] = 0.0
  // cumulHardCS[1] = 1st excitable level (distant longitudinal only)
  // cumulHardCS[2] = 1st excitable level (distant longitudinal + transverse)
  // cumulHardCS[3] = 1st excitable level (distant longitudinal + transverse + close)
  // cumulHardCS[4] = 1st excitable level (all channels) + 2nd excitable level (distant long only)
  // etc.
  // This is used for sampling the atomic level which is ionised and the channel of the
  // interaction.
  //
  // For each index iFill of the cumulHardCS vector,
  // nbOfLevel[iFill] contains the current excitable atomic level and
  // typeOfInteraction[iFill] contains the current interaction channel, with the legenda:
  //   1 = distant longitudinal interaction
  //   2 = distant transverse interaction
  //   3 = close collision
  //   4 = close collision with outer shells (in this case nbOfLevel < 0 --> no binding energy)


  G4int nOscil = ionizationEnergy->find(Z)->second->size();
  G4double totalHardCS = 0.0;
  G4double involvedElectrons = 0.0;
  for (G4int i=0;i<nOscil;i++){
    G4double wi = (*(resonanceEnergy->find(Z)->second))[i];
    G4int occupNb = (G4int) (*(occupationNumber->find(Z)->second))[i];
    //Distant excitations
    if (wi>cutoffEnergy && wi<kinEnergy)
      {
	if (wi>(1e-6*kinEnergy)){
	  G4double cpp=std::sqrt((kinEnergy-wi)*(kinEnergy-wi+2.0*electron_mass_c2));
	  qm->push_back(std::sqrt((cp-cpp)*(cp-cpp)+ electron_mass_c2 * electron_mass_c2)-electron_mass_c2);
	}
	else
	  {
	    qm->push_back(wi*wi/(beta2+2.0*electron_mass_c2));
	  }
	//verificare che quando arriva qui il vettore ha SEMPRE l'i-esimo elemento
	if ((*qm)[i] < wi)
	  {
	    
	    G4double distantLongitCS =  occupNb*std::log(wi*((*qm)[i]+2.0*electron_mass_c2)/
					 ((*qm)[i]*(wi+2.0*electron_mass_c2)))/wi;
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(1); //distant longitudinal
	    nbOfLevel->push_back((G4double) i); //only excitable level are counted 
	    totalHardCS += distantLongitCS;
	    
	    G4double distantTransvCS = occupNb*distantTransvCS0/wi;
	    
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(2); //distant tranverse
	    nbOfLevel->push_back((G4double) i);
	    totalHardCS += distantTransvCS;
	  }
      }
    else 
      qm->push_back(wi);
      
    //close collisions
    if(wi < wmaxc)
      {
	if (wi < cutoffEnergy) {
	  involvedElectrons += occupNb;
	}
	else
	  {
	    rl=wi/kinEnergy;
	    rl1=1.0-rl;
	    closeCS = occupNb*(((1.0/rl)-1.0)+bha1*std::log(rl)+bha2*rl1
			       + (bha3/2.0)*((rl*rl)-1.0)
			       + (bha4/3.0)*(1.0-(rl*rl*rl)))/kinEnergy;
	    cumulHardCS->push_back(totalHardCS);
	    typeOfInteraction->push_back(3); //close
	    nbOfLevel->push_back((G4double) i);
	    totalHardCS += closeCS;
	  }
      }
  } // loop on the levels
  
  cumulHardCS->push_back(totalHardCS);
  typeOfInteraction->push_back(4); //close interaction with outer shells
  nbOfLevel->push_back(-1.0);
  totalHardCS += involvedElectrons*closeCS0;
  cumulHardCS->push_back(totalHardCS); //this is the final value of the totalHardCS

  if (totalHardCS < 1e-30) {
    kineticEnergy1=kinEnergy;
    cosThetaPrimary=1.0;
    energySecondary=0.0;
    cosThetaSecondary=0.0;
    iOsc=-1;
    delete qm;
    delete cumulHardCS;
    delete typeOfInteraction;
    delete nbOfLevel;
    return;
  }

  //Selection of the active oscillator on the basis of the cumulative cross sections
  G4double TST = totalHardCS*G4UniformRand();
  G4int is=0;
  G4int js= nbOfLevel->size();
  do{
    G4int it=(is+js)/2;
    if (TST > (*cumulHardCS)[it]) is=it;
    if (TST <= (*cumulHardCS)[it]) js=it;
  }while((js-is) > 1);

  G4double UII=0.0;
  G4double rkc=cutoffEnergy/kinEnergy;
  G4double dde;
  G4int kks;

  G4int sampledInteraction = (*typeOfInteraction)[is];
  iOsc = (G4int) (*nbOfLevel)[is];

  //Generates the final state according to the sampled level and 
  //interaction channel
  
  if (sampledInteraction == 1)  //Hard distant longitudinal collisions
    {
      dde= (*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=kinEnergy-dde;
      G4double qs=(*qm)[iOsc]/(1.0+((*qm)[iOsc]/(2.0*electron_mass_c2)));
      G4double q=qs/(std::pow((qs/dde)*(1.0+(0.5*dde/electron_mass_c2)),G4UniformRand())-(0.5*qs/electron_mass_c2));
      G4double qtrev = q*(q+2.0*electron_mass_c2);
      G4double cpps = kineticEnergy1*(kineticEnergy1+2.0*electron_mass_c2);
      cosThetaPrimary = (cpps+cps-qtrev)/(2.0*cp*std::sqrt(cpps));
      if (cosThetaPrimary>1.0) cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      if (kks>4) 
	energySecondary=dde;
	
      else
	energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
      cosThetaSecondary = 0.5*(dde*(kinEnergy+rb-dde)+qtrev)/std::sqrt(cps*qtrev);
      if (cosThetaSecondary>1.0) cosThetaSecondary=1.0;
    }
  else if (sampledInteraction == 2)  //Hard distant transverse collisions
    {
      dde=(*(resonanceEnergy->find(Z)->second))[iOsc];
      kineticEnergy1=kinEnergy-dde;
      cosThetaPrimary=1.0;
      //Energy and emission angle of the delta ray
      kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
      if (kks>4)
	{
	  energySecondary=dde;
	}
      else
	{
	  energySecondary=dde-(*(ionizationEnergy->find(Z)->second))[iOsc];
	}
      cosThetaSecondary = 1.0;
    }
  else if (sampledInteraction == 3 || sampledInteraction == 4) //Close interaction
    {
      if (sampledInteraction == 4) //interaction with inner shells
	{
	  UII=0.0;
	  rkc = cutoffEnergy/kinEnergy;
	  iOsc = -1;
	}
      else
	{
	  kks = (G4int) (*(shellFlag->find(Z)->second))[iOsc];
	  if (kks > 4) {
	    UII=0.0;
	  }
	  else
	    {
	      UII = (*(ionizationEnergy->find(Z)->second))[iOsc];
	    }
	  rkc = (*(resonanceEnergy->find(Z)->second))[iOsc]/kinEnergy;
	}
      G4double phi,rk;
      do{
	rk=rkc/(1.0-G4UniformRand()*(1.0-rkc));
	phi = 1.0-rk*(bha1-rk*(bha2-rk*(bha3-bha4*rk)));
      }while ( G4UniformRand() > phi);
      //Energy and scattering angle (primary electron);
      kineticEnergy1 = kinEnergy*(1.0-rk);
      cosThetaPrimary = std::sqrt(kineticEnergy1*rb/(kinEnergy*(rb-(rk*kinEnergy))));
      //Energy and scattering angle of the delta ray
      energySecondary = kinEnergy-kineticEnergy1-UII;
      cosThetaSecondary = std::sqrt(rk*kinEnergy*rb/(kinEnergy*(rk*kinEnergy+2.0*electron_mass_c2)));
    }      
  else
    {
      G4String excep = "G4PenelopeIonisationModel - Error in the calculation of the final state";
      G4Exception(excep);
    }

  delete qm;
  delete cumulHardCS;
  delete typeOfInteraction;
  delete nbOfLevel;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4PenelopeIonisationModel::CalculateCrossSectionsRatio(G4double kinEnergy,
						       G4double cutoffEnergy,
						       G4int Z, G4double electronVolumeDensity,
						       const G4ParticleDefinition* theParticle)
{
  //Constants
  G4double gamma = 1.0+kinEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double constant = pi*classic_electr_radius*classic_electr_radius*2.0*electron_mass_c2/beta2;
  G4double delta = CalculateDeltaFermi(kinEnergy,Z,electronVolumeDensity);
  G4int nbOsc = (G4int) resonanceEnergy->find(Z)->second->size();
  G4double softCS = 0.0;
  G4double hardCS = 0.0;
  for (G4int i=0;i<nbOsc;i++){
    G4double resEnergy = (*(resonanceEnergy->find(Z)->second))[i];
    std::pair<G4double,G4double> SoftAndHardXS(0.,0.);
    if (theParticle == G4Electron::Electron())
      SoftAndHardXS = CrossSectionsRatioForElectrons(kinEnergy,resEnergy,delta,cutoffEnergy);
    else if (theParticle == G4Positron::Positron())
      SoftAndHardXS = CrossSectionsRatioForPositrons(kinEnergy,resEnergy,delta,cutoffEnergy);      
    G4double occupNb = (*(occupationNumber->find(Z)->second))[i];
    softCS += occupNb*constant*SoftAndHardXS.first;
    hardCS += occupNb*constant*SoftAndHardXS.second;
  }
  G4double ratio = 0.0;
  
  if (softCS+hardCS) ratio = (hardCS)/(softCS+hardCS);
  return ratio;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::pair<G4double,G4double> 
G4PenelopeIonisationModel::CrossSectionsRatioForElectrons(G4double kineticEnergy,
							  G4double resEnergy,
							  G4double densityCorrection,
							  G4double cutoffEnergy)
{
  std::pair<G4double,G4double> theResult(0.,0.);
  if (kineticEnergy < resEnergy) return theResult;

  //Calculate constants
  G4double gamma = 1.0+kineticEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double cps = kineticEnergy*(kineticEnergy+2.0*electron_mass_c2);
  G4double amol = (kineticEnergy/(kineticEnergy+electron_mass_c2)) * (kineticEnergy/(kineticEnergy+electron_mass_c2)) ;
  G4double hardCont = 0.0;
  G4double softCont = 0.0;
    
  //Distant interactions
  G4double cp1s = (kineticEnergy-resEnergy)*(kineticEnergy-resEnergy+2.0*electron_mass_c2);
  G4double cp1 = std::sqrt(cp1s);
  G4double cp = std::sqrt(cps);
  G4double sdLong=0.0, sdTrans = 0.0, sdDist=0.0;

  //Distant longitudinal interactions
  G4double qm = 0.0;

  if (resEnergy > kineticEnergy*(1e-6))
    {
      qm = std::sqrt((cp-cp1)*(cp-cp1)+(electron_mass_c2*electron_mass_c2))-electron_mass_c2;
    }
  else
    {
      qm = resEnergy*resEnergy/(beta2*2.0*electron_mass_c2);
      qm = qm*(1.0-0.5*qm/electron_mass_c2);
    }

  if (qm < resEnergy)
    {
      sdLong = std::log(resEnergy*(qm+2.0*electron_mass_c2)/(qm*(resEnergy+2.0*electron_mass_c2)));
    }
  else
    {
      sdLong = 0.0;
    }
  
  if (sdLong > 0) {
    sdTrans = std::max(std::log(gamma2)-beta2-densityCorrection,0.0);
    sdDist = sdTrans + sdLong;
    if (cutoffEnergy > resEnergy) 
      {
	softCont = sdDist/resEnergy;
      }
    else
      {
	hardCont = sdDist/resEnergy;
      }
  }

  // Close collisions (Moeller's cross section)
  G4double wl = std::max(cutoffEnergy,resEnergy);
  G4double wu = 0.5*kineticEnergy;
 
  if (wl < (wu-1*eV)) 
    { 
      hardCont += (1.0/(kineticEnergy-wu))-(1.0/(kineticEnergy-wl))
	- (1.0/wu)+(1.0/wl)
	+ (1.0-amol)*std::log(((kineticEnergy-wu)*wl)/((kineticEnergy-wl)*wu))/kineticEnergy
	+ amol*(wu-wl)/(kineticEnergy*kineticEnergy);
      wu=wl;
    }

  wl = resEnergy;
  if (wl > (wu-1*eV)) 
    {
      theResult.first = softCont;
      theResult.second = hardCont;
      return theResult;
    }
  softCont += (1.0/(kineticEnergy-wu))-(1.0/(kineticEnergy-wl))
	- (1.0/wu)+(1.0/wl)
	+ (1.0-amol)*std::log(((kineticEnergy-wu)*wl)/((kineticEnergy-wl)*wu))/kineticEnergy
	+ amol*(wu-wl)/(kineticEnergy*kineticEnergy);
  theResult.first = softCont;
  theResult.second = hardCont;
  return theResult;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::pair<G4double,G4double> 
G4PenelopeIonisationModel::CrossSectionsRatioForPositrons(G4double kineticEnergy,
							  G4double resEnergy,
							  G4double densityCorrection,
							  G4double cutoffEnergy)
{

  std::pair<G4double,G4double> theResult(0.,0.);
  if (kineticEnergy < resEnergy) return theResult;

  //Calculate constants
  G4double gamma = 1.0+kineticEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double cps = kineticEnergy*(kineticEnergy+2.0*electron_mass_c2);
  G4double amol = (kineticEnergy/(kineticEnergy+electron_mass_c2)) * (kineticEnergy/(kineticEnergy+electron_mass_c2)) ;
  G4double help = (gamma+1.0)*(gamma+1.0);
  G4double bha1 = amol*(2.0*help-1.0)/(gamma2-1.0);
  G4double bha2 = amol*(3.0+1.0/help);
  G4double bha3 = amol*2.0*gamma*(gamma-1.0)/help;
  G4double bha4 = amol*(gamma-1.0)*(gamma-1.0)/help;
  G4double hardCont = 0.0;
  G4double softCont = 0.0;

  //Distant interactions
  G4double cp1s = (kineticEnergy-resEnergy)*(kineticEnergy-resEnergy+2.0*electron_mass_c2);
  G4double cp1 = std::sqrt(cp1s);
  G4double cp = std::sqrt(cps);
  G4double sdLong=0.0, sdTrans = 0.0, sdDist=0.0;

  //Distant longitudinal interactions
  G4double qm = 0.0;

  if (resEnergy > kineticEnergy*(1e-6))
    {
      qm = std::sqrt((cp-cp1)*(cp-cp1)+(electron_mass_c2*electron_mass_c2))-electron_mass_c2;
    }
  else
    {
      qm = resEnergy*resEnergy/(beta2*2.0*electron_mass_c2);
      qm = qm*(1.0-0.5*qm/electron_mass_c2);
    }

  if (qm < resEnergy)
    {
      sdLong = std::log(resEnergy*(qm+2.0*electron_mass_c2)/(qm*(resEnergy+2.0*electron_mass_c2)));
    }
  else
    {
      sdLong = 0.0;
    }
  
  if (sdLong > 0) {
    sdTrans = std::max(std::log(gamma2)-beta2-densityCorrection,0.0);
    sdDist = sdTrans + sdLong;
    if (cutoffEnergy > resEnergy) 
      {
	softCont = sdDist/resEnergy;
      }
    else
      {
	hardCont = sdDist/resEnergy;
      }
  }


  // Close collisions (Bhabha's cross section)
  G4double wl = std::max(cutoffEnergy,resEnergy);
  G4double wu = kineticEnergy;
 
  if (wl < (wu-1*eV)) {
    hardCont += (1.0/wl)-(1.0/wu)-bha1*std::log(wu/wl)/kineticEnergy
      + bha2*(wu-wl)/(kineticEnergy*kineticEnergy) 
      -bha3*((wu*wu)-(wl*wl))/(2.0*kineticEnergy*kineticEnergy*kineticEnergy)
      + bha4*((wu*wu*wu)-(wl*wl*wl))/(3.0*kineticEnergy*kineticEnergy*kineticEnergy*kineticEnergy);
    wu=wl;
  }
  wl = resEnergy;
  if (wl > (wu-1*eV)) 
    {  
      theResult.first = softCont;
      theResult.second = hardCont;
      return theResult;
    }
  softCont += (1.0/wl)-(1.0/wu)-bha1*std::log(wu/wl)/kineticEnergy
    + bha2*(wu-wl)/(kineticEnergy*kineticEnergy) 
    -bha3*((wu*wu)-(wl*wl))/(2.0*kineticEnergy*kineticEnergy*kineticEnergy)
    + bha4*((wu*wu*wu)-(wl*wl*wl))/(3.0*kineticEnergy*kineticEnergy*kineticEnergy*kineticEnergy);
  

  theResult.first = softCont;
  theResult.second = hardCont;
  return theResult;

}

G4double G4PenelopeIonisationModel::ComputeStoppingPowerForElectrons(G4double kinEnergy,
                                            		             G4double cutEnergy,
					                             G4double deltaFermi,
					                             G4double resEnergy)
{
   //Calculate constants
  G4double gamma = 1.0+kinEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double cps = kinEnergy*(kinEnergy+2.0*electron_mass_c2);
  G4double amol = (gamma-1.0)*(gamma-1.0)/gamma2;
  G4double sPower = 0.0;
  if (kinEnergy < resEnergy) return sPower;

  //Distant interactions
  G4double cp1s = (kinEnergy-resEnergy)*(kinEnergy-resEnergy+2.0*electron_mass_c2);
  G4double cp1 = std::sqrt(cp1s);
  G4double cp = std::sqrt(cps);
  G4double sdLong=0.0, sdTrans = 0.0, sdDist=0.0;

  //Distant longitudinal interactions
  G4double qm = 0.0;

  if (resEnergy > kinEnergy*(1e-6))
    {
      qm = std::sqrt((cp-cp1)*(cp-cp1)+(electron_mass_c2*electron_mass_c2))-electron_mass_c2;
    }
  else
    {
      qm = resEnergy*resEnergy/(beta2*2.0*electron_mass_c2);
      qm = qm*(1.0-0.5*qm/electron_mass_c2);
    }

  if (qm < resEnergy)
    sdLong = std::log(resEnergy*(qm+2.0*electron_mass_c2)/(qm*(resEnergy+2.0*electron_mass_c2)));
  else
    sdLong = 0.0;
  
  if (sdLong > 0) {
    sdTrans = std::max(std::log(gamma2)-beta2-deltaFermi,0.0);
    sdDist = sdTrans + sdLong;
    if (cutEnergy > resEnergy) sPower = sdDist;
  }


  // Close collisions (Moeller's cross section)
  G4double wl = std::max(cutEnergy,resEnergy);
  G4double wu = 0.5*kinEnergy;
 
  if (wl < (wu-1*eV)) wu=wl;
  wl = resEnergy;
  if (wl > (wu-1*eV)) return sPower;
  sPower += std::log(wu/wl)+(kinEnergy/(kinEnergy-wu))-(kinEnergy/(kinEnergy-wl))
    + (2.0 - amol)*std::log((kinEnergy-wu)/(kinEnergy-wl))
    + amol*((wu*wu)-(wl*wl))/(2.0*kinEnergy*kinEnergy);
  return sPower;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeIonisationModel::ComputeStoppingPowerForPositrons(G4double kinEnergy,
                                            		             G4double cutEnergy,
					                             G4double deltaFermi,
					                             G4double resEnergy)
{
 //Calculate constants
  G4double gamma = 1.0+kinEnergy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta2 = (gamma2-1.0)/gamma2;
  G4double cps = kinEnergy*(kinEnergy+2.0*electron_mass_c2);
  G4double amol = (kinEnergy/(kinEnergy+electron_mass_c2)) * (kinEnergy/(kinEnergy+electron_mass_c2));
  G4double help = (gamma+1.0)*(gamma+1.0);
  G4double bha1 = amol*(2.0*help-1.0)/(gamma2-1.0);
  G4double bha2 = amol*(3.0+1.0/help);
  G4double bha3 = amol*2.0*gamma*(gamma-1.0)/help;
  G4double bha4 = amol*(gamma-1.0)*(gamma-1.0)/help;

  G4double sPower = 0.0;
  if (kinEnergy < resEnergy) return sPower;

  //Distant interactions
  G4double cp1s = (kinEnergy-resEnergy)*(kinEnergy-resEnergy+2.0*electron_mass_c2);
  G4double cp1 = std::sqrt(cp1s);
  G4double cp = std::sqrt(cps);
  G4double sdLong=0.0, sdTrans = 0.0, sdDist=0.0;

  //Distant longitudinal interactions
  G4double qm = 0.0;

  if (resEnergy > kinEnergy*(1e-6))
    {
      qm = std::sqrt((cp-cp1)*(cp-cp1)+(electron_mass_c2*electron_mass_c2))-electron_mass_c2;
    }
  else
    {
      qm = resEnergy*resEnergy/(beta2*2.0*electron_mass_c2);
      qm = qm*(1.0-0.5*qm/electron_mass_c2);
    }

  if (qm < resEnergy)
    sdLong = std::log(resEnergy*(qm+2.0*electron_mass_c2)/(qm*(resEnergy+2.0*electron_mass_c2)));
  else
    sdLong = 0.0;
    
  if (sdLong > 0) {
    sdTrans = std::max(std::log(gamma2)-beta2-deltaFermi,0.0);
    sdDist = sdTrans + sdLong;
    if (cutEnergy > resEnergy) sPower = sdDist;
  }


  // Close collisions (Bhabha's cross section)
  G4double wl = std::max(cutEnergy,resEnergy);
  G4double wu = kinEnergy;
 
  if (wl < (wu-1*eV)) wu=wl;
  wl = resEnergy;
  if (wl > (wu-1*eV)) return sPower;
  sPower += std::log(wu/wl)-bha1*(wu-wl)/kinEnergy
    + bha2*((wu*wu)-(wl*wl))/(2.0*kinEnergy*kinEnergy)
    - bha3*((wu*wu*wu)-(wl*wl*wl))/(3.0*kinEnergy*kinEnergy*kinEnergy)
    + bha4*((wu*wu*wu*wu)-(wl*wl*wl*wl))/(4.0*kinEnergy*kinEnergy*kinEnergy*kinEnergy);
  return sPower;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

/* Notice: the methods here above are only temporary. They will become obsolete in a while */

#include "G4VDataSetAlgorithm.hh"
#include "G4LinLogLogInterpolation.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4EMDataSet.hh"

std::vector<G4VEMDataSet*>* 
G4PenelopeIonisationModel::BuildCrossSectionTable(const G4ParticleDefinition* theParticle)
{
  std::vector<G4VEMDataSet*>* set = new std::vector<G4VEMDataSet*>;

  size_t nOfBins = 200;
  //Temporary vector, a quick way to produce a log-spaced energy grid
  G4PhysicsLogVector* theLogVector = new G4PhysicsLogVector(LowEnergyLimit(),
							    HighEnergyLimit(), 
							    nOfBins);
  G4DataVector* energies;
  G4DataVector* cs;
  
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();


  for (size_t m=0; m<numOfCouples; m++) 
    {
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
      const G4Material* material= couple->GetMaterial();
      const G4ElementVector* elementVector = material->GetElementVector();
      const G4double* nAtomsPerVolume = material->GetAtomicNumDensityVector();
      G4double electronVolumeDensity = 
	material->GetTotNbOfElectPerVolume();  //electron density
      G4int nElements = material->GetNumberOfElements();
  
      G4double tcut  = (*(theCoupleTable->GetEnergyCutsVector(1)))[couple->GetIndex()];

      G4VDataSetAlgorithm* algo =  new G4LinLogLogInterpolation();

      G4VEMDataSet* setForMat = new G4CompositeEMDataSet(algo,1.,1.);

      for (G4int i=0; i<nElements; i++) 
	{
	  G4int iZ = (G4int) (*elementVector)[i]->GetZ();
	  energies = new G4DataVector;
	  cs = new G4DataVector;
	  G4double density = nAtomsPerVolume[i];	  
	  for (size_t bin=0; bin<nOfBins; bin++) 
	    {
	      G4double e = theLogVector->GetLowEdgeEnergy(bin);
	      energies->push_back(e);
	      G4double value = 0.0;	      
	      if(e > tcut) 
		{
		  G4double ratio = CalculateCrossSectionsRatio(e,tcut,iZ,
							       electronVolumeDensity,
							       theParticle);
		  value = crossSectionHandler->FindValue(iZ,e)*ratio*density;
		}
	      cs->push_back(value);
	    }
      G4VDataSetAlgorithm* algo1 = algo->Clone();
      G4VEMDataSet* elSet = new G4EMDataSet(i,energies,cs,algo1,1.,1.);
      setForMat->AddComponent(elSet);
    }
    set->push_back(setForMat);
  }
  delete theLogVector;
  return set;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

G4int G4PenelopeIonisationModel::SampleRandomAtom(const G4MaterialCutsCouple* couple,
                                                     G4double e) const
{
  // Select randomly an element within the material, according to the weight
  // determined by the cross sections in the data set

  const G4Material* material = couple->GetMaterial();
  G4int nElements = material->GetNumberOfElements();

  // Special case: the material consists of one element
  if (nElements == 1)
    {
      G4int Z = (G4int) material->GetZ();
      return Z;
    }

  // Composite material
  const G4ElementVector* elementVector = material->GetElementVector();
  size_t materialIndex = couple->GetIndex();

  G4VEMDataSet* materialSet = (*theXSTable)[materialIndex];
  G4double materialCrossSection0 = 0.0;
  G4DataVector cross;
  cross.clear();
  for ( G4int i=0; i < nElements; i++ )
    {
      G4double cr = materialSet->GetComponent(i)->FindValue(e);
      materialCrossSection0 += cr;
      cross.push_back(materialCrossSection0);
    }

  G4double random = G4UniformRand() * materialCrossSection0;

  for (G4int k=0 ; k < nElements ; k++ )
    {
      if (random <= cross[k]) return (G4int) (*elementVector)[k]->GetZ();
    }
  // It should never get here
  return 0;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeIonisationModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4PenelopeIonisationModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
