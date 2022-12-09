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
// 08 Jan 2010   L Pandola  First implementation
// 01 Feb 2011   L Pandola  Suppress fake energy-violation warning when Auger
//                          is active.
//                          Make sure that fluorescence/Auger is generated
//                          only if above threshold
// 25 May 2011   L Pandola  Renamed (make v2008 as default Penelope)
// 10 Jun 2011   L Pandola  Migrate atomic deexcitation interface
// 07 Oct 2011   L Pandola  Bug fix (potential violation of energy conservation)
// 27 Sep 2013   L Pandola  Migrate to MT paradigm, only master model manages
//                          tables.
// 02 Oct 2013   L Pandola  Rewrite sampling algorithm of SelectRandomShell()
//                          to improve CPU performances
//

#include "G4PenelopePhotoElectricModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4AutoLock.hh"
#include "G4LossTableManager.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4int G4PenelopePhotoElectricModel::fMaxZ;
G4PhysicsTable* G4PenelopePhotoElectricModel::fLogAtomicShellXS[] = {nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopePhotoElectricModel::G4PenelopePhotoElectricModel(const G4ParticleDefinition* part,
							   const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),fParticle(nullptr),
   fAtomDeexcitation(nullptr),fIsInitialised(false),fLocalTable(false)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //

  if (part)
    SetParticle(part);

  fVerboseLevel= 0;
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

G4PenelopePhotoElectricModel::~G4PenelopePhotoElectricModel()
{
  if (IsMaster() || fLocalTable)
    {
      for(G4int i=0; i<=fMaxZ; ++i) 
	{
	  if(fLogAtomicShellXS[i]) { 
	    fLogAtomicShellXS[i]->clearAndDestroy();
	    delete fLogAtomicShellXS[i];
	    fLogAtomicShellXS[i] = nullptr;
	  }
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::Initialise(const G4ParticleDefinition* particle,
					      const G4DataVector& cuts)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopePhotoElectricModel::Initialise()" << G4endl;

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  //Issue warning if the AtomicDeexcitation has not been declared
  if (!fAtomDeexcitation)
    {
      G4cout << G4endl;
      G4cout << "WARNING from G4PenelopePhotoElectricModel " << G4endl;
      G4cout << "Atomic de-excitation module is not instantiated, so there will not be ";
      G4cout << "any fluorescence/Auger emission." << G4endl;
      G4cout << "Please make sure this is intended" << G4endl;
    }

  SetParticle(particle);

  //Only the master model creates/fills/destroys the tables
  if (IsMaster() && particle == fParticle)
    {
      G4ProductionCutsTable* theCoupleTable =
	G4ProductionCutsTable::GetProductionCutsTable();

      for (G4int i=0;i<(G4int)theCoupleTable->GetTableSize();++i)
	{
	  const G4Material* material =
	    theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
	  const G4ElementVector* theElementVector = material->GetElementVector();

	  for (std::size_t j=0;j<material->GetNumberOfElements();++j)
	    {
	      G4int iZ = theElementVector->at(j)->GetZasInt();
	      //read data files only in the master
	      if (!fLogAtomicShellXS[iZ])
		ReadDataFile(iZ);
	    }
	}

      InitialiseElementSelectors(particle,cuts);

      if (fVerboseLevel > 0) {
	G4cout << "Penelope Photo-Electric model v2008 is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / MeV << " MeV - "
	       << HighEnergyLimit() / GeV << " GeV";
      }
    }

  if(fIsInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  fIsInitialised = true;

}

void G4PenelopePhotoElectricModel::InitialiseLocal(const G4ParticleDefinition* part,
						     G4VEmModel *masterModel)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopePhotoElectricModel::InitialiseLocal()" << G4endl;
  //
  //Check that particle matches: one might have multiple master models (e.g.
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      SetElementSelectors(masterModel->GetElementSelectors());

      //Get the const table pointers from the master to the workers
      const G4PenelopePhotoElectricModel* theModel =
	static_cast<G4PenelopePhotoElectricModel*> (masterModel);
      for(G4int i=0; i<=fMaxZ; ++i) 
	fLogAtomicShellXS[i] = theModel->fLogAtomicShellXS[i];
      //Same verbosity for all workers, as the master
      fVerboseLevel = theModel->fVerboseLevel;
    }

 return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
namespace { G4Mutex  PenelopePhotoElectricModelMutex = G4MUTEX_INITIALIZER; }
G4double G4PenelopePhotoElectricModel::ComputeCrossSectionPerAtom(
								  const G4ParticleDefinition*,
								  G4double energy,
								  G4double Z, G4double,
								  G4double, G4double)
{
  //
  // Penelope model v2008
  //
  if (fVerboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopePhotoElectricModel" << G4endl;

  G4int iZ = G4int(Z);

  if (!fLogAtomicShellXS[iZ])
    {
      //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
      //not filled up. This can happen in a UnitTest or via G4EmCalculator
      if (fVerboseLevel > 0)
	{
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to retrieve the shell cross section table for Z=" << iZ << G4endl;
	  ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	  G4Exception("G4PenelopePhotoElectricModel::ComputeCrossSectionPerAtom()",
		      "em2038",JustWarning,ed);
	}
      //protect file reading via autolock
      G4AutoLock lock(&PenelopePhotoElectricModelMutex);
      ReadDataFile(iZ);
      lock.unlock();
    }

  G4double cross = 0;
  G4PhysicsTable* theTable =  fLogAtomicShellXS[iZ];
  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[0];

   if (!totalXSLog)
     {
       G4Exception("G4PenelopePhotoElectricModel::ComputeCrossSectionPerAtom()",
		   "em2039",FatalException,
		   "Unable to retrieve the total cross section table");
       return 0;
     }
   G4double logene = G4Log(energy);
   G4double logXS = totalXSLog->Value(logene);
   cross = G4Exp(logXS);

  if (fVerboseLevel > 2)
    G4cout << "Photoelectric cross section at " << energy/MeV << " MeV for Z=" << Z <<
      " = " << cross/barn << " barn" << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						     const G4MaterialCutsCouple* couple,
						     const G4DynamicParticle* aDynamicGamma,
						     G4double,
						     G4double)
{
  //
  // Photoelectric effect, Penelope model v2008
  //
  // The target atom and the target shell are sampled according to the Livermore
  // database
  //  D.E. Cullen et al., Report UCRL-50400 (1989)
  // The angular distribution of the electron in the final state is sampled
  // according to the Sauter distribution from
  //  F. Sauter, Ann. Phys. 11 (1931) 454
  // The energy of the final electron is given by the initial photon energy minus
  // the binding energy. Fluorescence de-excitation is subsequently produced
  // (to fill the vacancy) according to the general Geant4 G4DeexcitationManager:
  //  J. Stepanek, Comp. Phys. Comm. 1206 pp 1-1-9 (1997)

  if (fVerboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopePhotoElectricModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();

  // always kill primary
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->SetProposedKineticEnergy(0.);

  if (photonEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return ;
    }

  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  if (fVerboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;

  // atom can be selected efficiently if element selectors are initialised
  const G4Element* anElement =
    SelectRandomAtom(couple,G4Gamma::GammaDefinition(),photonEnergy);
  G4int Z = anElement->GetZasInt();
  if (fVerboseLevel > 2)
    G4cout << "Selected " << anElement->GetName() << G4endl;

  // Select the ionised shell in the current atom according to shell cross sections
  //shellIndex = 0 --> K shell
  //             1-3 --> L shells
  //             4-8 --> M shells
  //             9 --> outer shells cumulatively
  //
  std::size_t shellIndex = SelectRandomShell(Z,photonEnergy);

  if (fVerboseLevel > 2)
    G4cout << "Selected shell " << shellIndex << " of element " << anElement->GetName() << G4endl;

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();

  //The number of shell cross section possibly reported in the Penelope database
  //might be different from the number of shells in the G4AtomicTransitionManager
  //(namely, Penelope may contain more shell, especially for very light elements).
  //In order to avoid a warning message from the G4AtomicTransitionManager, I
  //add this protection. Results are anyway changed, because when G4AtomicTransitionManager
  //has a shellID>maxID, it sets the shellID to the last valid shell.
  std::size_t numberOfShells = (std::size_t) transitionManager->NumberOfShells(Z);
  if (shellIndex >= numberOfShells)
    shellIndex = numberOfShells-1;

  const G4AtomicShell* shell = fTransitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();

  //Penelope considers only K, L and M shells. Cross sections of outer shells are
  //not included in the Penelope database. If SelectRandomShell() returns
  //shellIndex = 9, it means that an outer shell was ionized. In this case the
  //Penelope recipe is to set bindingEnergy = 0 (the energy is entirely assigned
  //to the electron) and to disregard fluorescence.
  if (shellIndex == 9)
    bindingEnergy = 0.*eV;

  G4double localEnergyDeposit = 0.0;
  G4double cosTheta = 1.0;

  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;

  // There may be cases where the binding energy of the selected shell is > photon energy
  // In such cases do not generate secondaries
  if (eKineticEnergy > 0.)
    {
      // The electron is created
      // Direction sampled from the Sauter distribution
      cosTheta = SampleElectronDirection(eKineticEnergy);
      G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
      G4double phi = twopi * G4UniformRand() ;
      G4double dirx = sinTheta * std::cos(phi);
      G4double diry = sinTheta * std::sin(phi);
      G4double dirz = cosTheta ;
      G4ThreeVector electronDirection(dirx,diry,dirz); //electron direction
      electronDirection.rotateUz(photonDirection);
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   electronDirection,
							   eKineticEnergy);
      fvect->push_back(electron);
    }
  else
    bindingEnergy = photonEnergy;

  G4double energyInFluorescence = 0; //testing purposes
  G4double energyInAuger = 0; //testing purposes

  //Now, take care of fluorescence, if required. According to the Penelope
  //recipe, I have to skip fluoresence completely if shellIndex == 9
  //(= sampling of a shell outer than K,L,M)
  if (fAtomDeexcitation && shellIndex<9)
    {
      G4int index = couple->GetIndex();
      if (fAtomDeexcitation->CheckDeexcitationActiveRegion(index))
	{
	  std::size_t nBefore = fvect->size();
	  fAtomDeexcitation->GenerateParticles(fvect,shell,Z,index);
	  std::size_t nAfter = fvect->size();

	  if (nAfter > nBefore) //actual production of fluorescence
	    {
	      for (std::size_t j=nBefore;j<nAfter;++j) //loop on products
		{
		  G4double itsEnergy = ((*fvect)[j])->GetKineticEnergy();
		  if (itsEnergy < bindingEnergy) // valid secondary, generate it
		    {
		      bindingEnergy -= itsEnergy;
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

  //Residual energy is deposited locally
  localEnergyDeposit += bindingEnergy;

  if (localEnergyDeposit < 0) //Should not be: issue a G4Exception (warning)
    {
      G4Exception("G4PenelopePhotoElectricModel::SampleSecondaries()",
		  "em2099",JustWarning,"WARNING: Negative local energy deposit");
      localEnergyDeposit = 0;
    }

  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (fVerboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopePhotoElectric" << G4endl;
      G4cout << "Selected shell: " << WriteTargetShell(shellIndex) << " of element " <<
	anElement->GetName() << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      if (eKineticEnergy)
	G4cout << "Outgoing electron " << eKineticEnergy/keV << " keV" << G4endl;
      if (energyInFluorescence)
	G4cout << "Fluorescence x-rays: " << energyInFluorescence/keV << " keV" << G4endl;
      if (energyInAuger)
	G4cout << "Auger electrons: " << energyInAuger/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " <<
	(eKineticEnergy+energyInFluorescence+localEnergyDeposit+energyInAuger)/keV <<
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (fVerboseLevel > 0)
    {
      G4double energyDiff =
	std::fabs(eKineticEnergy+energyInFluorescence+localEnergyDeposit+energyInAuger-photonEnergy);
      if (energyDiff > 0.05*keV)
	{
	  G4cout << "Warning from G4PenelopePhotoElectric: problem with energy conservation: " <<
	    (eKineticEnergy+energyInFluorescence+localEnergyDeposit+energyInAuger)/keV
		 << " keV (final) vs. " <<
	    photonEnergy/keV << " keV (initial)" << G4endl;
	  G4cout << "-----------------------------------------------------------" << G4endl;
	  G4cout << "Energy balance from G4PenelopePhotoElectric" << G4endl;
	  G4cout << "Selected shell: " << WriteTargetShell(shellIndex) << " of element " <<
	    anElement->GetName() << G4endl;
	  G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
	  G4cout << "-----------------------------------------------------------" << G4endl;
	  if (eKineticEnergy)
	    G4cout << "Outgoing electron " << eKineticEnergy/keV << " keV" << G4endl;
	  if (energyInFluorescence)
	    G4cout << "Fluorescence x-rays: " << energyInFluorescence/keV << " keV" << G4endl;
	  if (energyInAuger)
	    G4cout << "Auger electrons: " << energyInAuger/keV << " keV" << G4endl;
	  G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
	  G4cout << "Total final state: " <<
	    (eKineticEnergy+energyInFluorescence+localEnergyDeposit+energyInAuger)/keV <<
	    " keV" << G4endl;
	  G4cout << "-----------------------------------------------------------" << G4endl;
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopePhotoElectricModel::SampleElectronDirection(G4double energy)
{
  G4double costheta = 1.0;
  if (energy>1*GeV) return costheta;

  //1) initialize energy-dependent variables
  // Variable naming according to Eq. (2.24) of Penelope Manual
  // (pag. 44)
  G4double gamma = 1.0 + energy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta = std::sqrt((gamma2-1.0)/gamma2);

  // ac corresponds to "A" of Eq. (2.31)
  //
  G4double ac = (1.0/beta) - 1.0;
  G4double a1 = 0.5*beta*gamma*(gamma-1.0)*(gamma-2.0);
  G4double a2 = ac + 2.0;
  G4double gtmax = 2.0*(a1 + 1.0/ac);

  G4double tsam = 0;
  G4double gtr = 0;

  //2) sampling. Eq. (2.31) of Penelope Manual
  // tsam = 1-std::cos(theta)
  // gtr = rejection function according to Eq. (2.28)
  do{
    G4double rand = G4UniformRand();
    tsam = 2.0*ac * (2.0*rand + a2*std::sqrt(rand)) / (a2*a2 - 4.0*rand);
    gtr = (2.0 - tsam) * (a1 + 1.0/(ac+tsam));
  }while(G4UniformRand()*gtmax > gtr);
  costheta = 1.0-tsam;

  return costheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::ReadDataFile(G4int Z)
{
  if (!IsMaster())
      //Should not be here!
    G4Exception("G4PenelopePhotoElectricModel::ReadDataFile()",
		"em0100",FatalException,"Worker thread in this method");

  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopePhotoElectricModel::ReadDataFile()" << G4endl;
      G4cout << "Going to read PhotoElectric data files for Z=" << Z << G4endl;
    }

    const char* path = G4FindDataDir("G4LEDATA");
    if(!path)
    {
      G4String excep = "G4PenelopePhotoElectricModel - G4LEDATA environment variable not set!";
      G4Exception("G4PenelopePhotoElectricModel::ReadDataFile()",
		  "em0006",FatalException,excep);
      return;
    }

  /*
    Read the cross section file
  */
  std::ostringstream ost;
  if (Z>9)
    ost << path << "/penelope/photoelectric/pdgph" << Z << ".p08";
  else
    ost << path << "/penelope/photoelectric/pdgph0" << Z << ".p08";
  std::ifstream file(ost.str().c_str());
  if (!file.is_open())
    {
      G4String excep = "G4PenelopePhotoElectricModel - data file " + G4String(ost.str()) + " not found!";
      G4Exception("G4PenelopePhotoElectricModel::ReadDataFile()",
		  "em0003",FatalException,excep);
    }
  //I have to know in advance how many points are in the data list
  //to initialize the G4PhysicsFreeVector()
  std::size_t ndata=0;
  G4String line;
  while( getline(file, line) )
    ndata++;
  ndata -= 1;
  //G4cout << "Found: " << ndata << " lines" << G4endl;

  file.clear();
  file.close();
  file.open(ost.str().c_str());

  G4int readZ =0;
  std::size_t nShells= 0;
  file >> readZ >> nShells;

  if (fVerboseLevel > 3)
    G4cout << "Element Z=" << Z << " , nShells = " << nShells << G4endl;

  //check the right file is opened.
  if (readZ != Z || nShells <= 0 || nShells > 50) //protect nShell against large values
    {
      G4ExceptionDescription ed;
      ed << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception("G4PenelopePhotoElectricModel::ReadDataFile()",
		  "em0005",FatalException,ed);
      return;
    }
  G4PhysicsTable* thePhysicsTable = new G4PhysicsTable();
  
  //the table has to contain nShell+1 G4PhysicsFreeVectors,
  //(theTable)[0] --> total cross section
  //(theTable)[ishell] --> cross section for shell (ishell-1)

  //reserve space for the vectors
  //everything is log-log
  for (std::size_t i=0;i<nShells+1;++i)
    thePhysicsTable->push_back(new G4PhysicsFreeVector(ndata));

  std::size_t k =0;
  for (k=0;k<ndata && !file.eof();++k)
    {
      G4double energy = 0;
      G4double aValue = 0;
      file >> energy ;
      energy *= eV;
      G4double logene = G4Log(energy);
      //loop on the columns
      for (std::size_t i=0;i<nShells+1;++i)
	{
	  file >> aValue;
	  aValue *= barn;
	  G4PhysicsFreeVector* theVec = (G4PhysicsFreeVector*) ((*thePhysicsTable)[i]);
	  if (aValue < 1e-40*cm2) //protection against log(0)
	    aValue = 1e-40*cm2;
	  theVec->PutValue(k,logene,G4Log(aValue));
	}
    }

  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopePhotoElectricModel: read " << k << " points for element Z = "
	     << Z << G4endl;
    }

  fLogAtomicShellXS[Z] = thePhysicsTable;

  file.close();
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4PenelopePhotoElectricModel::GetNumberOfShellXS(G4int Z)
{
  if (!IsMaster())
    //Should not be here!
    G4Exception("G4PenelopePhotoElectricModel::GetNumberOfShellXS()",
		"em0100",FatalException,"Worker thread in this method");

  //read data files
  if (!fLogAtomicShellXS[Z])
    ReadDataFile(Z);
  //now it should be ok
  if (!fLogAtomicShellXS[Z])
     {
       G4ExceptionDescription ed;
       ed << "Cannot find shell cross section data for Z=" << Z << G4endl;
       G4Exception("G4PenelopePhotoElectricModel::GetNumberOfShellXS()",
		   "em2038",FatalException,ed);
     }
  //one vector is allocated for the _total_ cross section
  std::size_t nEntries = fLogAtomicShellXS[Z]->entries();
  return  (nEntries-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopePhotoElectricModel::GetShellCrossSection(G4int Z,std::size_t shellID,G4double energy)
{
  //this forces also the loading of the data
  std::size_t entries = GetNumberOfShellXS(Z);

  if (shellID >= entries)
    {
      G4cout << "Element Z=" << Z << " has data for " << entries << " shells only" << G4endl;
      G4cout << "so shellID should be from 0 to " << entries-1 << G4endl;
      return 0;
    }

  G4PhysicsTable* theTable =  fLogAtomicShellXS[Z];
  //[0] is the total XS, shellID is in the element [shellID+1]
  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[shellID+1];

  if (!totalXSLog)
     {
       G4Exception("G4PenelopePhotoElectricModel::GetShellCrossSection()",
		   "em2039",FatalException,
		   "Unable to retrieve the total cross section table");
       return 0;
     }
   G4double logene = G4Log(energy);
   G4double logXS = totalXSLog->Value(logene);
   G4double cross = G4Exp(logXS);
   if (cross < 2e-40*cm2) cross = 0;
   return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4PenelopePhotoElectricModel::WriteTargetShell(std::size_t shellID)
{
  G4String theShell = "outer shell";
  if (shellID == 0)
    theShell = "K";
  else if (shellID == 1)
    theShell = "L1";
  else if (shellID == 2)
    theShell = "L2";
  else if (shellID == 3)
    theShell = "L3";
  else if (shellID == 4)
    theShell = "M1";
  else if (shellID == 5)
    theShell = "M2";
  else if (shellID == 6)
    theShell = "M3";
  else if (shellID == 7)
    theShell = "M4";
  else if (shellID == 8)
    theShell = "M5";

  return theShell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopePhotoElectricModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

std::size_t G4PenelopePhotoElectricModel::SelectRandomShell(G4int Z,G4double energy)
{
  G4double logEnergy = G4Log(energy);

  //Check if data have been read (it should be!)
  if (!fLogAtomicShellXS[Z])
     {
       G4ExceptionDescription ed;
       ed << "Cannot find shell cross section data for Z=" << Z << G4endl;
       G4Exception("G4PenelopePhotoElectricModel::SelectRandomShell()",
		   "em2038",FatalException,ed);
     }

  G4PhysicsTable* theTable =  fLogAtomicShellXS[Z];

  G4double sum = 0;
  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[0];
  G4double logXS = totalXSLog->Value(logEnergy);
  G4double totalXS = G4Exp(logXS);

  //Notice: totalXS is the total cross section and it does *not* correspond to
  //the sum of partialXS's, since these include only K, L and M shells.
  //
  // Therefore, here one have to consider the possibility of ionisation of
  // an outer shell. Conventionally, it is indicated with id=10 in Penelope
  //
  G4double random = G4UniformRand()*totalXS;

  for (std::size_t k=1;k<theTable->entries();++k)
    {
      //Add one shell
      G4PhysicsFreeVector* partialXSLog = (G4PhysicsFreeVector*) (*theTable)[k];
      G4double logXSLocal = partialXSLog->Value(logEnergy);
      G4double partialXS = G4Exp(logXSLocal);
      sum += partialXS;
      if (random <= sum)
	return k-1;
    }
  //none of the shells K, L, M: return outer shell
  return 9;
}
