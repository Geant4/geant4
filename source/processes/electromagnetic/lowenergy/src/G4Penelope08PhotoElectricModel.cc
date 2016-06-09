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
// $Id: G4Penelope08PhotoElectricModel.cc,v 1.6 2010-12-15 10:26:41 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 08 Jan 2010   L Pandola  First implementation

#include "G4Penelope08PhotoElectricModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4Penelope08PhotoElectricModel::G4Penelope08PhotoElectricModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),logAtomicShellXS(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //by default the model will inkove the atomic deexcitation
  SetDeexcitationFlag(true);  
  ActivateAuger(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4Penelope08PhotoElectricModel::~G4Penelope08PhotoElectricModel()
{  
  std::map <const G4int,G4PhysicsTable*>::iterator i;
  if (logAtomicShellXS)
    {
      for (i=logAtomicShellXS->begin();i != logAtomicShellXS->end();i++)
	{
	  G4PhysicsTable* tab = i->second;
	  tab->clearAndDestroy();
	  delete tab;
	}
    }
  delete logAtomicShellXS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08PhotoElectricModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling  G4Penelope08PhotoElectricModel::Initialise()" << G4endl;

  // logAtomicShellXS is created only once, since it is  never cleared
  if (!logAtomicShellXS)
    logAtomicShellXS = new std::map<const G4int,G4PhysicsTable*>;

  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 0) { 
    G4cout << "Penelope Photo-Electric model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08PhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  //
  // Penelope model. 
  //

  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4Penelope08PhotoElectricModel" << G4endl;

  G4int iZ = (G4int) Z;

  //read data files
  if (!logAtomicShellXS->count(iZ))
    ReadDataFile(iZ);
  //now it should be ok
  if (!logAtomicShellXS->count(iZ))
     {
       G4cout << "Problem in G4Penelope08PhotoElectricModel::ComputeCrossSectionPerAtom"
              << G4endl;
       G4Exception();
     }

  G4double cross = 0;

  G4PhysicsTable* theTable =  logAtomicShellXS->find(iZ)->second;
  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[0];

   if (!totalXSLog)
     {
       G4cout << "Problem in G4Penelope08PhotoElectricModel::ComputeCrossSectionPerAtom"
         << G4endl;
       G4Exception();
       return 0;
     }
   G4double logene = std::log(energy);
   G4double logXS = totalXSLog->Value(logene);
   cross = std::exp(logXS);
 
  if (verboseLevel > 2)
    G4cout << "Photoelectric cross section at " << energy/MeV << " MeV for Z=" << Z <<
      " = " << cross/barn << " barn" << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08PhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						       const G4MaterialCutsCouple* couple,
						       const G4DynamicParticle* aDynamicGamma,
						       G4double,
						       G4double)
{
  //
  // Photoelectric effect, Penelope model
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

  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4Penelope08PhotoElectricModel" << G4endl;

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
  if (verboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;

  // atom can be selected efficiently if element selectors are initialised
  const G4Element* anElement =
    SelectRandomAtom(couple,G4Gamma::GammaDefinition(),photonEnergy);
  G4int Z = (G4int) anElement->GetZ();
  if (verboseLevel > 2)
    G4cout << "Selected " << anElement->GetName() << G4endl;
  
  // Select the ionised shell in the current atom according to shell cross sections
  //shellIndex = 0 --> K shell
  //             1-3 --> L shells
  //             4-8 --> M shells
  //             9 --> outer shells cumulatively
  //
  size_t shellIndex = SelectRandomShell(Z,photonEnergy);

  if (verboseLevel > 2)
    G4cout << "Selected shell " << shellIndex << " of element " << anElement->GetName() << G4endl;

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();

  //The number of shell cross section possibly reported in the Penelope database 
  //might be different from the number of shells in the G4AtomicTransitionManager
  //(namely, Penelope may contain more shell, especially for very light elements).
  //In order to avoid a warning message from the G4AtomicTransitionManager, I 
  //add this protection. Results are anyway changed, because when G4AtomicTransitionManager
  //has a shellID>maxID, it sets the shellID to the last valid shell. 
  size_t numberOfShells = (size_t) transitionManager->NumberOfShells(Z);
  if (shellIndex >= numberOfShells)
    shellIndex = numberOfShells-1;

  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();

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
    {
      bindingEnergy = photonEnergy;
    }

  G4double energyInFluorescence = 0; //testing purposes

  //Now, take care of fluorescence, if required
  if(DeexcitationFlag() && Z > 5) 
    {
      const G4ProductionCutsTable* theCoupleTable=
	G4ProductionCutsTable::GetProductionCutsTable();
      size_t indx = couple->GetIndex();
      G4double cutG = (*(theCoupleTable->GetEnergyCutsVector(0)))[indx];
      G4double cutE = (*(theCoupleTable->GetEnergyCutsVector(1)))[indx];
      
      // Protection to avoid generating photons in the unphysical case of 
      // shell binding energy > photon energy
      if (bindingEnergy > cutG || bindingEnergy > cutE)
	{
	  deexcitationManager.SetCutForSecondaryPhotons(cutG);
	  deexcitationManager.SetCutForAugerElectrons(cutE);
	  std::vector<G4DynamicParticle*>* photonVector = 
	    deexcitationManager.GenerateParticles(Z,shellId); 
	  //Check for secondaries
          if(photonVector) 
	    {
	      for (size_t k=0; k< photonVector->size(); k++)
		{
		  G4DynamicParticle* aPhoton = (*photonVector)[k];
		  if (aPhoton)
		    {
		      G4double itsEnergy = aPhoton->GetKineticEnergy();
		      if (itsEnergy <= bindingEnergy)
			{
			  if(aPhoton->GetDefinition() == G4Gamma::Gamma())
			    energyInFluorescence += itsEnergy;
			  bindingEnergy -= itsEnergy;
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
  //Residual energy is deposited locally
  localEnergyDeposit += bindingEnergy;
      
  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING - "
	     << "G4Penelope08PhotoElectric::PostStepDoIt - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit = 0;
    }

  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4Penelope08PhotoElectric" << G4endl;
      G4cout << "Selected shell: " << WriteTargetShell(shellIndex) << " of element " << 
	anElement->GetName() << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      if (eKineticEnergy)
	G4cout << "Outgoing electron " << eKineticEnergy/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (eKineticEnergy+energyInFluorescence+localEnergyDeposit)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = 
	std::fabs(eKineticEnergy+energyInFluorescence+localEnergyDeposit-photonEnergy);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4Penelope08PhotoElectric: problem with energy conservation: " << 
	  (eKineticEnergy+energyInFluorescence+localEnergyDeposit)/keV 
	       << " keV (final) vs. " << 
	  photonEnergy/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4Penelope08PhotoElectricModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4Penelope08PhotoElectricModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08PhotoElectricModel::SampleElectronDirection(G4double energy)
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

void G4Penelope08PhotoElectricModel::ReadDataFile(G4int Z)
{
  if (verboseLevel > 2)
    {
      G4cout << "G4Penelope08PhotoElectricModel::ReadDataFile()" << G4endl;
      G4cout << "Going to read PhotoElectric data files for Z=" << Z << G4endl;
    }
 
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4Penelope08PhotoElectricModel - G4LEDATA environment variable not set!";
      G4Exception(excep);
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
      G4String excep = "G4Penelope08PhotoElectricModel - data file " + G4String(ost.str()) + " not found!";
      G4Exception(excep);
    }
  //I have to know in advance how many points are in the data list
  //to initialize the G4PhysicsFreeVector()
  size_t ndata=0;
  G4String line;
  while( getline(file, line) )
    ndata++;
  ndata -= 1;
  //G4cout << "Found: " << ndata << " lines" << G4endl;

  file.clear();
  file.close();
  file.open(ost.str().c_str());

  G4int readZ =0;
  size_t nShells= 0;
  file >> readZ >> nShells;

  if (verboseLevel > 3)
    G4cout << "Element Z=" << Z << " , nShells = " << nShells << G4endl;

  //check the right file is opened.
  if (readZ != Z || nShells <= 0)
    {
      G4cout << "G4Penelope08PhotoElectricModel::ReadDataFile()" << G4endl;
      G4cout << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception();
      return;
    }
  G4PhysicsTable* thePhysicsTable = new G4PhysicsTable();

  //the table has to contain nShell+1 G4PhysicsFreeVectors, 
  //(theTable)[0] --> total cross section
  //(theTable)[ishell] --> cross section for shell (ishell-1)

  //reserve space for the vectors
  //everything is log-log
  for (size_t i=0;i<nShells+1;i++)
    thePhysicsTable->push_back(new G4PhysicsFreeVector(ndata));

  size_t k =0;
  for (k=0;k<ndata && !file.eof();k++)
    {
      G4double energy = 0;
      G4double aValue = 0;
      file >> energy ;
      energy *= eV;
      G4double logene = std::log(energy);
      //loop on the columns
      for (size_t i=0;i<nShells+1;i++)
	{
	  file >> aValue;
	  aValue *= barn;
	  G4PhysicsFreeVector* theVec = (G4PhysicsFreeVector*) ((*thePhysicsTable)[i]);	 
	  if (aValue < 1e-40*cm2) //protection against log(0)
	    aValue = 1e-40*cm2;
	  theVec->PutValue(k,logene,std::log(aValue));
	}
    }

  if (verboseLevel > 2)
    {
      G4cout << "G4Penelope08PhotoElectricModel: read " << k << " points for element Z = " 
	     << Z << G4endl;
    }

  logAtomicShellXS->insert(std::make_pair(Z,thePhysicsTable));
 
  file.close();
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

size_t G4Penelope08PhotoElectricModel::SelectRandomShell(G4int Z,G4double energy)
{
  G4double logEnergy = std::log(energy);

  //Check if data have been read (it should be!)
  if (!logAtomicShellXS->count(Z))
     {
       G4cout << "Problem in G4Penelope08PhotoElectricModel::SelectRandomShell" << G4endl;
       G4cout << "Cannot find data for Z=" << Z << G4endl;
       G4Exception();
     }

  size_t shellIndex = 0;
 
  G4PhysicsTable* theTable =  logAtomicShellXS->find(Z)->second;

  G4DataVector* tempVector = new G4DataVector();

  G4double sum = 0;
  //loop on shell partial XS, retrieve the value for the given energy and store on 
  //a temporary vector
  tempVector->push_back(sum); //first element is zero

  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[0];
  G4double logXS = totalXSLog->Value(logEnergy);
  G4double totalXS = std::exp(logXS);
					   
  //Notice: totalXS is the total cross section and it does *not* correspond to 
  //the sum of partialXS's, since these include only K, L and M shells.
  //
  // Therefore, here one have to consider the possibility of ionisation of 
  // an outer shell. Conventionally, it is indicated with id=10 in Penelope
  //
  
  for (size_t k=1;k<theTable->entries();k++)
    {
      G4PhysicsFreeVector* partialXSLog = (G4PhysicsFreeVector*) (*theTable)[k];
      G4double logXS = partialXSLog->Value(logEnergy);
      G4double partialXS = std::exp(logXS);
      sum += partialXS;
      tempVector->push_back(sum);     
    }

  tempVector->push_back(totalXS); //last element

  G4double random = G4UniformRand()*totalXS; 

  /*
  for (size_t i=0;i<tempVector->size(); i++)
    G4cout << i << " " << (*tempVector)[i]/totalXS << G4endl;
  */
  
  //locate bin of tempVector
  //Now one has to sample according to the elements in tempVector
  //This gives the left edge of the interval...
  size_t lowerBound = 0;
  size_t upperBound = tempVector->size()-1; 
  while (lowerBound <= upperBound)
   {
     size_t midBin = (lowerBound + upperBound)/2;
     if( random < (*tempVector)[midBin])
       upperBound = midBin-1; 
     else
       lowerBound = midBin+1; 
   }
 
  shellIndex = upperBound;

  delete tempVector;
  return shellIndex;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

size_t G4Penelope08PhotoElectricModel::GetNumberOfShellXS(G4int Z)
{
  //read data files
  if (!logAtomicShellXS->count(Z))
    ReadDataFile(Z);
  //now it should be ok
  if (!logAtomicShellXS->count(Z))
     {
       G4cout << "Problem in G4Penelope08PhotoElectricModel::GetNumberOfShellXS()"
              << G4endl;
       G4Exception();
     }
  //one vector is allocated for the _total_ cross section
  size_t nEntries = logAtomicShellXS->find(Z)->second->entries();
  return  (nEntries-1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4Penelope08PhotoElectricModel::GetShellCrossSection(G4int Z,size_t shellID,G4double energy)
{
  //this forces also the loading of the data
  size_t entries = GetNumberOfShellXS(Z);

  if (shellID >= entries)
    {
      G4cout << "Element Z=" << Z << " has data for " << entries << " shells only" << G4endl;
      G4cout << "so shellID should be from 0 to " << entries-1 << G4endl;
      return 0;
    }
  
  G4PhysicsTable* theTable =  logAtomicShellXS->find(Z)->second;
  //[0] is the total XS, shellID is in the element [shellID+1]
  G4PhysicsFreeVector* totalXSLog = (G4PhysicsFreeVector*) (*theTable)[shellID+1];
 
  if (!totalXSLog)
     {
       G4cout << "Problem in G4Penelope08PhotoElectricModel::GetShellCrossSection()"
         << G4endl;
       G4Exception();
       return 0;
     }
   G4double logene = std::log(energy);
   G4double logXS = totalXSLog->Value(logene);
   G4double cross = std::exp(logXS);
   if (cross < 2e-40*cm2) cross = 0;
   return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String G4Penelope08PhotoElectricModel::WriteTargetShell(size_t shellID)
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
