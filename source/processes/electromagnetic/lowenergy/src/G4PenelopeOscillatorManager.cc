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
// Authors: Luciano Pandola (luciano.pandola at lngs.infn.it)
//
// History:
// -----------
//
//  03 Dec 2009  First working version, Luciano Pandola
//  16 Feb 2010  Added methods to store also total Z and A for the
//               molecule, Luciano Pandola
//  19 Feb 2010  Scale the Hartree factors in the Compton Oscillator
//               table by (1/fine_structure_const), since the models use
//               always the ratio (hartreeFactor/fine_structure_const)
//  16 Mar 2010  Added methods to calculate and store mean exc energy
//               and plasma energy (used for Ionisation). L Pandola
//  18 Mar 2010  Added method to retrieve number of atoms per
//               molecule. L. Pandola
//  06 Sep 2011  Override the local Penelope database and use the main
//               G4AtomicDeexcitation database to retrieve the shell
//               binding energies. L. Pandola
//  15 Mar 2012  Added method to retrieve number of atom of given Z per
//               molecule. Restore the original Penelope database for levels
//               below 100 eV. L. Pandola
//
// -------------------------------------------------------------------

#include "G4PenelopeOscillatorManager.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Material.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4PenelopeOscillatorManager* G4PenelopeOscillatorManager::instance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillatorManager::G4PenelopeOscillatorManager() :
  fOscillatorStoreIonisation(nullptr),fOscillatorStoreCompton(nullptr),
  fAtomicNumber(nullptr),fAtomicMass(nullptr),fExcitationEnergy(nullptr),
  fPlasmaSquared(nullptr),fAtomsPerMolecule(nullptr),
  fAtomTablePerMolecule(nullptr)
{
  fReadElementData = false;
  for (G4int i=0;i<5;i++)
    {
      for (G4int j=0;j<2000;j++)
	fElementData[i][j] = 0.;
    }
  fVerbosityLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillatorManager::~G4PenelopeOscillatorManager()
{
  Clear();
  delete instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillatorManager* G4PenelopeOscillatorManager::GetOscillatorManager()
{
  if (!instance)
    instance = new G4PenelopeOscillatorManager();
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeOscillatorManager::Clear()
{
  if (fVerbosityLevel > 1)
    G4cout << " G4PenelopeOscillatorManager::Clear() - Clean Oscillator Tables" << G4endl;

  //Clean up OscillatorStoreIonisation
  for (auto& item : (*fOscillatorStoreIonisation))
    {
      G4PenelopeOscillatorTable* table = item.second;
      if (table)
	{
	  for (std::size_t k=0;k<table->size();++k) //clean individual oscillators
	    {
	      if ((*table)[k])
		delete ((*table)[k]);
	    }
	  delete table;
	}
    }
  delete fOscillatorStoreIonisation;

  //Clean up OscillatorStoreCompton
  for (auto& item : (*fOscillatorStoreCompton))
    {
      G4PenelopeOscillatorTable* table = item.second;
      if (table)
	{
	  for (std::size_t k=0;k<table->size();++k) //clean individual oscillators
	    {
	      if ((*table)[k])
		delete ((*table)[k]);
	    }
	  delete table;
	}
    }
  delete fOscillatorStoreCompton;

  if (fAtomicMass) delete fAtomicMass;
  if (fAtomicNumber) delete fAtomicNumber;
  if (fExcitationEnergy) delete fExcitationEnergy;
  if (fPlasmaSquared) delete fPlasmaSquared;
  if (fAtomsPerMolecule) delete fAtomsPerMolecule;
  if (fAtomTablePerMolecule) delete fAtomTablePerMolecule;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeOscillatorManager::Dump(const G4Material* material)
{
  G4PenelopeOscillatorTable* theTable = GetOscillatorTableIonisation(material);
  if (!theTable)
    {
      G4cout << " G4PenelopeOscillatorManager::Dump " << G4endl;
      G4cout << "Problem in retrieving the Ionisation Oscillator Table for " 
	     << material->GetName() << G4endl;
      return;
    }
  G4cout << "*********************************************************************" << G4endl;
  G4cout << " Penelope Oscillator Table Ionisation for " << material->GetName() << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "The table contains " << theTable->size() << " oscillators " << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  if (theTable->size() < 10)
    for (std::size_t k=0;k<theTable->size();++k)
      {
	G4cout << "Oscillator # " << k << " Z = " << (*theTable)[k]->GetParentZ() <<
	  " Shell Flag = " << (*theTable)[k]->GetShellFlag() <<
	  " Parent shell ID = " << (*theTable)[k]->GetParentShellID() << G4endl;
	G4cout << "Ionisation energy = " << (*theTable)[k]->GetIonisationEnergy()/eV << " eV" << G4endl;
	G4cout << "Occupation number = " << (*theTable)[k]->GetOscillatorStrength() << G4endl;
	G4cout << "Resonance energy = " << (*theTable)[k]->GetResonanceEnergy()/eV << " eV" << G4endl;
	G4cout << "Cufoff resonance energy = " <<
		(*theTable)[k]->GetCutoffRecoilResonantEnergy()/eV << " eV" << G4endl;
	G4cout << "*********************************************************************" << G4endl;
      }
  for (std::size_t k=0;k<theTable->size();++k)
    {
      G4cout << k << " " << (*theTable)[k]->GetOscillatorStrength() << " " <<
	(*theTable)[k]->GetIonisationEnergy()/eV << " " 
	     << (*theTable)[k]->GetResonanceEnergy()/eV << " " <<
	(*theTable)[k]->GetParentZ() << " " << (*theTable)[k]->GetShellFlag() << " " <<
	(*theTable)[k]->GetParentShellID() << G4endl;
    }
  G4cout << "*********************************************************************" << G4endl;

  //Compton table
  theTable = GetOscillatorTableCompton(material);
  if (!theTable)
    {
      G4cout << " G4PenelopeOscillatorManager::Dump " << G4endl;
      G4cout << "Problem in retrieving the Compton Oscillator Table for " << 
	material->GetName() << G4endl;
      return;
    }
  G4cout << "*********************************************************************" << G4endl;
  G4cout << " Penelope Oscillator Table Compton for " << material->GetName() << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "The table contains " << theTable->size() << " oscillators " << G4endl;
  G4cout << "*********************************************************************" << G4endl;
  if (theTable->size() < 10)
    for (std::size_t k=0;k<theTable->size();++k)
      {
	G4cout << "Oscillator # " << k << " Z = " << (*theTable)[k]->GetParentZ() <<
	  " Shell Flag = " << (*theTable)[k]->GetShellFlag() <<
	   " Parent shell ID = " << (*theTable)[k]->GetParentShellID() << G4endl;
	G4cout << "Compton index = " << (*theTable)[k]->GetHartreeFactor() << G4endl;
	G4cout << "Ionisation energy = " << (*theTable)[k]->GetIonisationEnergy()/eV << " eV" << G4endl;
	G4cout << "Occupation number = " << (*theTable)[k]->GetOscillatorStrength() << G4endl;
	G4cout << "*********************************************************************" << G4endl;
      }
  for (std::size_t k=0;k<theTable->size();++k)
    {
      G4cout << k << " " << (*theTable)[k]->GetOscillatorStrength() << " " <<
	(*theTable)[k]->GetIonisationEnergy()/eV << " " << (*theTable)[k]->GetHartreeFactor() << " " <<
	(*theTable)[k]->GetParentZ() << " " << (*theTable)[k]->GetShellFlag() << " " <<
	(*theTable)[k]->GetParentShellID() << G4endl;
    }
  G4cout << "*********************************************************************" << G4endl;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeOscillatorManager::CheckForTablesCreated()
{
  //Tables should be created at the same time, since they are both filled
  //simultaneously
  if (!fOscillatorStoreIonisation)
    {
      fOscillatorStoreIonisation = new std::map<const G4Material*,G4PenelopeOscillatorTable*>;
      if (!fReadElementData)
	ReadElementData();
      if (!fOscillatorStoreIonisation)
	//It should be ok now
	G4Exception("G4PenelopeOscillatorManager::GetOscillatorTableIonisation()",
		    "em2034",FatalException,
		    "Problem in allocating the Oscillator Store for Ionisation");
    }

  if (!fOscillatorStoreCompton)
    {
      fOscillatorStoreCompton = new std::map<const G4Material*,G4PenelopeOscillatorTable*>;
      if (!fReadElementData)
	ReadElementData();
      if (!fOscillatorStoreCompton)
	//It should be ok now
	G4Exception("G4PenelopeOscillatorManager::GetOscillatorTableIonisation()",
		    "em2034",FatalException,
		    "Problem in allocating the Oscillator Store for Compton");
    }

  if (!fAtomicNumber)
    fAtomicNumber = new std::map<const G4Material*,G4double>;
  if (!fAtomicMass)
    fAtomicMass = new std::map<const G4Material*,G4double>;
  if (!fExcitationEnergy)
    fExcitationEnergy = new std::map<const G4Material*,G4double>;
  if (!fPlasmaSquared)
    fPlasmaSquared = new std::map<const G4Material*,G4double>;
  if (!fAtomsPerMolecule)
    fAtomsPerMolecule = new std::map<const G4Material*,G4double>;
  if (!fAtomTablePerMolecule)
    fAtomTablePerMolecule = new std::map< std::pair<const G4Material*,G4int>, G4double>;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeOscillatorManager::GetTotalZ(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fAtomicNumber->count(mat))
    return fAtomicNumber->find(mat)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fAtomicNumber->count(mat))
    return fAtomicNumber->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetTotalZ() " << G4endl;
      G4cout << "Impossible to retrieve the total Z for " << mat->GetName() << G4endl;
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeOscillatorManager::GetTotalA(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fAtomicMass->count(mat))
    return fAtomicMass->find(mat)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fAtomicMass->count(mat))
    return fAtomicMass->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetTotalA() " << G4endl;
      G4cout << "Impossible to retrieve the total A for " << mat->GetName() << G4endl;
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillatorTable* G4PenelopeOscillatorManager::GetOscillatorTableIonisation(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fOscillatorStoreIonisation->count(mat))
    {
      //Ok, it exists
      return fOscillatorStoreIonisation->find(mat)->second;
    }

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fOscillatorStoreIonisation->count(mat))
    return fOscillatorStoreIonisation->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetOscillatorTableIonisation() " << G4endl;
      G4cout << "Impossible to create ionisation oscillator table for " << mat->GetName() << G4endl;
      return nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator* G4PenelopeOscillatorManager::GetOscillatorIonisation(const G4Material* material,
									   G4int index)
{
  G4PenelopeOscillatorTable* theTable = GetOscillatorTableIonisation(material);
  if (((std::size_t)index) < theTable->size())
    return (*theTable)[index];
  else
    {
      G4cout << "WARNING: Ionisation table for material " << material->GetName() << " has " <<
	theTable->size() << " oscillators" << G4endl;
      G4cout << "Oscillator #" << index << " cannot be retrieved" << G4endl;
      G4cout << "Returning null pointer" << G4endl;
      return nullptr;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillatorTable* G4PenelopeOscillatorManager::GetOscillatorTableCompton(const G4Material* mat)
{
  // (1) First time, create fOscillatorStore and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fOscillatorStoreCompton->count(mat))
    {
      //Ok, it exists
      return fOscillatorStoreCompton->find(mat)->second;
    }

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fOscillatorStoreCompton->count(mat))
    return fOscillatorStoreCompton->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetOscillatorTableCompton() " << G4endl;
      G4cout << "Impossible to create Compton oscillator table for " << mat->GetName() << G4endl;
      return nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeOscillator* G4PenelopeOscillatorManager::GetOscillatorCompton(const G4Material* material,
									G4int index)
{
  G4PenelopeOscillatorTable* theTable = GetOscillatorTableCompton(material);
  if (((std::size_t)index) < theTable->size())
    return (*theTable)[index];
  else
    {
      G4cout << "WARNING: Compton table for material " << material->GetName() << " has " <<
	theTable->size() << " oscillators" << G4endl;
      G4cout << "Oscillator #" << index << " cannot be retrieved" << G4endl;
      G4cout << "Returning null pointer" << G4endl;
      return nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeOscillatorManager::BuildOscillatorTable(const G4Material* material)
{
  //THIS CORRESPONDS TO THE ROUTINE PEMATW of PENELOPE

  G4double meanAtomExcitationEnergy[99] = {19.2*eV, 41.8*eV, 40.0*eV, 63.7*eV, 76.0*eV, 81.0*eV,
					   82.0*eV, 95.0*eV,115.0*eV,137.0*eV,149.0*eV,156.0*eV,
					   166.0*eV,
					   173.0*eV,173.0*eV,180.0*eV,174.0*eV,188.0*eV,190.0*eV,191.0*eV,
					   216.0*eV,233.0*eV,245.0*eV,257.0*eV,272.0*eV,286.0*eV,297.0*eV,
					   311.0*eV,322.0*eV,330.0*eV,334.0*eV,350.0*eV,347.0*eV,348.0*eV,
					   343.0*eV,352.0*eV,363.0*eV,366.0*eV,379.0*eV,393.0*eV,417.0*eV,
					   424.0*eV,428.0*eV,441.0*eV,449.0*eV,470.0*eV,470.0*eV,469.0*eV,
					   488.0*eV,488.0*eV,487.0*eV,485.0*eV,491.0*eV,482.0*eV,488.0*eV,
					   491.0*eV,501.0*eV,523.0*eV,535.0*eV,546.0*eV,560.0*eV,574.0*eV,
					   580.0*eV,591.0*eV,614.0*eV,628.0*eV,650.0*eV,658.0*eV,674.0*eV,
					   684.0*eV,694.0*eV,705.0*eV,718.0*eV,727.0*eV,736.0*eV,746.0*eV,
					   757.0*eV,790.0*eV,790.0*eV,800.0*eV,810.0*eV,823.0*eV,823.0*eV,
					   830.0*eV,825.0*eV,794.0*eV,827.0*eV,826.0*eV,841.0*eV,847.0*eV,
					   878.0*eV,890.0*eV,902.0*eV,921.0*eV,934.0*eV,939.0*eV,952.0*eV,
					   966.0*eV,980.0*eV};

  if (fVerbosityLevel > 0)
    G4cout << "Going to build Oscillator Table for " << material->GetName() << G4endl;

  G4int nElements = (G4int)material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  //At the moment, there's no way in Geant4 to know if a material
  //is defined with atom numbers or fraction of weigth
  const G4double* fractionVector = material->GetFractionVector();

  //Take always the composition by fraction of mass. For the composition by
  //atoms: it is calculated by Geant4 but with some rounding to integers
  G4double totalZ = 0;
  G4double totalMolecularWeight = 0;
  G4double meanExcitationEnergy = 0;

  std::vector<G4double> *StechiometricFactors = new std::vector<G4double>;

  for (G4int i=0;i<nElements;i++)
    {
      //G4int iZ = (G4int) (*elementVector)[i]->GetZ();
      G4double fraction = fractionVector[i];
      G4double atomicWeigth = (*elementVector)[i]->GetAtomicMassAmu();
      StechiometricFactors->push_back(fraction/atomicWeigth);
    }
  //Find max
  G4double MaxStechiometricFactor = 0.;
  for (G4int i=0;i<nElements;i++)
    {
      if ((*StechiometricFactors)[i] > MaxStechiometricFactor)
	MaxStechiometricFactor = (*StechiometricFactors)[i];
    }
  if (MaxStechiometricFactor<1e-16)
    {
      G4ExceptionDescription ed;
      ed << "Problem with the mass composition of " << material->GetName() << G4endl;
      ed << "MaxStechiometricFactor = " << MaxStechiometricFactor << G4endl;
      G4Exception("G4PenelopeOscillatorManager::BuildOscillatorTable()",
		  "em2035",FatalException,ed);
    }
  //Normalize
  for (G4int i=0;i<nElements;++i)
    (*StechiometricFactors)[i] /=  MaxStechiometricFactor;

  // Equivalent atoms per molecule
  G4double theatomsPerMolecule = 0;
  for (G4int i=0;i<nElements;i++)
    theatomsPerMolecule += (*StechiometricFactors)[i];
  G4double moleculeDensity =
    material->GetTotNbOfAtomsPerVolume()/theatomsPerMolecule; //molecules per unit volume

  if (fVerbosityLevel > 1)
    {
      for (std::size_t i=0;i<StechiometricFactors->size();++i)
	{
	  G4cout << "Element " << (*elementVector)[i]->GetSymbol() << " (Z = " <<
	    (*elementVector)[i]->GetZasInt() << ") --> " <<
	    (*StechiometricFactors)[i] << " atoms/molecule " << G4endl;
	}
    }

  for (G4int i=0;i<nElements;++i)
    {
      G4int iZ = (*elementVector)[i]->GetZasInt();
      totalZ += iZ * (*StechiometricFactors)[i];
      totalMolecularWeight += (*elementVector)[i]->GetAtomicMassAmu() * (*StechiometricFactors)[i];
      meanExcitationEnergy += iZ*G4Log(meanAtomExcitationEnergy[iZ-1])*(*StechiometricFactors)[i];
      std::pair<const G4Material*,G4int> theKey = std::make_pair(material,iZ);
      if (!fAtomTablePerMolecule->count(theKey))
	fAtomTablePerMolecule->insert(std::make_pair(theKey,(*StechiometricFactors)[i]));
    }
  meanExcitationEnergy = G4Exp(meanExcitationEnergy/totalZ);

  fAtomicNumber->insert(std::make_pair(material,totalZ));
  fAtomicMass->insert(std::make_pair(material,totalMolecularWeight));
  fExcitationEnergy->insert(std::make_pair(material,meanExcitationEnergy));
  fAtomsPerMolecule->insert(std::make_pair(material,theatomsPerMolecule));

  if (fVerbosityLevel > 1)
    {
      G4cout << "Calculated mean excitation energy for " << material->GetName() <<
	" = " << meanExcitationEnergy/eV << " eV" << G4endl;
    }

  std::vector<G4PenelopeOscillator> *helper = new std::vector<G4PenelopeOscillator>;

  //First Oscillator: conduction band. Tentativaly assumed to consist of valence electrons (each
  //atom contributes a number of electrons equal to its lowest chemical valence)
  G4PenelopeOscillator newOsc;
  newOsc.SetOscillatorStrength(0.);
  newOsc.SetIonisationEnergy(0*eV);
  newOsc.SetHartreeFactor(0);
  newOsc.SetParentZ(0);
  newOsc.SetShellFlag(30);
  newOsc.SetParentShellID(30); //does not correspond to any "real" level
  helper->push_back(newOsc);

  //Load elements and oscillators
  for (G4int k=0;k<nElements;k++)
    {
      G4double Z = (*elementVector)[k]->GetZ();
      G4bool finished = false;
      for (G4int i=0;i<2000 && !finished;i++)
	{
	  /*
	    fElementData[0][i] = Z;
	    fElementData[1][i] = shellCode;
	    fElementData[2][i] = occupationNumber;
	    fElementData[3][i] = ionisationEnergy;
	    fElementData[4][i] = hartreeProfile;
	  */
	  if (fElementData[0][i] == Z)
	    {
	      G4int shellID = (G4int) fElementData[1][i];
	      G4double occup = fElementData[2][i];
	      if (shellID > 0)
		{

		  if (std::fabs(occup) > 0)
		    {
		      G4PenelopeOscillator newOscLocal;
		      newOscLocal.SetOscillatorStrength(std::fabs(occup)*(*StechiometricFactors)[k]);
		      newOscLocal.SetIonisationEnergy(fElementData[3][i]);
		      newOscLocal.SetHartreeFactor(fElementData[4][i]/fine_structure_const);
		      newOscLocal.SetParentZ(fElementData[0][i]);
		      //keep track of the origianl shell level
		      newOscLocal.SetParentShellID((G4int)fElementData[1][i]);
		      //register only K, L and M shells. Outer shells all grouped with
		      //shellIndex = 30
		      if (fElementData[0][i] > 6 && fElementData[1][i] < 10)
			newOscLocal.SetShellFlag(((G4int)fElementData[1][i]));
		      else
			newOscLocal.SetShellFlag(30);
		      helper->push_back(newOscLocal);
		      if (occup < 0)
			{
			  G4double ff = (*helper)[0].GetOscillatorStrength();
			  ff += std::fabs(occup)*(*StechiometricFactors)[k];
			  (*helper)[0].SetOscillatorStrength(ff);
			}
		    }
		}
	    }
	  if (fElementData[0][i] > Z)
	    finished = true;
	}
    }

  delete StechiometricFactors;

  //NOW: sort oscillators according to increasing ionisation energy
  //Notice: it works because helper is a vector of _object_, not a
  //vector to _pointers_
  std::sort(helper->begin(),helper->end());

  // Plasma energy and conduction band excitation
  static const G4double RydbergEnergy = 13.60569*eV;
  G4double Omega = std::sqrt(4*pi*moleculeDensity*totalZ*Bohr_radius)*Bohr_radius*2.0*RydbergEnergy;
  G4double conductionStrength = (*helper)[0].GetOscillatorStrength();
  G4double plasmaEnergy = Omega*std::sqrt(conductionStrength/totalZ);

  fPlasmaSquared->insert(std::make_pair(material,Omega*Omega));

  G4bool isAConductor = false;
  G4int nullOsc = 0;

  if (fVerbosityLevel > 1)
    {
      G4cout << "Estimated oscillator strength and energy of plasmon: " <<
	conductionStrength << " and " << plasmaEnergy/eV << " eV" << G4endl;
    }

  if (conductionStrength < 0.01 || plasmaEnergy<1.0*eV) //this is an insulator
    {
      if (fVerbosityLevel >1 )
	G4cout << material->GetName() << " is an insulator " << G4endl;
      //remove conduction band oscillator
      helper->erase(helper->begin());
    }
  else //this is a conductor, Outer shells moved to conduction band
    {
      if (fVerbosityLevel >1 )
	G4cout << material->GetName() << " is a conductor " << G4endl;
      isAConductor = true;
      //copy the conduction strength.. The number is going to change.
      G4double conductionStrengthCopy = conductionStrength;
      G4bool quit = false;
      for (std::size_t i = 1; i<helper->size() && !quit ;++i)
	{
	  G4double oscStre = (*helper)[i].GetOscillatorStrength();
	  //loop is repeated over here
	  if (oscStre < conductionStrengthCopy)
	    {
	      conductionStrengthCopy = conductionStrengthCopy-oscStre;
	      (*helper)[i].SetOscillatorStrength(0.);
	      nullOsc++;
	    }
	  else //this is passed only once - no goto -
	    {
	      quit = true;
	      (*helper)[i].SetOscillatorStrength(oscStre-conductionStrengthCopy);
	      if (std::fabs((*helper)[i].GetOscillatorStrength()) < 1e-12)
		{
		  conductionStrength += (*helper)[i].GetOscillatorStrength();
		  (*helper)[i].SetOscillatorStrength(0.);
		  nullOsc++;
		}
	    }
	}
      //Update conduction band
      (*helper)[0].SetOscillatorStrength(conductionStrength);
      (*helper)[0].SetIonisationEnergy(0.);
      (*helper)[0].SetResonanceEnergy(plasmaEnergy);
      G4double hartree = 0.75/std::sqrt(3.0*pi*pi*moleculeDensity*
					Bohr_radius*Bohr_radius*Bohr_radius*conductionStrength);
      (*helper)[0].SetHartreeFactor(hartree/fine_structure_const);
  }

  //Check f-sum rule
  G4double sum = 0;
  for (std::size_t i=0;i<helper->size();++i)
    {
      sum += (*helper)[i].GetOscillatorStrength();
    }
  if (std::fabs(sum-totalZ) > (1e-6*totalZ))
    {
      G4ExceptionDescription ed;
      ed << "Inconsistent oscillator data for " << material->GetName() << G4endl;
      ed << sum << " " << totalZ << G4endl;
      G4Exception("G4PenelopeOscillatorManager::BuildOscillatorTable()",
		  "em2036",FatalException,ed);
    }
  if (std::fabs(sum-totalZ) > (1e-12*totalZ))
    {
      G4double fact = totalZ/sum;
      for (std::size_t i=0;i<helper->size();++i)
	{
	  G4double ff = (*helper)[i].GetOscillatorStrength()*fact;
	  (*helper)[i].SetOscillatorStrength(ff);
	}
    }

   //Remove null items
  for (G4int k=0;k<nullOsc;k++)
    {
      G4bool exit=false;
      for (std::size_t i=0;i<helper->size() && !exit;++i)
	{
	  if (std::fabs((*helper)[i].GetOscillatorStrength()) < 1e-12)
	    {
	      helper->erase(helper->begin()+i);
	      exit = true;
	    }
	}
    }

  //Sternheimer's adjustment factor
  G4double adjustmentFactor = 0;
  if (helper->size() > 1)
    {
      G4double TST = totalZ*G4Log(meanExcitationEnergy/eV);
      G4double AALow = 0.1;
      G4double AAHigh = 10.;
      do
	{
	  adjustmentFactor = (AALow+AAHigh)*0.5;
	  G4double sumLocal = 0;
	  for (std::size_t i=0;i<helper->size();++i)
	    {
	      if (i == 0 && isAConductor)
		{
		  G4double resEne = (*helper)[i].GetResonanceEnergy();
		  sumLocal += (*helper)[i].GetOscillatorStrength()*G4Log(resEne/eV);
		}
	      else
		{
		  G4double ionEne = (*helper)[i].GetIonisationEnergy();
		  G4double oscStre = (*helper)[i].GetOscillatorStrength();
		  G4double WI2 = (adjustmentFactor*adjustmentFactor*ionEne*ionEne) +
		    2./3.*(oscStre/totalZ)*Omega*Omega;
		  G4double resEne = std::sqrt(WI2);
		  (*helper)[i].SetResonanceEnergy(resEne);
		  sumLocal +=  (*helper)[i].GetOscillatorStrength()*G4Log(resEne/eV);
		}
	    }
	  if (sumLocal < TST)
	    AALow = adjustmentFactor;
	  else
	    AAHigh = adjustmentFactor;
	  if (fVerbosityLevel > 3)
	    G4cout << "Sternheimer's adjustment factor loops: " << AALow << " " << AAHigh << " " <<
	      adjustmentFactor << " " << TST << " " <<
	      sumLocal << G4endl;
	}while((AAHigh-AALow)>(1e-14*adjustmentFactor));
    }
  else
    {
      G4double ionEne = (*helper)[0].GetIonisationEnergy();
      (*helper)[0].SetIonisationEnergy(std::fabs(ionEne));
      (*helper)[0].SetResonanceEnergy(meanExcitationEnergy);
    }
  if (fVerbosityLevel > 1)
    {
      G4cout << "Sternheimer's adjustment factor: " << adjustmentFactor << G4endl;
    }

  //Check again for data consistency
  G4double xcheck = (*helper)[0].GetOscillatorStrength()*G4Log((*helper)[0].GetResonanceEnergy());
  G4double TST = (*helper)[0].GetOscillatorStrength();
  for (std::size_t i=1;i<helper->size();++i)
    {
      xcheck += (*helper)[i].GetOscillatorStrength()*G4Log((*helper)[i].GetResonanceEnergy());
      TST += (*helper)[i].GetOscillatorStrength();
    }
  if (std::fabs(TST-totalZ)>1e-8*totalZ)
    {
      G4ExceptionDescription ed;
      ed << "Inconsistent oscillator data " << G4endl;
      ed << TST << " " << totalZ << G4endl;
      G4Exception("G4PenelopeOscillatorManager::BuildOscillatorTable()",
		  "em2036",FatalException,ed);
    }
  xcheck = G4Exp(xcheck/totalZ);
  if (std::fabs(xcheck-meanExcitationEnergy) > 1e-8*meanExcitationEnergy)
    {
      G4ExceptionDescription ed;
      ed << "Error in Sterheimer factor calculation " << G4endl;
      ed << xcheck/eV << " " << meanExcitationEnergy/eV << G4endl;
      G4Exception("G4PenelopeOscillatorManager::BuildOscillatorTable()",
		  "em2037",FatalException,ed);
    }

  //Selection of the lowest ionisation energy for inner shells. Only the K, L and M shells with
  //ionisation energy less than the N1 shell of the heaviest element in the material are considered as
  //inner shells. As a results, the inner/outer shell character of an atomic shell depends on the
  //composition of the material.
  G4double Zmax = 0;
  for (G4int k=0;k<nElements;k++)
    {
      G4double Z = (*elementVector)[k]->GetZ();
      if (Z>Zmax) Zmax = Z;
    }
  //Find N1 level of the heaviest element (if any).
  G4bool found = false;
  G4double cutEnergy = 50*eV;
  for (std::size_t i=0;i<helper->size() && !found;++i)
    {
      G4double Z = (*helper)[i].GetParentZ();
      G4int shID = (*helper)[i].GetParentShellID(); //look for the N1 level
      if (shID == 10 && Z == Zmax)
	{
	  found = true;
	  if ((*helper)[i].GetIonisationEnergy() > cutEnergy)
	    cutEnergy = (*helper)[i].GetIonisationEnergy();
	}
    }
  //Make that cutEnergy cannot be higher than 250 eV, namely the fluorescence level by
  //Geant4
  G4double lowEnergyLimitForFluorescence = 250*eV;
  cutEnergy = std::min(cutEnergy,lowEnergyLimitForFluorescence);

  if (fVerbosityLevel > 1)
      G4cout << "Cutoff energy: " << cutEnergy/eV << " eV" << G4endl;
  //
  //Copy helper in the oscillatorTable for Ionisation
  //
  //Oscillator table Ionisation for the material
  G4PenelopeOscillatorTable* theTable = new G4PenelopeOscillatorTable(); //vector of oscillator
  G4PenelopeOscillatorResEnergyComparator comparator;
  std::sort(helper->begin(),helper->end(),comparator);

  //COPY THE HELPER (vector of object) to theTable (vector of Pointers).
  for (std::size_t i=0;i<helper->size();++i)
    {
      //copy content --> one may need it later (e.g. to fill another table, with variations)
      G4PenelopeOscillator* theOsc = new G4PenelopeOscillator((*helper)[i]);
      theTable->push_back(theOsc);
    }

  //Oscillators of outer shells with resonance energies differing by a factor less than
  //Rgroup are grouped as a single oscillator
  G4double Rgroup = 1.05;
  std::size_t Nost = theTable->size();

  std::size_t firstIndex = (isAConductor) ? 1 : 0; //for conductors, skip conduction oscillator
  G4bool loopAgain = false;
  G4int nLoops = 0;
  G4int removedLevels = 0;
  do
    {
      loopAgain = false;
      nLoops++;
      if (Nost>firstIndex+1)
	{	  
	  removedLevels = 0;
	  for (std::size_t i=firstIndex;i<theTable->size()-1;++i)
	    {
	      G4bool skipLoop = false;
	      G4int shellFlag = (*theTable)[i]->GetShellFlag();
	      G4double ionEne = (*theTable)[i]->GetIonisationEnergy();
	      G4double resEne = (*theTable)[i]->GetResonanceEnergy();
	      G4double resEnePlus1 = (*theTable)[i+1]->GetResonanceEnergy();
	      G4double oscStre = (*theTable)[i]->GetOscillatorStrength();
	      G4double oscStrePlus1 = (*theTable)[i+1]->GetOscillatorStrength();
	      //if (shellFlag < 10 && ionEne>cutEnergy) in Penelope
	      if (ionEne>cutEnergy) //remove condition that shellFlag < 10!
		skipLoop = true;
	      if (resEne<1.0*eV || resEnePlus1<1.0*eV)
		skipLoop = true;
	      if (resEnePlus1 > Rgroup*resEne)
		skipLoop = true;
	      if (!skipLoop)
		{
		  G4double newRes = G4Exp((oscStre*G4Log(resEne)+
					      oscStrePlus1*G4Log(resEnePlus1))
					     /(oscStre+oscStrePlus1));
		  (*theTable)[i]->SetResonanceEnergy(newRes);
		  G4double newIon = (oscStre*ionEne+
				     oscStrePlus1*(*theTable)[i+1]->GetIonisationEnergy())/
		    (oscStre+oscStrePlus1);
		  (*theTable)[i]->SetIonisationEnergy(newIon);
		  G4double newStre = oscStre+oscStrePlus1;
		  (*theTable)[i]->SetOscillatorStrength(newStre);
		  G4double newHartree = (oscStre*(*theTable)[i]->GetHartreeFactor()+
					 oscStrePlus1*(*theTable)[i+1]->GetHartreeFactor())/
		    (oscStre+oscStrePlus1);
		  (*theTable)[i]->SetHartreeFactor(newHartree);
		  if ((*theTable)[i]->GetParentZ() != (*theTable)[i+1]->GetParentZ())
		    (*theTable)[i]->SetParentZ(0.);
		  if (shellFlag < 10 || (*theTable)[i+1]->GetShellFlag() < 10)
		    {
		      G4int newFlag = std::min(shellFlag,(*theTable)[i+1]->GetShellFlag());
		      (*theTable)[i]->SetShellFlag(newFlag);
		    }
		  else
		    (*theTable)[i]->SetShellFlag(30);
		  //We've lost anyway the track of the original level
		  (*theTable)[i]->SetParentShellID((*theTable)[i]->GetShellFlag());


		  if (i<theTable->size()-2)
		    {
		      for (std::size_t ii=i+1;ii<theTable->size()-1;++ii)
			(*theTable)[ii] = (*theTable)[ii+1];
		    }
		  //G4cout << theTable->size() << G4endl;
		  theTable->erase(theTable->begin()+theTable->size()-1); //delete last element
		  removedLevels++;
		}
	    }
	}
      if (removedLevels)
	{
	  Nost -= removedLevels;
	  loopAgain = true;
	}
      if (Rgroup < 1.414213 || Nost > 64)
	{
	  Rgroup = Rgroup*Rgroup;
	  loopAgain = true;
	}
      //Add protection against infinite loops here
      if (nLoops > 100 && !removedLevels)
	loopAgain = false;
    }while(loopAgain);

  if (fVerbosityLevel > 1)
    {
      G4cout << "Final grouping factor for Ionisation: " << Rgroup << G4endl;
    }

  //Final Electron/Positron model parameters
  for (std::size_t i=0;i<theTable->size();++i)
    {
      //Set cutoff recoil energy for the resonant mode
      G4double ionEne = (*theTable)[i]->GetIonisationEnergy();
      if (ionEne < 1e-3*eV)
	{
	  G4double resEne = (*theTable)[i]->GetResonanceEnergy();
	  (*theTable)[i]->SetIonisationEnergy(0.*eV);
	  (*theTable)[i]->SetCutoffRecoilResonantEnergy(resEne);
	}
      else
	(*theTable)[i]->SetCutoffRecoilResonantEnergy(ionEne);
    }

  //Last step
  fOscillatorStoreIonisation->insert(std::make_pair(material,theTable));

  /******************************************
    SAME FOR COMPTON
  ******************************************/
  //
  //Copy helper in the oscillatorTable for Compton
  //
  //Oscillator table Ionisation for the material
  G4PenelopeOscillatorTable* theTableC = new G4PenelopeOscillatorTable(); //vector of oscillator
  //order by ionisation energy
  std::sort(helper->begin(),helper->end());
  //COPY THE HELPER (vector of object) to theTable (vector of Pointers).
  for (std::size_t i=0;i<helper->size();++i)
    {
      //copy content --> one may need it later (e.g. to fill another table, with variations)
      G4PenelopeOscillator* theOsc = new G4PenelopeOscillator((*helper)[i]);
      theTableC->push_back(theOsc);
    }
  //Oscillators of outer shells with resonance energies differing by a factor less than
  //Rgroup are grouped as a single oscillator
  Rgroup = 1.5;
  Nost = theTableC->size();

  firstIndex = (isAConductor) ? 1 : 0; //for conductors, skip conduction oscillator
  loopAgain = false;
  removedLevels = 0;
  do
    {
      nLoops++;
      loopAgain = false;
      if (Nost>firstIndex+1)
	{
	  removedLevels = 0;
	  for (std::size_t i=firstIndex;i<theTableC->size()-1;++i)
	    {
	      G4bool skipLoop = false;
	      G4double ionEne = (*theTableC)[i]->GetIonisationEnergy();
	      G4double ionEnePlus1 = (*theTableC)[i+1]->GetIonisationEnergy();
	      G4double oscStre = (*theTableC)[i]->GetOscillatorStrength();
	      G4double oscStrePlus1 = (*theTableC)[i+1]->GetOscillatorStrength();
	      //if (shellFlag < 10 && ionEne>cutEnergy) in Penelope
	      if (ionEne>cutEnergy)
		skipLoop = true;
	      if (ionEne<1.0*eV || ionEnePlus1<1.0*eV)
		skipLoop = true;
	      if (ionEnePlus1 > Rgroup*ionEne)
		skipLoop = true;

	      if (!skipLoop)
		{
		  G4double newIon = (oscStre*ionEne+
				     oscStrePlus1*ionEnePlus1)/
		    (oscStre+oscStrePlus1);
		  (*theTableC)[i]->SetIonisationEnergy(newIon);
		  G4double newStre = oscStre+oscStrePlus1;
		  (*theTableC)[i]->SetOscillatorStrength(newStre);
		  G4double newHartree = (oscStre*(*theTableC)[i]->GetHartreeFactor()+
					 oscStrePlus1*(*theTableC)[i+1]->GetHartreeFactor())/
		    (oscStre+oscStrePlus1);
		  (*theTableC)[i]->SetHartreeFactor(newHartree);
		  if ((*theTableC)[i]->GetParentZ() != (*theTableC)[i+1]->GetParentZ())
		    (*theTableC)[i]->SetParentZ(0.);
		  (*theTableC)[i]->SetShellFlag(30);
		  (*theTableC)[i]->SetParentShellID((*theTableC)[i]->GetShellFlag());

		  if (i<theTableC->size()-2)
		    {
		      for (std::size_t ii=i+1;ii<theTableC->size()-1;++ii)
			(*theTableC)[ii] = (*theTableC)[ii+1];
		    }
		  theTableC->erase(theTableC->begin()+theTableC->size()-1); //delete last element
		  removedLevels++;
		}
	    }
	}
      if (removedLevels)
	{
	  Nost -= removedLevels;
	  loopAgain = true;
	}
      if (Rgroup < 2.0 || Nost > 64)
	{
	  Rgroup = Rgroup*Rgroup;
	  loopAgain = true;
	}
      //Add protection against infinite loops here
      if (nLoops > 100 && !removedLevels)
	loopAgain = false;
    }while(loopAgain);


   if (fVerbosityLevel > 1)
    {
      G4cout << "Final grouping factor for Compton: " << Rgroup << G4endl;
    }

   //Last step
   fOscillatorStoreCompton->insert(std::make_pair(material,theTableC));
   
   //CLEAN UP theHelper and its content
   delete helper;
   if (fVerbosityLevel > 1)
     Dump(material);
   
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeOscillatorManager::ReadElementData()
{
  if (fVerbosityLevel > 0)
    {
      G4cout << "G4PenelopeOscillatorManager::ReadElementData()" << G4endl;
      G4cout << "Going to read Element Data" << G4endl;
    }
    const char* path = G4FindDataDir("G4LEDATA");
    if(!path)
    {
      G4String excep = "G4PenelopeOscillatorManager - G4LEDATA environment variable not set!";
      G4Exception("G4PenelopeOscillatorManager::ReadElementData()",
		  "em0006",FatalException,excep);
      return;
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/pdatconf.p08";
  std::ifstream file(pathFile);

  if (!file.is_open())
    {
      G4String excep = "G4PenelopeOscillatorManager - data file " + pathFile + " not found!";
      G4Exception("G4PenelopeOscillatorManager::ReadElementData()",
		  "em0003",FatalException,excep);
    }

  G4AtomicTransitionManager* theTransitionManager =
    G4AtomicTransitionManager::Instance();
  theTransitionManager->Initialise();

  //Read header (22 lines)
  G4String theHeader;
  for (G4int iline=0;iline<22;iline++)
    getline(file,theHeader);
  //Done
  G4int Z=0;
  G4int shellCode = 0;
  G4String shellId = "NULL";
  G4int occupationNumber = 0;
  G4double ionisationEnergy = 0.0*eV;
  G4double hartreeProfile = 0.;
  G4int shellCounter = 0;
  G4int oldZ = -1;
  G4int numberOfShells = 0;
  //Start reading data
  for (G4int i=0;!file.eof();i++)
    {
      file >> Z >> shellCode >> shellId >> occupationNumber >> ionisationEnergy >> hartreeProfile;
      if (Z>0 && i<2000)
	{
	  fElementData[0][i] = Z;
	  fElementData[1][i] = shellCode;
	  fElementData[2][i] = occupationNumber;
	  //reset things
	  if (Z != oldZ)
	    {
	      shellCounter = 0;
	      oldZ = Z;
	      numberOfShells = theTransitionManager->NumberOfShells(Z);
	    }
	  G4double bindingEnergy = -1*eV;
	  if (shellCounter<numberOfShells)
	    {
	      G4AtomicShell* shell = theTransitionManager->Shell(Z,shellCounter);
	      bindingEnergy = shell->BindingEnergy();
	    }
	  //Valid level found in the G4AtomicTransition database: keep it, otherwise use
	  //the ionisation energy found in the Penelope database
	  fElementData[3][i] = (bindingEnergy>100*eV) ? bindingEnergy : ionisationEnergy*eV;
	  fElementData[4][i] = hartreeProfile;
	  shellCounter++;
	}
    }
  file.close();

  if (fVerbosityLevel > 1)
    {
      G4cout << "G4PenelopeOscillatorManager::ReadElementData(): Data file read" << G4endl;
    }
  fReadElementData = true;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PenelopeOscillatorManager::GetMeanExcitationEnergy(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fExcitationEnergy->count(mat))
    return fExcitationEnergy->find(mat)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fExcitationEnergy->count(mat))
    return fExcitationEnergy->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetMolecularExcitationEnergy() " << G4endl;
      G4cout << "Impossible to retrieve the excitation energy for  " << mat->GetName() << G4endl;
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4PenelopeOscillatorManager::GetPlasmaEnergySquared(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fPlasmaSquared->count(mat))
    return fPlasmaSquared->find(mat)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fPlasmaSquared->count(mat))
    return fPlasmaSquared->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetPlasmaEnergySquared() " << G4endl;
      G4cout << "Impossible to retrieve the plasma energy for  " << mat->GetName() << G4endl;
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeOscillatorManager::GetAtomsPerMolecule(const G4Material* mat)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material has been already included
  if (fAtomsPerMolecule->count(mat))
    return fAtomsPerMolecule->find(mat)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fAtomsPerMolecule->count(mat))
    return fAtomsPerMolecule->find(mat)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetAtomsPerMolecule() " << G4endl;
      G4cout << "Impossible to retrieve the number of atoms per molecule for  "
	     << mat->GetName() << G4endl;
      return 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeOscillatorManager::GetNumberOfZAtomsPerMolecule(const G4Material* mat,G4int Z)
{
  // (1) First time, create fOscillatorStores and read data
  CheckForTablesCreated();

  // (2) Check if the material/Z couple has been already included
  std::pair<const G4Material*,G4int> theKey = std::make_pair(mat,Z);
  if (fAtomTablePerMolecule->count(theKey))
    return fAtomTablePerMolecule->find(theKey)->second;

  // (3) If we are here, it means that we have to create the table for the material
  BuildOscillatorTable(mat);

  // (4) now, the oscillator store should be ok
  if (fAtomTablePerMolecule->count(theKey))
    return fAtomTablePerMolecule->find(theKey)->second;
  else
    {
      G4cout << "G4PenelopeOscillatorManager::GetAtomsPerMolecule() " << G4endl;
      G4cout << "Impossible to retrieve the number of atoms per molecule for Z = "
	     << Z << " in material " << mat->GetName() << G4endl;
      return 0;
    }
}
