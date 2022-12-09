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
// 03 Dec 2009   L Pandola    First implementation
// 25 May 2011   L.Pandola    Renamed (make v2008 as default Penelope)
// 19 Sep 2013   L.Pandola    Migration to MT
//

#include "G4PenelopeRayleighModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PenelopeSamplingData.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4AutoLock.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4int G4PenelopeRayleighModel::fMaxZ;
G4PhysicsFreeVector* G4PenelopeRayleighModel::fLogAtomicCrossSection[] = {nullptr};
G4PhysicsFreeVector* G4PenelopeRayleighModel::fAtomicFormFactor[] = {nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeRayleighModel::G4PenelopeRayleighModel(const G4ParticleDefinition* part,
						 const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),fParticle(nullptr),
   fLogFormFactorTable(nullptr),fPMaxTable(nullptr),fSamplingTable(nullptr),
   fIsInitialised(false),fLocalTable(false)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);

  if (part)
    SetParticle(part);

  //
  fVerboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //build the energy grid. It is the same for all materials
  G4double logenergy = G4Log(fIntrinsicLowEnergyLimit/2.);
  G4double logmaxenergy = G4Log(1.5*fIntrinsicHighEnergyLimit);
  //finer grid below 160 keV
  G4double logtransitionenergy = G4Log(160*keV);
  G4double logfactor1 = G4Log(10.)/250.;
  G4double logfactor2 = logfactor1*10;
  fLogEnergyGridPMax.push_back(logenergy);
  do{
    if (logenergy < logtransitionenergy)
      logenergy += logfactor1;
    else
      logenergy += logfactor2;
    fLogEnergyGridPMax.push_back(logenergy);
  }while (logenergy < logmaxenergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeRayleighModel::~G4PenelopeRayleighModel()
{
  if (IsMaster() || fLocalTable)
    {
      
      for(G4int i=0; i<=fMaxZ; ++i) 
	{
	  if(fLogAtomicCrossSection[i]) 
	    { 
	      delete fLogAtomicCrossSection[i];
	      fLogAtomicCrossSection[i] = nullptr;
	    }
	  if(fAtomicFormFactor[i])
	    {
	      delete fAtomicFormFactor[i];
	      fAtomicFormFactor[i] = nullptr;
	    }
	}
      ClearTables();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PenelopeRayleighModel::ClearTables()
{
   if (fLogFormFactorTable)
     {
       for (auto& item : (*fLogFormFactorTable))
	 if (item.second) delete item.second;
       delete fLogFormFactorTable;
       fLogFormFactorTable = nullptr; //zero explicitly
     }
   if (fPMaxTable)
     {
       for (auto& item : (*fPMaxTable))
	 if (item.second) delete item.second;
       delete fPMaxTable;
       fPMaxTable = nullptr; //zero explicitly
     }
   if (fSamplingTable)
     {
       for (auto& item : (*fSamplingTable))
	 if (item.second) delete item.second;
       delete fSamplingTable;
       fSamplingTable = nullptr; //zero explicitly
     }
   return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::Initialise(const G4ParticleDefinition* part,
					 const G4DataVector& )
{
  if (fVerboseLevel > 3)
    G4cout << "Calling G4PenelopeRayleighModel::Initialise()" << G4endl;

  SetParticle(part);

  //Only the master model creates/fills/destroys the tables
  if (IsMaster() && part == fParticle)
    {
      //clear tables depending on materials, not the atomic ones
      ClearTables();

      if (fVerboseLevel > 3)
	G4cout << "Calling G4PenelopeRayleighModel::Initialise() [master]" << G4endl;

      //create new tables  
      if (!fLogFormFactorTable)
	fLogFormFactorTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
      if (!fPMaxTable)
	fPMaxTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
      if (!fSamplingTable)
	fSamplingTable = new std::map<const G4Material*,G4PenelopeSamplingData*>;

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
	      if (!fLogAtomicCrossSection[iZ])
		ReadDataFile(iZ);
	    }

	  //1) If the table has not been built for the material, do it!
	  if (!fLogFormFactorTable->count(material))
	    BuildFormFactorTable(material);

	  //2) retrieve or build the sampling table
	  if (!(fSamplingTable->count(material)))
	    InitializeSamplingAlgorithm(material);

	  //3) retrieve or build the pMax data
	  if (!fPMaxTable->count(material))
	    GetPMaxTable(material);
	}

      if (fVerboseLevel > 1) {
	G4cout << "Penelope Rayleigh model v2008 is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / keV << " keV - "
	       << HighEnergyLimit() / GeV << " GeV"
	       << G4endl;
      }
    }

  if(fIsInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  fIsInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::InitialiseLocal(const G4ParticleDefinition* part,
						     G4VEmModel *masterModel)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopeRayleighModel::InitialiseLocal()" << G4endl;
  //
  //Check that particle matches: one might have multiple master models (e.g.
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeRayleighModel* theModel =
	static_cast<G4PenelopeRayleighModel*> (masterModel);

      //Copy pointers to the data tables
      for(G4int i=0; i<=fMaxZ; ++i) 
	{
	  fLogAtomicCrossSection[i] = theModel->fLogAtomicCrossSection[i];
	  fAtomicFormFactor[i] = theModel->fAtomicFormFactor[i];
	}
      fLogFormFactorTable = theModel->fLogFormFactorTable;
      fPMaxTable = theModel->fPMaxTable;
      fSamplingTable = theModel->fSamplingTable;

      //copy the G4DataVector with the grid
      fLogQSquareGrid = theModel->fLogQSquareGrid;

      //Same verbosity for all workers, as the master
      fVerboseLevel = theModel->fVerboseLevel;
    }

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
namespace { G4Mutex  PenelopeRayleighModelMutex = G4MUTEX_INITIALIZER; }
G4double G4PenelopeRayleighModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
							     G4double energy,
							     G4double Z,
							     G4double,
							     G4double,
							     G4double)
{
  // Cross section of Rayleigh scattering in Penelope v2008 is calculated by the EPDL97
  // tabulation, Cuellen et al. (1997), with non-relativistic form factors from Hubbel
  // et al. J. Phys. Chem. Ref. Data 4 (1975) 471; Erratum ibid. 6 (1977) 615.

   if (fVerboseLevel > 3)
    G4cout << "Calling CrossSectionPerAtom() of G4PenelopeRayleighModel" << G4endl;

   G4int iZ = G4int(Z);

   if (!fLogAtomicCrossSection[iZ])
     {
       //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
       //not filled up. This can happen in a UnitTest or via G4EmCalculator
       if (fVerboseLevel > 0)
	{
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to retrieve the cross section table for Z=" << iZ << G4endl;
	  ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	  G4Exception("G4PenelopeRayleighModel::ComputeCrossSectionPerAtom()",
		      "em2040",JustWarning,ed);
	}
       //protect file reading via autolock
       G4AutoLock lock(&PenelopeRayleighModelMutex);
       ReadDataFile(iZ);
       lock.unlock();
     }

   G4double cross = 0;
   G4PhysicsFreeVector* atom = fLogAtomicCrossSection[iZ];
   if (!atom)
     {
       G4ExceptionDescription ed;
       ed << "Unable to find Z=" << iZ << " in the atomic cross section table" << G4endl;
       G4Exception("G4PenelopeRayleighModel::ComputeCrossSectionPerAtom()",
		   "em2041",FatalException,ed);
       return 0;
     }
   G4double logene = G4Log(energy);
   G4double logXS = atom->Value(logene);
   cross = G4Exp(logXS);

   if (fVerboseLevel > 2)
     {
       G4cout << "Rayleigh cross section at " << energy/keV << " keV for Z=" << Z <<
	 " = " << cross/barn << " barn" << G4endl;
     }
   return cross;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4PenelopeRayleighModel::BuildFormFactorTable(const G4Material* material)
{
  /*
    1) get composition and equivalent molecular density
  */
  std::size_t nElements = material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();
  const G4double* fractionVector = material->GetFractionVector();

  std::vector<G4double> *StechiometricFactors = new std::vector<G4double>;
  for (std::size_t i=0;i<nElements;++i)
    {
      G4double fraction = fractionVector[i];
      G4double atomicWeigth = (*elementVector)[i]->GetA()/(g/mole);
      StechiometricFactors->push_back(fraction/atomicWeigth);
    }
  //Find max
  G4double MaxStechiometricFactor = 0.;
  for (std::size_t i=0;i<nElements;++i)
    {
      if ((*StechiometricFactors)[i] > MaxStechiometricFactor)
        MaxStechiometricFactor = (*StechiometricFactors)[i];
    }
  if (MaxStechiometricFactor<1e-16)
    {
      G4ExceptionDescription ed;
      ed << "Inconsistent data of atomic composition for " <<
	material->GetName() << G4endl;
      G4Exception("G4PenelopeRayleighModel::BuildFormFactorTable()",
		  "em2042",FatalException,ed);
    }
  //Normalize
  for (std::size_t i=0;i<nElements;++i)
    (*StechiometricFactors)[i] /=  MaxStechiometricFactor;

  /*
    CREATE THE FORM FACTOR TABLE
  */
  G4PhysicsFreeVector* theFFVec = new G4PhysicsFreeVector(fLogQSquareGrid.size(),/*spline=*/true);

  for (std::size_t k=0;k<fLogQSquareGrid.size();++k)
    {
      G4double ff2 = 0; //squared form factor
      for (std::size_t i=0;i<nElements;++i)
	{
	  G4int iZ = (*elementVector)[i]->GetZasInt();
	  G4PhysicsFreeVector* theAtomVec = fAtomicFormFactor[iZ];
	  G4double f = (*theAtomVec)[k]; //the q-grid is always the same
	  ff2 += f*f*(*StechiometricFactors)[i];
	}
      if (ff2)
	theFFVec->PutValue(k,fLogQSquareGrid[k],G4Log(ff2)); //NOTICE: THIS IS log(Q^2) vs. log(F^2)
    }
  theFFVec->FillSecondDerivatives(); //vector is filled!
  fLogFormFactorTable->insert(std::make_pair(material,theFFVec));

  delete StechiometricFactors;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* ,
						const G4MaterialCutsCouple* couple,
						const G4DynamicParticle* aDynamicGamma,
						G4double,
						G4double)
{
  // Sampling of the Rayleigh final state (namely, scattering angle of the photon)
  // from the Penelope2008 model. The scattering angle is sampled from the atomic
  // cross section dOmega/d(cosTheta) from Born ("Atomic Phyisics", 1969), disregarding
  // anomalous scattering effects. The Form Factor F(Q) function which appears in the
  // analytical cross section is retrieved via the method GetFSquared(); atomic data
  // are tabulated for F(Q). Form factor for compounds is calculated according to
  // the additivity rule. The sampling from the F(Q) is made via a Rational Inverse
  // Transform with Aliasing (RITA) algorithm; RITA parameters are calculated once
  // for each material and managed by G4PenelopeSamplingData objects.
  // The sampling algorithm (rejection method) has efficiency 67% at low energy, and
  // increases with energy. For E=100 keV the efficiency is 100% and 86% for
  // hydrogen and uranium, respectively.

  if (fVerboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeRayleighModel" << G4endl;

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  if (photonEnergy0 <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
      return ;
    }

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  const G4Material* theMat = couple->GetMaterial();

  //1) Verify if tables are ready
  //Either Initialize() was not called, or we are in a slave and InitializeLocal() was
  //not invoked
  if (!fPMaxTable || !fSamplingTable || !fLogFormFactorTable)
    {
      //create a **thread-local** version of the table. Used only for G4EmCalculator and
      //Unit Tests
      fLocalTable = true;
      if (!fLogFormFactorTable)
	fLogFormFactorTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
      if (!fPMaxTable)
	fPMaxTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
      if (!fSamplingTable)
	fSamplingTable = new std::map<const G4Material*,G4PenelopeSamplingData*>;
    }

  if (!fSamplingTable->count(theMat))
    {
      //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
      //not filled up. This can happen in a UnitTest
      if (fVerboseLevel > 0)
	{
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to find the fSamplingTable data for " <<
	    theMat->GetName() << G4endl;
	  ed << "This can happen only in Unit Tests" << G4endl;
	  G4Exception("G4PenelopeRayleighModel::SampleSecondaries()",
		      "em2019",JustWarning,ed);
	}
      const G4ElementVector* theElementVector = theMat->GetElementVector();
      //protect file reading via autolock
      G4AutoLock lock(&PenelopeRayleighModelMutex);
      for (std::size_t j=0;j<theMat->GetNumberOfElements();++j)
	{
	  G4int iZ = theElementVector->at(j)->GetZasInt();
	  if (!fLogAtomicCrossSection[iZ])
	    {
	      lock.lock();
	      ReadDataFile(iZ);
	      lock.unlock();
	    }
	}
      lock.lock();
      //1) If the table has not been built for the material, do it!
      if (!fLogFormFactorTable->count(theMat))
	BuildFormFactorTable(theMat);

      //2) retrieve or build the sampling table
      if (!(fSamplingTable->count(theMat)))
	InitializeSamplingAlgorithm(theMat);

      //3) retrieve or build the pMax data
      if (!fPMaxTable->count(theMat))
	GetPMaxTable(theMat);
      lock.unlock();
    }

  //Ok, restart the job
  G4PenelopeSamplingData* theDataTable = fSamplingTable->find(theMat)->second;
  G4PhysicsFreeVector* thePMax = fPMaxTable->find(theMat)->second;

  G4double cosTheta = 1.0;

  //OK, ready to go!
  G4double qmax = 2.0*photonEnergy0/electron_mass_c2; //this is non-dimensional now

  if (qmax < 1e-10) //very low momentum transfer
    {
      G4bool loopAgain=false;
      do
	{
	  loopAgain = false;
	  cosTheta = 1.0-2.0*G4UniformRand();
	  G4double G = 0.5*(1+cosTheta*cosTheta);
	  if (G4UniformRand()>G)
	    loopAgain = true;
	}while(loopAgain);
    }
  else //larger momentum transfer
    {
      std::size_t nData = theDataTable->GetNumberOfStoredPoints();
      G4double LastQ2inTheTable = theDataTable->GetX(nData-1);
      G4double q2max = std::min(qmax*qmax,LastQ2inTheTable);

      G4bool loopAgain = false;
      G4double MaxPValue = thePMax->Value(photonEnergy0);
      G4double xx=0;

      //Sampling by rejection method. The rejection function is
      //G = 0.5*(1+cos^2(theta))
      //
      do{
	loopAgain = false;
	G4double RandomMax = G4UniformRand()*MaxPValue;
	xx = theDataTable->SampleValue(RandomMax);
	//xx is a random value of q^2 in (0,q2max),sampled according to
	//F(Q^2) via the RITA algorithm
	if (xx > q2max)
	  loopAgain = true;
	cosTheta = 1.0-2.0*xx/q2max;
	G4double G = 0.5*(1+cosTheta*cosTheta);
	if (G4UniformRand()>G)
	  loopAgain = true;
      }while(loopAgain);
    }

  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);

  // Scattered photon angles. ( Z - axis along the parent photon)
  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*std::cos(phi);
  G4double dirY = sinTheta*std::sin(phi);
  G4double dirZ = cosTheta;

  // Update G4VParticleChange for the scattered photon
  G4ThreeVector photonDirection1(dirX, dirY, dirZ);

  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;
  fParticleChange->SetProposedKineticEnergy(photonEnergy0) ;

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::ReadDataFile(const G4int Z)
{
  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopeRayleighModel::ReadDataFile()" << G4endl;
      G4cout << "Going to read Rayleigh data files for Z=" << Z << G4endl;
    }
    const char* path = G4FindDataDir("G4LEDATA");
    if(!path)
    {
      G4String excep = "G4LEDATA environment variable not set!";
      G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		  "em0006",FatalException,excep);
      return;
    }

  /*
    Read first the cross section file
  */
  std::ostringstream ost;
  if (Z>9)
    ost << path << "/penelope/rayleigh/pdgra" << Z << ".p08";
  else
    ost << path << "/penelope/rayleigh/pdgra0" << Z << ".p08";
  std::ifstream file(ost.str().c_str());
  if (!file.is_open())
    {
      G4String excep = "Data file " + G4String(ost.str()) + " not found!";
      G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		  "em0003",FatalException,excep);
    }
  G4int readZ =0;
  std::size_t nPoints= 0;
  file >> readZ >> nPoints;
  //check the right file is opened.
  if (readZ != Z || nPoints <= 0 || nPoints >= 5000)
    {
      G4ExceptionDescription ed;
      ed << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		  "em0005",FatalException,ed);
      return;
    }

  fLogAtomicCrossSection[Z] = new G4PhysicsFreeVector((std::size_t)nPoints);
  G4double ene=0,f1=0,f2=0,xs=0;
  for (std::size_t i=0;i<nPoints;++i)
    {
      file >> ene >> f1 >> f2 >> xs;
      //dimensional quantities
      ene *= eV;
      xs *= cm2;
      fLogAtomicCrossSection[Z]->PutValue(i,G4Log(ene),G4Log(xs));
      if (file.eof() && i != (nPoints-1)) //file ended too early
	{
	  G4ExceptionDescription ed ;
	  ed << "Corrupted data file for Z=" << Z << G4endl;
	  ed << "Found less than " << nPoints << "entries " <<G4endl;
	  G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		      "em0005",FatalException,ed);
	}
    }
  file.close();

  /*
    Then read the form factor file
  */
  std::ostringstream ost2;
  if (Z>9)
    ost2 << path << "/penelope/rayleigh/pdaff" << Z << ".p08";
  else
    ost2 << path << "/penelope/rayleigh/pdaff0" << Z << ".p08";
  file.open(ost2.str().c_str());
  if (!file.is_open())
    {
      G4String excep = "Data file " + G4String(ost2.str()) + " not found!";
      G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		  "em0003",FatalException,excep);
    }
  file >> readZ >> nPoints;
  //check the right file is opened.
  if (readZ != Z || nPoints <= 0 || nPoints >= 5000)
    {
      G4ExceptionDescription ed;
      ed << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		  "em0005",FatalException,ed);
      return;
    }
  fAtomicFormFactor[Z] = new G4PhysicsFreeVector((std::size_t)nPoints);
  G4double q=0,ff=0,incoh=0;
  G4bool fillQGrid = false;
  //fill this vector only the first time.
  if (!fLogQSquareGrid.size())
    fillQGrid = true;
  for (std::size_t i=0;i<nPoints;++i)
    {
      file >> q >> ff >> incoh;
      //q and ff are dimensionless (q is in units of (m_e*c)
      fAtomicFormFactor[Z]->PutValue(i,q,ff);
      if (fillQGrid)
	{
	  fLogQSquareGrid.push_back(2.0*G4Log(q));
	}
      if (file.eof() && i != (nPoints-1)) //file ended too early
	{
	  G4ExceptionDescription ed;
	  ed << "Corrupted data file for Z=" << Z << G4endl;
	  ed << "Found less than " << nPoints << "entries " <<G4endl;
	  G4Exception("G4PenelopeRayleighModel::ReadDataFile()",
		      "em0005",FatalException,ed);
	}
    }
  file.close();
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeRayleighModel::GetFSquared(const G4Material* mat, const G4double QSquared)
{
  G4double f2 = 0;
  //Input value QSquared could be zero: protect the log() below against
  //the FPE exception
  //If Q<1e-10, set Q to 1e-10
  G4double logQSquared = (QSquared>1e-10) ? G4Log(QSquared) : -23.;
  //last value of the table
  G4double maxlogQ2 = fLogQSquareGrid[fLogQSquareGrid.size()-1];

  //now it should  be all right
  G4PhysicsFreeVector* theVec = fLogFormFactorTable->find(mat)->second;

  if (!theVec)
    {
      G4ExceptionDescription ed;
      ed << "Unable to retrieve F squared table for " << mat->GetName() << G4endl;
      G4Exception("G4PenelopeRayleighModel::GetFSquared()",
		  "em2046",FatalException,ed);
      return 0;
    }
  if (logQSquared < -20) // Q < 1e-9
    {
      G4double logf2 = (*theVec)[0]; //first value of the table
      f2 = G4Exp(logf2);
    }
  else if (logQSquared > maxlogQ2)
    f2 =0;
  else
    {
      //log(Q^2) vs. log(F^2)
      G4double logf2 = theVec->Value(logQSquared);
      f2 = G4Exp(logf2);

    }
  if (fVerboseLevel > 3)
    {
      G4cout << "G4PenelopeRayleighModel::GetFSquared() in " << mat->GetName() << G4endl;
      G4cout << "Q^2 = " <<  QSquared << " (units of 1/(m_e*c); F^2 = " << f2 << G4endl;
    }
  return f2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::InitializeSamplingAlgorithm(const G4Material* mat)
{
  G4double q2min = 0;
  G4double q2max = 0;
  const std::size_t np = 150; //hard-coded in Penelope
  //G4cout << "Init N= " << fLogQSquareGrid.size() << G4endl;
  for (std::size_t i=1;i<fLogQSquareGrid.size();++i)
    {
      G4double Q2 = G4Exp(fLogQSquareGrid[i]);
      if (GetFSquared(mat,Q2) >  1e-35)
	{
	  q2max = G4Exp(fLogQSquareGrid[i-1]);
	}
      //G4cout << "Q2= " << Q2 << " q2max= " << q2max << G4endl;
    }

  std::size_t nReducedPoints = np/4;

  //check for errors
  if (np < 16)
    {
      G4Exception("G4PenelopeRayleighModel::InitializeSamplingAlgorithm()",
		  "em2047",FatalException,
		  "Too few points to initialize the sampling algorithm");
    }
  if (q2min > (q2max-1e-10))
    {
      G4cout << "q2min= " << q2min << " q2max= " << q2max << G4endl;
      G4Exception("G4PenelopeRayleighModel::InitializeSamplingAlgorithm()",
		  "em2048",FatalException,
		  "Too narrow grid to initialize the sampling algorithm");
    }

  //This is subroutine RITAI0 of Penelope
  //Create an object of type G4PenelopeRayleighSamplingData --> store in a map::Material*

  //temporary vectors --> Then everything is stored in G4PenelopeSamplingData
  G4DataVector* x = new G4DataVector();

  /*******************************************************************************
    Start with a grid of NUNIF points uniformly spaced in the interval q2min,q2max
  ********************************************************************************/
  std::size_t NUNIF = std::min(std::max(((std::size_t)8),nReducedPoints),np/2);
  const G4int nip = 51; //hard-coded in Penelope

  G4double dx = (q2max-q2min)/((G4double) NUNIF-1);
  x->push_back(q2min);
  for (std::size_t i=1;i<NUNIF-1;++i)
    {
      G4double app = q2min + i*dx;
      x->push_back(app); //increase
    }
  x->push_back(q2max);

  if (fVerboseLevel> 3)
    G4cout << "Vector x has " << x->size() << " points, while NUNIF = " << NUNIF << G4endl;

  //Allocate temporary storage vectors
  G4DataVector* area = new G4DataVector();
  G4DataVector* a = new G4DataVector();
  G4DataVector* b = new G4DataVector();
  G4DataVector* c = new G4DataVector();
  G4DataVector* err = new G4DataVector();

  for (std::size_t i=0;i<NUNIF-1;++i) //build all points but the last
    {
      //Temporary vectors for this loop
      G4DataVector* pdfi = new G4DataVector();
      G4DataVector* pdfih = new G4DataVector();
      G4DataVector* sumi = new G4DataVector();

      G4double dxi = ((*x)[i+1]-(*x)[i])/(G4double (nip-1));
      G4double pdfmax = 0;
      for (G4int k=0;k<nip;k++)
	{
	  G4double xik = (*x)[i]+k*dxi;
	  G4double pdfk = std::max(GetFSquared(mat,xik),0.);
	  pdfi->push_back(pdfk);
	  pdfmax = std::max(pdfmax,pdfk);
	  if (k < (nip-1))
	    {
	      G4double xih = xik + 0.5*dxi;
	      G4double pdfIK = std::max(GetFSquared(mat,xih),0.);
	      pdfih->push_back(pdfIK);
	      pdfmax = std::max(pdfmax,pdfIK);
	    }
	}

      //Simpson's integration
      G4double cons = dxi*0.5*(1./3.);
      sumi->push_back(0.);
      for (G4int k=1;k<nip;k++)
	{
	  G4double previous = (*sumi)[k-1];
	  G4double next = previous + cons*((*pdfi)[k-1]+4.0*(*pdfih)[k-1]+(*pdfi)[k]);
	  sumi->push_back(next);
	}

      G4double lastIntegral = (*sumi)[sumi->size()-1];
      area->push_back(lastIntegral);
      //Normalize cumulative function
      G4double factor = 1.0/lastIntegral;
      for (std::size_t k=0;k<sumi->size();++k)
	(*sumi)[k] *= factor;

      //When the PDF vanishes at one of the interval end points, its value is modified
      if ((*pdfi)[0] < 1e-35)
	(*pdfi)[0] = 1e-5*pdfmax;
      if ((*pdfi)[pdfi->size()-1] < 1e-35)
	(*pdfi)[pdfi->size()-1] = 1e-5*pdfmax;

      G4double pli = (*pdfi)[0]*factor;
      G4double pui = (*pdfi)[pdfi->size()-1]*factor;
      G4double B_temp = 1.0-1.0/(pli*pui*dx*dx);
      G4double A_temp = (1.0/(pli*dx))-1.0-B_temp;
      G4double C_temp = 1.0+A_temp+B_temp;
      if (C_temp < 1e-35)
	{
	  a->push_back(0.);
	  b->push_back(0.);
	  c->push_back(1.);
	}
      else
	{
	  a->push_back(A_temp);
	  b->push_back(B_temp);
	  c->push_back(C_temp);
	}

      //OK, now get ERR(I), the integral of the absolute difference between the rational interpolation
      //and the true pdf, extended over the interval (X(I),X(I+1))
      G4int icase = 1; //loop code
      G4bool reLoop = false;
      err->push_back(0.);
      do
	{
	  reLoop = false;
	  (*err)[i] = 0.; //zero variable
	  for (G4int k=0;k<nip;k++)
	    {
	      G4double rr = (*sumi)[k];
	      G4double pap = (*area)[i]*(1.0+((*a)[i]+(*b)[i]*rr)*rr)*(1.0+((*a)[i]+(*b)[i]*rr)*rr)/
		((1.0-(*b)[i]*rr*rr)*(*c)[i]*((*x)[i+1]-(*x)[i]));
	      if (k == 0 || k == nip-1)
		(*err)[i] += 0.5*std::fabs(pap-(*pdfi)[k]);
	      else
		(*err)[i] += std::fabs(pap-(*pdfi)[k]);
	    }
	  (*err)[i] *= dxi;

	  //If err(I) is too large, the pdf is approximated by a uniform distribution
	  if ((*err)[i] > 0.1*(*area)[i] && icase == 1)
	    {
	      (*b)[i] = 0;
	      (*a)[i] = 0;
	      (*c)[i] = 1.;
	      icase = 2;
	      reLoop = true;
	    }
	}while(reLoop);
      delete pdfi;
      delete pdfih;
      delete sumi;
    } //end of first loop over i

  //Now assign last point
  (*x)[x->size()-1] = q2max;
  a->push_back(0.);
  b->push_back(0.);
  c->push_back(0.);
  err->push_back(0.);
  area->push_back(0.);

  if (x->size() != NUNIF || a->size() != NUNIF ||
      err->size() != NUNIF || area->size() != NUNIF)
    {
      G4ExceptionDescription ed;
      ed << "Problem in building the Table for Sampling: array dimensions do not match" << G4endl;
      G4Exception("G4PenelopeRayleighModel::InitializeSamplingAlgorithm()",
		  "em2049",FatalException,ed);
    }

  /*******************************************************************************
   New grid points are added by halving the sub-intervals with the largest absolute error
  This is done up to np=150 points in the grid
  ********************************************************************************/
  do
    {
      G4double maxError = 0.0;
      std::size_t iErrMax = 0;
      for (std::size_t i=0;i<err->size()-2;++i)
	{
	  //maxError is the lagest of the interval errors err[i]
	  if ((*err)[i] > maxError)
	    {
	      maxError = (*err)[i];
	      iErrMax = i;
	    }
	}

      //OK, now I have to insert one new point in the position iErrMax
      G4double newx = 0.5*((*x)[iErrMax]+(*x)[iErrMax+1]);

      x->insert(x->begin()+iErrMax+1,newx);
      //Add place-holders in the other vectors
      area->insert(area->begin()+iErrMax+1,0.);
      a->insert(a->begin()+iErrMax+1,0.);
      b->insert(b->begin()+iErrMax+1,0.);
      c->insert(c->begin()+iErrMax+1,0.);
      err->insert(err->begin()+iErrMax+1,0.);

      //Now calculate the other parameters
      for (std::size_t i=iErrMax;i<=iErrMax+1;++i)
	{
	  //Temporary vectors for this loop
	  G4DataVector* pdfi = new G4DataVector();
	  G4DataVector* pdfih = new G4DataVector();
	  G4DataVector* sumi = new G4DataVector();

	  G4double dxLocal = (*x)[i+1]-(*x)[i];
	  G4double dxi = ((*x)[i+1]-(*x)[i])/(G4double (nip-1));
	  G4double pdfmax = 0;
	  for (G4int k=0;k<nip;k++)
	    {
	      G4double xik = (*x)[i]+k*dxi;
	      G4double pdfk = std::max(GetFSquared(mat,xik),0.);
	      pdfi->push_back(pdfk);
	      pdfmax = std::max(pdfmax,pdfk);
	      if (k < (nip-1))
		{
		  G4double xih = xik + 0.5*dxi;
		  G4double pdfIK = std::max(GetFSquared(mat,xih),0.);
		  pdfih->push_back(pdfIK);
		  pdfmax = std::max(pdfmax,pdfIK);
		}
	    }

	  //Simpson's integration
	  G4double cons = dxi*0.5*(1./3.);
	  sumi->push_back(0.);
	  for (G4int k=1;k<nip;k++)
	    {
	      G4double previous = (*sumi)[k-1];
	      G4double next = previous + cons*((*pdfi)[k-1]+4.0*(*pdfih)[k-1]+(*pdfi)[k]);
	      sumi->push_back(next);
	    }
	  G4double lastIntegral = (*sumi)[sumi->size()-1];
	  (*area)[i] = lastIntegral;

	  //Normalize cumulative function
	  G4double factor = 1.0/lastIntegral;
	  for (std::size_t k=0;k<sumi->size();++k)
	    (*sumi)[k] *= factor;

	  //When the PDF vanishes at one of the interval end points, its value is modified
	  if ((*pdfi)[0] < 1e-35)
	    (*pdfi)[0] = 1e-5*pdfmax;
	  if ((*pdfi)[pdfi->size()-1] < 1e-35)
	    (*pdfi)[pdfi->size()-1] = 1e-5*pdfmax;

	  G4double pli = (*pdfi)[0]*factor;
	  G4double pui = (*pdfi)[pdfi->size()-1]*factor;
	  G4double B_temp = 1.0-1.0/(pli*pui*dxLocal*dxLocal);
	  G4double A_temp = (1.0/(pli*dxLocal))-1.0-B_temp;
	  G4double C_temp = 1.0+A_temp+B_temp;
	  if (C_temp < 1e-35)
	    {
	      (*a)[i]= 0.;
	      (*b)[i] = 0.;
	      (*c)[i] = 1;
	    }
	  else
	    {
	      (*a)[i]= A_temp;
	      (*b)[i] = B_temp;
	      (*c)[i] = C_temp;
	    }
	  //OK, now get ERR(I), the integral of the absolute difference between the rational interpolation
	  //and the true pdf, extended over the interval (X(I),X(I+1))
	  G4int icase = 1; //loop code
	  G4bool reLoop = false;
	  do
	    {
	      reLoop = false;
	      (*err)[i] = 0.; //zero variable
	      for (G4int k=0;k<nip;k++)
		{
		  G4double rr = (*sumi)[k];
		  G4double pap = (*area)[i]*(1.0+((*a)[i]+(*b)[i]*rr)*rr)*(1.0+((*a)[i]+(*b)[i]*rr)*rr)/
		    ((1.0-(*b)[i]*rr*rr)*(*c)[i]*((*x)[i+1]-(*x)[i]));
		  if (k == 0 || k == nip-1)
		    (*err)[i] += 0.5*std::fabs(pap-(*pdfi)[k]);
		  else
		    (*err)[i] += std::fabs(pap-(*pdfi)[k]);
		}
	      (*err)[i] *= dxi;

	      //If err(I) is too large, the pdf is approximated by a uniform distribution
	      if ((*err)[i] > 0.1*(*area)[i] && icase == 1)
		{
		  (*b)[i] = 0;
		  (*a)[i] = 0;
		  (*c)[i] = 1.;
		  icase = 2;
		  reLoop = true;
		}
	    }while(reLoop);
	  delete pdfi;
	  delete pdfih;
	  delete sumi;
	}
    }while(x->size() < np);

  if (x->size() != np || a->size() != np ||
      err->size() != np || area->size() != np)
    {
      G4Exception("G4PenelopeRayleighModel::InitializeSamplingAlgorithm()",
		  "em2050",FatalException,
		  "Problem in building the extended Table for Sampling: array dimensions do not match ");
    }

  /*******************************************************************************
   Renormalization
  ********************************************************************************/
  G4double ws = 0;
  for (std::size_t i=0;i<np-1;++i)
    ws += (*area)[i];
  ws = 1.0/ws;
  G4double errMax = 0;
  for (std::size_t i=0;i<np-1;++i)
    {
      (*area)[i] *= ws;
      (*err)[i] *= ws;
      errMax = std::max(errMax,(*err)[i]);
    }

  //Vector with the normalized cumulative distribution
  G4DataVector* PAC = new G4DataVector();
  PAC->push_back(0.);
  for (std::size_t i=0;i<np-1;++i)
    {
      G4double previous = (*PAC)[i];
      PAC->push_back(previous+(*area)[i]);
    }
  (*PAC)[PAC->size()-1] = 1.;

  /*******************************************************************************
  Pre-calculated limits for the initial binary search for subsequent sampling
  ********************************************************************************/
  std::vector<std::size_t> *ITTL = new std::vector<std::size_t>;
  std::vector<std::size_t> *ITTU = new std::vector<std::size_t>;

  //Just create place-holders
  for (std::size_t i=0;i<np;++i)
    {
      ITTL->push_back(0);
      ITTU->push_back(0);
    }

  G4double bin = 1.0/(np-1);
  (*ITTL)[0]=0;
  for (std::size_t i=1;i<(np-1);++i)
    {
      G4double ptst = i*bin;
      G4bool found = false;
      for (std::size_t j=(*ITTL)[i-1];j<np && !found;++j)
	{
	  if ((*PAC)[j] > ptst)
	    {
	      (*ITTL)[i] = j-1;
	      (*ITTU)[i-1] = j;
	      found = true;
	    }
	}
    }
  (*ITTU)[ITTU->size()-2] = ITTU->size()-1;
  (*ITTU)[ITTU->size()-1] = ITTU->size()-1;
  (*ITTL)[ITTL->size()-1] = ITTU->size()-2;

  if (ITTU->size() != np || ITTU->size() != np)
    {
      G4Exception("G4PenelopeRayleighModel::InitializeSamplingAlgorithm()",
		  "em2051",FatalException,
		  "Problem in building the Limit Tables for Sampling: array dimensions do not match");
    }

  /********************************************************************************
    Copy tables
  ********************************************************************************/
  G4PenelopeSamplingData* theTable = new G4PenelopeSamplingData(np);
  for (std::size_t i=0;i<np;++i)
    {
      theTable->AddPoint((*x)[i],(*PAC)[i],(*a)[i],(*b)[i],(*ITTL)[i],(*ITTU)[i]);
    }

  if (fVerboseLevel > 2)
    {
      G4cout << "*************************************************************************" <<
	G4endl;
      G4cout << "Sampling table for Penelope Rayleigh scattering in " << mat->GetName() << G4endl;
      theTable->DumpTable();
    }
  fSamplingTable->insert(std::make_pair(mat,theTable));

  //Clean up temporary vectors
  delete x;
  delete a;
  delete b;
  delete c;
  delete err;
  delete area;
  delete PAC;
  delete ITTL;
  delete ITTU;

  //DONE!
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::GetPMaxTable(const G4Material* mat)
{
  if (!fPMaxTable)
    {
      G4cout << "G4PenelopeRayleighModel::BuildPMaxTable" << G4endl;
      G4cout << "Going to instanziate the fPMaxTable !" << G4endl;
      G4cout << "That should _not_ be here! " << G4endl;
      fPMaxTable = new std::map<const G4Material*,G4PhysicsFreeVector*>;
    }
  //check if the table is already there
  if (fPMaxTable->count(mat))
    return;

  //otherwise build it
  if (!fSamplingTable)
    {
      G4Exception("G4PenelopeRayleighModel::GetPMaxTable()",
		  "em2052",FatalException,
		  "SamplingTable is not properly instantiated");
      return;
    }

  //This should not be: the sampling table is built before the p-table
  if (!fSamplingTable->count(mat))
    {
       G4ExceptionDescription ed;
       ed << "Sampling table for material " << mat->GetName() << " not found";
       G4Exception("G4PenelopeRayleighModel::GetPMaxTable()",
                  "em2052",FatalException,
                  ed);
       return;
    }

  G4PenelopeSamplingData *theTable = fSamplingTable->find(mat)->second;
  std::size_t tablePoints = theTable->GetNumberOfStoredPoints();

  std::size_t nOfEnergyPoints = fLogEnergyGridPMax.size();
  G4PhysicsFreeVector* theVec = new G4PhysicsFreeVector(nOfEnergyPoints);

  const std::size_t nip = 51; //hard-coded in Penelope

  for (std::size_t ie=0;ie<fLogEnergyGridPMax.size();++ie)
    {
      G4double energy = G4Exp(fLogEnergyGridPMax[ie]);
      G4double Qm = 2.0*energy/electron_mass_c2; //this is non-dimensional now
      G4double Qm2 = Qm*Qm;
      G4double firstQ2 = theTable->GetX(0);
      G4double lastQ2 = theTable->GetX(tablePoints-1);
      G4double thePMax = 0;

      if (Qm2 > firstQ2)
	{
	  if (Qm2 < lastQ2)
	    {
	      //bisection to look for the index of Qm
	      std::size_t lowerBound = 0;
	      std::size_t upperBound = tablePoints-1;
	      while (lowerBound <= upperBound)
		{
		  std::size_t midBin = (lowerBound + upperBound)/2;
		  if( Qm2 < theTable->GetX(midBin))
		    { upperBound = midBin-1; }
		  else
		    { lowerBound = midBin+1; }
		}
	      //upperBound is the output (but also lowerBounf --> should be the same!)
	      G4double Q1 = theTable->GetX(upperBound);
 	      G4double Q2 = Qm2;
	      G4double DQ = (Q2-Q1)/((G4double)(nip-1));
	      G4double theA = theTable->GetA(upperBound);
	      G4double theB = theTable->GetB(upperBound);
	      G4double thePAC = theTable->GetPAC(upperBound);
	      G4DataVector* fun = new G4DataVector();
	      for (std::size_t k=0;k<nip;++k)
		{
		  G4double qi = Q1 + k*DQ;
		  G4double tau = (qi-Q1)/
		    (theTable->GetX(upperBound+1)-Q1);
		  G4double con1 = 2.0*theB*tau;
		  G4double ci = 1.0+theA+theB;
		  G4double con2 = ci-theA*tau;
		  G4double etap = 0;
		  if (std::fabs(con1) > 1.0e-16*std::fabs(con2))
		    etap = con2*(1.0-std::sqrt(1.0-2.0*tau*con1/(con2*con2)))/con1;
		  else
		    etap = tau/con2;
		  G4double theFun = (theTable->GetPAC(upperBound+1)-thePAC)*
		    (1.0+(theA+theB*etap)*etap)*(1.0+(theA+theB*etap)*etap)/
		    ((1.0-theB*etap*etap)*ci*(theTable->GetX(upperBound+1)-Q1));
		  fun->push_back(theFun);
		}
	      //Now intergrate numerically the fun Cavalieri-Simpson's method
	      G4DataVector* sum = new G4DataVector;
	      G4double CONS = DQ*(1./12.);
	      G4double HCONS = 0.5*CONS;
	      sum->push_back(0.);
	      G4double secondPoint = (*sum)[0] +
		(5.0*(*fun)[0]+8.0*(*fun)[1]-(*fun)[2])*CONS;
	      sum->push_back(secondPoint);
	      for (std::size_t hh=2;hh<nip-1;++hh)
		{
		  G4double previous = (*sum)[hh-1];
		  G4double next = previous+(13.0*((*fun)[hh-1]+(*fun)[hh])-
					    (*fun)[hh+1]-(*fun)[hh-2])*HCONS;
		  sum->push_back(next);
		}
	      G4double last = (*sum)[nip-2]+(5.0*(*fun)[nip-1]+8.0*(*fun)[nip-2]-
					     (*fun)[nip-3])*CONS;
	      sum->push_back(last);
	      thePMax = thePAC + (*sum)[sum->size()-1]; //last point
	      delete fun;
	      delete sum;
	    }
	  else
	    {
	      thePMax = 1.0;
	    }
	}
      else
	{
	  thePMax = theTable->GetPAC(0);
	}

      //Write number in the table
      theVec->PutValue(ie,energy,thePMax);
  }

  fPMaxTable->insert(std::make_pair(mat,theVec));
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeRayleighModel::DumpFormFactorTable(const G4Material* mat)
{
  G4cout << "*****************************************************************" << G4endl;
  G4cout << "G4PenelopeRayleighModel: Form Factor Table for " << mat->GetName() << G4endl;
  //try to use the same format as Penelope-Fortran, namely Q (/m_e*c) and F
  G4cout <<  "Q/(m_e*c)                 F(Q)     " << G4endl;
  G4cout << "*****************************************************************" << G4endl;
  if (!fLogFormFactorTable->count(mat))
    BuildFormFactorTable(mat);

  G4PhysicsFreeVector* theVec = fLogFormFactorTable->find(mat)->second;
  for (std::size_t i=0;i<theVec->GetVectorLength();++i)
    {
      G4double logQ2 = theVec->GetLowEdgeEnergy(i);
      G4double Q = G4Exp(0.5*logQ2);
      G4double logF2 = (*theVec)[i];
      G4double F = G4Exp(0.5*logF2);
      G4cout << Q << "              " << F << G4endl;
    }
  //DONE
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeRayleighModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;
  }
}
