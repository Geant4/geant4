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
// 23 Nov 2010   L Pandola    First complete implementation, Penelope v2008
// 24 May 2011   L. Pandola   Renamed (make default Penelope)
// 13 Mar 2012   L. Pandola   Updated the interface for the angular generator
// 18 Jul 2012   L. Pandola   Migrate to the new interface of the angular generator, which
//                            now provides the G4ThreeVector and takes care of rotation
// 02 Oct 2013   L. Pandola   Migrated to MT
// 17 Oct 2013   L. Pandola   Partially revert the MT migration: the angular generator is
//                             kept as thread-local, and created/managed by the workers.
//

#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PenelopeBremsstrahlungFS.hh"
#include "G4PenelopeBremsstrahlungAngular.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PenelopeOscillatorManager.hh"
#include "G4PenelopeCrossSection.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
namespace { G4Mutex  PenelopeBremsstrahlungModelMutex = G4MUTEX_INITIALIZER; }

G4PenelopeBremsstrahlungModel::G4PenelopeBremsstrahlungModel(const G4ParticleDefinition* part,
							     const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),fParticle(nullptr),
   fPenelopeFSHelper(nullptr),fPenelopeAngular(nullptr),fEnergyGrid(nullptr),
   fXSTableElectron(nullptr),fXSTablePositron(nullptr),
   fIsInitialised(false),fLocalTable(false)

{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  nBins = 200;

  if (part)
    SetParticle(part);

  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
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

G4PenelopeBremsstrahlungModel::~G4PenelopeBremsstrahlungModel()
{
  if (IsMaster() || fLocalTable)
    {
      ClearTables();
      if (fPenelopeFSHelper)
	delete fPenelopeFSHelper;
    }
  // This is thread-local at the moment
  if (fPenelopeAngular)
    delete fPenelopeAngular;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::Initialise(const G4ParticleDefinition* part,
                                             const G4DataVector& theCuts)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling G4PenelopeBremsstrahlungModel::Initialise()" << G4endl;

  SetParticle(part);

  if (IsMaster() && part == fParticle)
    {
      if (!fPenelopeFSHelper)
	fPenelopeFSHelper = new G4PenelopeBremsstrahlungFS(fVerboseLevel);
      if (!fPenelopeAngular)
	fPenelopeAngular = new G4PenelopeBremsstrahlungAngular();
      //Clear and re-build the tables
      ClearTables();

      //forces the cleaning of tables, in this specific case
      if (fPenelopeAngular)
	fPenelopeAngular->Initialize();

      //Set the number of bins for the tables. 20 points per decade
      nBins = (std::size_t) (20*std::log10(HighEnergyLimit()/LowEnergyLimit()));
      nBins = std::max(nBins,(std::size_t)100);
      fEnergyGrid = new G4PhysicsLogVector(LowEnergyLimit(),
                                      HighEnergyLimit(),
                                      nBins-1); //one hidden bin is added

      fXSTableElectron = new
	std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;
      fXSTablePositron = new
	std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;

      G4ProductionCutsTable* theCoupleTable =
	G4ProductionCutsTable::GetProductionCutsTable();

      //Build tables for all materials
      for (G4int i=0;i<(G4int)theCoupleTable->GetTableSize();++i)
	{
	  const G4Material* theMat =
	    theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
	  //Forces the building of the helper tables
	  fPenelopeFSHelper->BuildScaledXSTable(theMat,theCuts.at(i),IsMaster());
	  fPenelopeAngular->PrepareTables(theMat,IsMaster());
	  BuildXSTable(theMat,theCuts.at(i));

	}

      if (fVerboseLevel > 2) {
	G4cout << "Penelope Bremsstrahlung model v2008 is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / keV << " keV - "
	       << HighEnergyLimit() / GeV << " GeV."
	       << G4endl;
      }
    }

  if(fIsInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  fIsInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::InitialiseLocal(const G4ParticleDefinition* part,
						    G4VEmModel *masterModel)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopeBremsstrahlungModel::InitialiseLocal()" << G4endl;
  //
  //Check that particle matches: one might have multiple master models (e.g.
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeBremsstrahlungModel* theModel =
	static_cast<G4PenelopeBremsstrahlungModel*> (masterModel);

      //Copy pointers to the data tables
      fEnergyGrid = theModel->fEnergyGrid;
      fXSTableElectron = theModel->fXSTableElectron;
      fXSTablePositron = theModel->fXSTablePositron;
      fPenelopeFSHelper = theModel->fPenelopeFSHelper;

      //created in each thread and initialized.
      if (!fPenelopeAngular)
	fPenelopeAngular = new G4PenelopeBremsstrahlungAngular();
      //forces the cleaning of tables, in this specific case
      if (fPenelopeAngular)
	fPenelopeAngular->Initialize();

      G4ProductionCutsTable* theCoupleTable =
	G4ProductionCutsTable::GetProductionCutsTable();
      //Build tables for all materials
      for (G4int i=0;i<(G4int)theCoupleTable->GetTableSize();++i)
	{
	  const G4Material* theMat =
	    theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial();
	  fPenelopeAngular->PrepareTables(theMat,IsMaster());
	}

      //copy the data
      nBins = theModel->nBins;

      //Same verbosity for all workers, as the master
      fVerboseLevel = theModel->fVerboseLevel;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::CrossSectionPerVolume(const G4Material* material,
                                           const G4ParticleDefinition* theParticle,
                                           G4double energy,
                                           G4double cutEnergy,
                                           G4double)
{
  //
  if (fVerboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4PenelopeBremsstrahlungModel" << G4endl;

  SetupForMaterial(theParticle, material, energy);
  G4double crossPerMolecule = 0.;

  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(theParticle,material,
                                                                cutEnergy);
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
    G4cout << "G4PenelopeBremsstrahlungModel " << G4endl;
    G4cout << "Mean free path for gamma emission > " << cutEnergy/keV << " keV at " <<
      energy/keV << " keV = " << (1./crossPerVolume)/mm << " mm" << G4endl;
  }

  return crossPerVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//This is a dummy method. Never inkoved by the tracking, it just issues
//a warning if one tries to get Cross Sections per Atom via the
//G4EmCalculator.
G4double G4PenelopeBremsstrahlungModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
								     G4double,
								     G4double,
								     G4double,
								     G4double,
								     G4double)
{
  G4cout << "*** G4PenelopeBremsstrahlungModel -- WARNING ***" << G4endl;
  G4cout << "Penelope Bremsstrahlung model v2008 does not calculate cross section _per atom_ " << G4endl;
  G4cout << "so the result is always zero. For physics values, please invoke " << G4endl;
  G4cout << "GetCrossSectionPerVolume() or GetMeanFreePath() via the G4EmCalculator" << G4endl;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::ComputeDEDXPerVolume(const G4Material* material,
							     const G4ParticleDefinition* theParticle,
							     G4double kineticEnergy,
							     G4double cutEnergy)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling ComputeDEDX() of G4PenelopeBremsstrahlungModel" << G4endl;

  G4PenelopeCrossSection* theXS = GetCrossSectionTableForCouple(theParticle,material,
                                                                cutEnergy);
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
      G4cout << "G4PenelopeBremsstrahlungModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " <<
        kineticEnergy/keV << " keV = " <<
        sPowerPerVolume/(keV/mm) << " keV/mm" << G4endl;
    }
  return sPowerPerVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::SampleSecondaries(std::vector<G4DynamicParticle*>*fvect,
						      const G4MaterialCutsCouple* couple,
						      const G4DynamicParticle* aDynamicParticle,
						      G4double cutG,
						      G4double)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4PenelopeBremsstrahlungModel" << G4endl;

  G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
  const G4Material* material = couple->GetMaterial();

  if (kineticEnergy <= fIntrinsicLowEnergyLimit)
   {
     fParticleChange->SetProposedKineticEnergy(0.);
     fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy);
     return ;
   }

  G4ParticleMomentum particleDirection0 = aDynamicParticle->GetMomentumDirection();
  //This is the momentum
  G4ThreeVector initialMomentum =  aDynamicParticle->GetMomentum();

  //Not enough energy to produce a secondary! Return with nothing happened
  if (kineticEnergy < cutG)
    return;

  if (fVerboseLevel > 3)
    G4cout << "Going to sample gamma energy for: " <<material->GetName() << " " <<
      "energy = " << kineticEnergy/keV << ", cut = " << cutG/keV << G4endl;

   //Sample gamma's energy according to the spectrum
  G4double gammaEnergy =
    fPenelopeFSHelper->SampleGammaEnergy(kineticEnergy,material,cutG);

  if (fVerboseLevel > 3)
    G4cout << "Sampled gamma energy: " << gammaEnergy/keV << " keV" << G4endl;

  //Now sample the direction for the Gamma. Notice that the rotation is already done
  //Z is unused here, I plug 0. The information is in the material pointer
   G4ThreeVector gammaDirection1 =
     fPenelopeAngular->SampleDirection(aDynamicParticle,gammaEnergy,0,material);

  if (fVerboseLevel > 3)
    G4cout << "Sampled cosTheta for e-: " << gammaDirection1.cosTheta() << G4endl;

  G4double residualPrimaryEnergy = kineticEnergy-gammaEnergy;
  if (residualPrimaryEnergy < 0)
    {
      //Ok we have a problem, all energy goes with the gamma
      gammaEnergy += residualPrimaryEnergy;
      residualPrimaryEnergy = 0.0;
    }

  //Produce final state according to momentum conservation
  G4ThreeVector particleDirection1 = initialMomentum - gammaEnergy*gammaDirection1;
  particleDirection1 = particleDirection1.unit(); //normalize

  //Update the primary particle
  if (residualPrimaryEnergy > 0.)
    {
      fParticleChange->ProposeMomentumDirection(particleDirection1);
      fParticleChange->SetProposedKineticEnergy(residualPrimaryEnergy);
    }
  else
    {
      fParticleChange->SetProposedKineticEnergy(0.);
    }

  //Now produce the photon
  G4DynamicParticle* theGamma = new G4DynamicParticle(G4Gamma::Gamma(),
                                                      gammaDirection1,
                                                      gammaEnergy);
  fvect->push_back(theGamma);

  if (fVerboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeBremsstrahlung" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << residualPrimaryEnergy/keV << " keV" << G4endl;
      G4cout << "Bremsstrahlung photon " << gammaEnergy/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (residualPrimaryEnergy+gammaEnergy)/keV
             << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }

  if (fVerboseLevel > 0)
    {
      G4double energyDiff = std::fabs(residualPrimaryEnergy+gammaEnergy-kineticEnergy);
      if (energyDiff > 0.05*keV)
        G4cout << "Warning from G4PenelopeBremsstrahlung: problem with energy conservation: "
	       <<
          (residualPrimaryEnergy+gammaEnergy)/keV <<
          " keV (final) vs. " <<
          kineticEnergy/keV << " keV (initial)" << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::ClearTables()
{
  if (!IsMaster() && !fLocalTable)
    //Should not be here!
    G4Exception("G4PenelopeBremsstrahlungModel::ClearTables()",
		"em0100",FatalException,"Worker thread in this method");

  if (fXSTableElectron)
    {
      for (auto& item : (*fXSTableElectron))        
	delete item.second;        
      delete fXSTableElectron;
      fXSTableElectron = nullptr;
    }
  if (fXSTablePositron)
    {
      for (auto& item : (*fXSTablePositron))                
	delete item.second;    
      delete fXSTablePositron;
      fXSTablePositron = nullptr;
    }
  /*
  if (fEnergyGrid)
    delete fEnergyGrid;
  */
  if (fPenelopeFSHelper)
    fPenelopeFSHelper->ClearTables(IsMaster());

  if (fVerboseLevel > 2)
    G4cout << "G4PenelopeBremsstrahlungModel: cleared tables" << G4endl;

 return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::MinEnergyCut(const G4ParticleDefinition*,
						       const G4MaterialCutsCouple*)
{
  return fIntrinsicLowEnergyLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungModel::BuildXSTable(const G4Material* mat,G4double cut)
{
  if (!IsMaster() && !fLocalTable)
    //Should not be here!
    G4Exception("G4PenelopeBremsstrahlungModel::BuildXSTable()",
		"em0100",FatalException,"Worker thread in this method");

  //The key of the map
  std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);

  //The table already exists
  if (fXSTableElectron->count(theKey) && fXSTablePositron->count(theKey))
    return;

  //
  //This method fills the G4PenelopeCrossSection containers for electrons or positrons
  //and for the given material/cut couple.
  //Equivalent of subroutines EBRaT and PINaT of Penelope
  //
  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopeBremsstrahlungModel: going to build cross section table " << G4endl;
      G4cout << "for e+/e- in " << mat->GetName() << " for Ecut(gamma)= " <<
	cut/keV << " keV " << G4endl;
    }

  //Tables have been already created (checked by GetCrossSectionTableForCouple)
  if (fEnergyGrid->GetVectorLength() != nBins)
    {
      G4ExceptionDescription ed;
      ed << "Energy Grid looks not initialized" << G4endl;
      ed << nBins << " " << fEnergyGrid->GetVectorLength() << G4endl;
      G4Exception("G4PenelopeBremsstrahlungModel::BuildXSTable()",
		  "em2016",FatalException,ed);
    }

  G4PenelopeCrossSection* XSEntry = new G4PenelopeCrossSection(nBins);
  G4PenelopeCrossSection* XSEntry2 = new G4PenelopeCrossSection(nBins);
  const G4PhysicsTable* table = fPenelopeFSHelper->GetScaledXSTable(mat,cut);

  //loop on the energy grid
  for (std::size_t bin=0;bin<nBins;++bin)
    {
       G4double energy = fEnergyGrid->GetLowEdgeEnergy(bin);
       G4double XH0=0, XH1=0, XH2=0;
       G4double XS0=0, XS1=0, XS2=0;

       //Global xs factor
       G4double fact = fPenelopeFSHelper->GetEffectiveZSquared(mat)*
	 ((energy+electron_mass_c2)*(energy+electron_mass_c2)/
	  (energy*(energy+2.0*electron_mass_c2)));

       G4double restrictedCut = cut/energy;

       //Now I need the dSigma/dX profile - interpolated on energy - for
       //the 32-point x grid. Interpolation is log-log
       std::size_t nBinsX = fPenelopeFSHelper->GetNBinsX();
       G4double* tempData = new G4double[nBinsX];
       G4double logene = G4Log(energy);
       for (std::size_t ix=0;ix<nBinsX;++ix)
	 {
	   //find dSigma/dx for the given E. X belongs to the 32-point grid.
	   G4double val = (*table)[ix]->Value(logene);
	   tempData[ix] = G4Exp(val); //back to the real value!
	 }

       G4double XH0A = 0.;
       if (restrictedCut <= 1) //calculate only if we are above threshold!
	 XH0A = fPenelopeFSHelper->GetMomentumIntegral(tempData,1.0,-1) -
	   fPenelopeFSHelper->GetMomentumIntegral(tempData,restrictedCut,-1);
       G4double XS1A = fPenelopeFSHelper->GetMomentumIntegral(tempData,
							      restrictedCut,0);
       G4double XS2A = fPenelopeFSHelper->GetMomentumIntegral(tempData,
							      restrictedCut,1);
       G4double XH1A=0, XH2A=0;
       if (restrictedCut <=1)
	 {
	   XH1A = fPenelopeFSHelper->GetMomentumIntegral(tempData,1.0,0) -
	     XS1A;
	   XH2A = fPenelopeFSHelper->GetMomentumIntegral(tempData,1.0,1) -
	     XS2A;
	 }
       delete[] tempData;

       XH0 = XH0A*fact;
       XS1 = XS1A*fact*energy;
       XH1 = XH1A*fact*energy;
       XS2 = XS2A*fact*energy*energy;
       XH2 = XH2A*fact*energy*energy;

       XSEntry->AddCrossSectionPoint(bin,energy,XH0,XH1,XH2,XS0,XS1,XS2);

       //take care of positrons
       G4double posCorrection = GetPositronXSCorrection(mat,energy);
       XSEntry2->AddCrossSectionPoint(bin,energy,XH0*posCorrection,
				      XH1*posCorrection,
				      XH2*posCorrection,
				      XS0,
				      XS1*posCorrection,
				      XS2*posCorrection);
    }

  //Insert in the appropriate table
  fXSTableElectron->insert(std::make_pair(theKey,XSEntry));
  fXSTablePositron->insert(std::make_pair(theKey,XSEntry2));

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeCrossSection*
G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple(const G4ParticleDefinition* part,
							     const G4Material* mat,
							     G4double cut)
{
  if (part != G4Electron::Electron() && part != G4Positron::Positron())
    {
      G4ExceptionDescription ed;
      ed << "Invalid particle: " << part->GetParticleName() << G4endl;
      G4Exception("G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple()",
		  "em0001",FatalException,ed);
      return nullptr;
    }

  if (part == G4Electron::Electron())
    {
      //Either Initialize() was not called, or we are in a slave and InitializeLocal() was
      //not invoked
      if (!fXSTableElectron)
        {
	  //create a **thread-local** version of the table. Used only for G4EmCalculator and
	  //Unit Tests
          G4String excep = "The Cross Section Table for e- was not initialized correctly!";
          G4Exception("G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple()",
		      "em2013",JustWarning,excep);
	  fLocalTable = true;
	  fXSTableElectron = new
	    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;
	  if (!fEnergyGrid)
	    fEnergyGrid = new G4PhysicsLogVector(LowEnergyLimit(),
						HighEnergyLimit(),
						nBins-1); //one hidden bin is added
	  if (!fPenelopeFSHelper)
 	    fPenelopeFSHelper = new G4PenelopeBremsstrahlungFS(fVerboseLevel);
        }
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (fXSTableElectron->count(theKey)) //table already built
        return fXSTableElectron->find(theKey)->second;
      else
	{
	  //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
	  //not filled up. This can happen in a UnitTest or via G4EmCalculator
	  if (fVerboseLevel > 0)
	    {
	      //G4Exception (warning) is issued only in verbose mode
	      G4ExceptionDescription ed;
	      ed << "Unable to find e- table for " << mat->GetName() << " at Ecut(gamma)= "
		 << cut/keV << " keV " << G4endl;
	      ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple()",
			  "em2009",JustWarning,ed);
	    }
	  //protect file reading via autolock
	  G4AutoLock lock(&PenelopeBremsstrahlungModelMutex);
          fPenelopeFSHelper->BuildScaledXSTable(mat,cut,true); //pretend to be a master
	  BuildXSTable(mat,cut);
	  lock.unlock();
	  //now it should be ok
	  return fXSTableElectron->find(theKey)->second;
	}
    }
  if (part == G4Positron::Positron())
    {
      //Either Initialize() was not called, or we are in a slave and InitializeLocal() was
      //not invoked
      if (!fXSTablePositron)
        {
	  G4String excep = "The Cross Section Table for e+ was not initialized correctly!";
          G4Exception("G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple()",
		      "em2013",JustWarning,excep);
	  fLocalTable = true;
	  fXSTablePositron = new
	    std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*>;
	  if (!fEnergyGrid)
	    fEnergyGrid = new G4PhysicsLogVector(LowEnergyLimit(),
						HighEnergyLimit(),
						nBins-1); //one hidden bin is added
	  if (!fPenelopeFSHelper)
            fPenelopeFSHelper = new G4PenelopeBremsstrahlungFS(fVerboseLevel);
        }
      std::pair<const G4Material*,G4double> theKey = std::make_pair(mat,cut);
      if (fXSTablePositron->count(theKey)) //table already built
        return fXSTablePositron->find(theKey)->second;
      else
        {
	  //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
	  //not filled up. This can happen in a UnitTest or via G4EmCalculator
	  if (fVerboseLevel > 0)
	    {
	      //Issue a G4Exception (warning) only in verbose mode
	      G4ExceptionDescription ed;
	      ed << "Unable to find e+ table for " << mat->GetName() << " at Ecut(gamma)= "
		 << cut/keV << " keV " << G4endl;
	      ed << "This can happen only in Unit Tests or via G4EmCalculator" << G4endl;
	      G4Exception("G4PenelopeBremsstrahlungModel::GetCrossSectionTableForCouple()",
			  "em2009",JustWarning,ed);
	    }
	  //protect file reading via autolock
	  G4AutoLock lock(&PenelopeBremsstrahlungModelMutex);
          fPenelopeFSHelper->BuildScaledXSTable(mat,cut,true); //pretend to be a master
	  BuildXSTable(mat,cut);
	  lock.unlock();
	  //now it should be ok
	  return fXSTablePositron->find(theKey)->second;
        }
    }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungModel::GetPositronXSCorrection(const G4Material* mat,
								  G4double energy)
{
  //The electron-to-positron correction factor is set equal to the ratio of the
  //radiative stopping powers for positrons and electrons, which has been calculated
  //by Kim et al. (1986) (cf. Berger and Seltzer, 1982). Here, it is used an
  //analytical approximation which reproduces the tabulated values with 0.5%
  //accuracy
  G4double t=G4Log(1.0+1e6*energy/
		      (electron_mass_c2*fPenelopeFSHelper->GetEffectiveZSquared(mat)));
  G4double corr = 1.0-G4Exp(-t*(1.2359e-1-t*(6.1274e-2-t*
					   (3.1516e-2-t*(7.7446e-3-t*(1.0595e-3-t*
								      (7.0568e-5-t*
								       1.8080e-6)))))));
  return corr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeBremsstrahlungModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;
  }
}
