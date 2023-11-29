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
// 13 Jan 2010   L Pandola    First implementation (updated to Penelope08)
// 24 May 2011   L Pandola    Renamed (make v2008 as default Penelope)
// 18 Sep 2013   L Pandola    Migration to MT paradigm. Only master model deals with
//                             data and creates shared tables
//

#include "G4PenelopeGammaConversionModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4AutoLock.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
const G4int G4PenelopeGammaConversionModel::fMaxZ;
G4PhysicsFreeVector* G4PenelopeGammaConversionModel::fLogAtomicCrossSection[] = {nullptr};
G4double G4PenelopeGammaConversionModel::fAtomicScreeningRadius[] = {0.,  //pad a zero, so to use fAtomicScreeningRadius[Z]
								     1.2281e+02,7.3167e+01,6.9228e+01,6.7301e+01,
								     6.4696e+01,6.1228e+01,5.7524e+01,5.4033e+01,
								     5.0787e+01,4.7851e+01,4.6373e+01,4.5401e+01,
								     4.4503e+01,4.3815e+01,4.3074e+01,4.2321e+01,
								     4.1586e+01,4.0953e+01,4.0524e+01,4.0256e+01,
								     3.9756e+01,3.9144e+01,3.8462e+01,3.7778e+01,
								     3.7174e+01,3.6663e+01,3.5986e+01,3.5317e+01,
								     3.4688e+01,3.4197e+01,3.3786e+01,3.3422e+01,
								     3.3068e+01,3.2740e+01,3.2438e+01,3.2143e+01,
								     3.1884e+01,3.1622e+01,3.1438e+01,3.1142e+01,
								     3.0950e+01,3.0758e+01,3.0561e+01,3.0285e+01,
								     3.0097e+01,2.9832e+01,2.9581e+01,2.9411e+01,
								     2.9247e+01,2.9085e+01,2.8930e+01,2.8721e+01,
								     2.8580e+01,2.8442e+01,2.8312e+01,2.8139e+01,
								     2.7973e+01,2.7819e+01,2.7675e+01,2.7496e+01,
								     2.7285e+01,2.7093e+01,2.6911e+01,2.6705e+01,
								     2.6516e+01,2.6304e+01,2.6108e+01,2.5929e+01,
								     2.5730e+01,2.5577e+01,2.5403e+01,2.5245e+01,
								     2.5100e+01,2.4941e+01,2.4790e+01,2.4655e+01,
								     2.4506e+01,2.4391e+01,2.4262e+01,2.4145e+01,
								     2.4039e+01,2.3922e+01,2.3813e+01,2.3712e+01,
								     2.3621e+01,2.3523e+01,2.3430e+01,2.3331e+01,
								     2.3238e+01,2.3139e+01,2.3048e+01,2.2967e+01,
								     2.2833e+01,2.2694e+01,2.2624e+01,2.2545e+01,
								     2.2446e+01,2.2358e+01,2.2264e+01};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeGammaConversionModel::G4PenelopeGammaConversionModel(const G4ParticleDefinition* part,
							       const G4String& nam)
  :G4VEmModel(nam),fParticleChange(nullptr),fParticle(nullptr),
   fEffectiveCharge(nullptr),fMaterialInvScreeningRadius(nullptr),
   fScreeningFunction(nullptr),fIsInitialised(false),fLocalTable(false)
{
  fIntrinsicLowEnergyLimit = 2.0*electron_mass_c2;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  fSmallEnergy = 1.1*MeV;

  if (part)
    SetParticle(part);

  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  fVerboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeGammaConversionModel::~G4PenelopeGammaConversionModel()
{
  //Delete shared tables, they exist only in the master model
  if (IsMaster() || fLocalTable)
    {
      for(G4int i=0; i<=fMaxZ; ++i) 
	{
	  if(fLogAtomicCrossSection[i]) { 
	    delete fLogAtomicCrossSection[i];
	    fLogAtomicCrossSection[i] = nullptr;
	  }
	}
      if (fEffectiveCharge)
	delete fEffectiveCharge;
      if (fMaterialInvScreeningRadius)
	delete fMaterialInvScreeningRadius;
      if (fScreeningFunction)
	delete fScreeningFunction;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeGammaConversionModel::Initialise(const G4ParticleDefinition* part,
						const G4DataVector&)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopeGammaConversionModel::Initialise()" << G4endl;

  SetParticle(part);

  //Only the master model creates/fills/destroys the tables
  if (IsMaster() && part == fParticle)
    {
      //delete old material data...
      if (fEffectiveCharge)
	{
	  delete fEffectiveCharge;
	  fEffectiveCharge = nullptr;
	}
      if (fMaterialInvScreeningRadius)
	{
	  delete fMaterialInvScreeningRadius;
	  fMaterialInvScreeningRadius = nullptr;
	}
      if (fScreeningFunction)
	{
	  delete fScreeningFunction;
	  fScreeningFunction = nullptr;
	}
      //and create new ones
      fEffectiveCharge = new std::map<const G4Material*,G4double>;
      fMaterialInvScreeningRadius = new std::map<const G4Material*,G4double>;
      fScreeningFunction = new std::map<const G4Material*,std::pair<G4double,G4double> >;

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
	      if (iZ <= fMaxZ &&  !fLogAtomicCrossSection[iZ])		
		ReadDataFile(iZ);
	    }

	  //check if material data are available
	  if (!fEffectiveCharge->count(material))
	    InitializeScreeningFunctions(material);
	}
      if (fVerboseLevel > 0) {
	G4cout << "Penelope Gamma Conversion model v2008 is initialized " << G4endl
	       << "Energy range: "
	       << LowEnergyLimit() / MeV << " MeV - "
	       << HighEnergyLimit() / GeV << " GeV"
	       << G4endl;
      }
    }
  if(fIsInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  fIsInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeGammaConversionModel::InitialiseLocal(const G4ParticleDefinition* part,
						     G4VEmModel *masterModel)
{
  if (fVerboseLevel > 3)
    G4cout << "Calling  G4PenelopeGammaConversionModel::InitialiseLocal()" << G4endl;
  //
  //Check that particle matches: one might have multiple master models (e.g.
  //for e+ and e-).
  //
  if (part == fParticle)
    {
      //Get the const table pointers from the master to the workers
      const G4PenelopeGammaConversionModel* theModel =
	static_cast<G4PenelopeGammaConversionModel*> (masterModel);

      //Copy pointers to the data tables
      fEffectiveCharge = theModel->fEffectiveCharge;
      fMaterialInvScreeningRadius = theModel->fMaterialInvScreeningRadius;
      fScreeningFunction = theModel->fScreeningFunction;      
      for(G4int i=0; i<=fMaxZ; ++i) 
	fLogAtomicCrossSection[i] = theModel->fLogAtomicCrossSection[i];

      //Same verbosity for all workers, as the master
      fVerboseLevel = theModel->fVerboseLevel;
    }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
namespace { G4Mutex  PenelopeGammaConversionModelMutex = G4MUTEX_INITIALIZER; }

G4double G4PenelopeGammaConversionModel::ComputeCrossSectionPerAtom(
								    const G4ParticleDefinition*,
								    G4double energy,
								    G4double Z, G4double,
								    G4double, G4double)
{
  //
  // Penelope model v2008.
  // Cross section (including triplet production) read from database and managed
  // through the G4CrossSectionHandler utility. Cross section data are from
  // M.J. Berger and J.H. Hubbel (XCOM), Report NBSIR 887-3598
  //

  if (energy < fIntrinsicLowEnergyLimit)
    return 0;

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
	   G4Exception("G4PenelopeGammaConversionModel::ComputeCrossSectionPerAtom()",
		       "em2018",JustWarning,ed);
	 }
       //protect file reading via autolock
       G4AutoLock lock(&PenelopeGammaConversionModelMutex);
       ReadDataFile(iZ);
       lock.unlock();
       fLocalTable = true;
     }
  G4double cs = 0;
  G4double logene = G4Log(energy);
  G4PhysicsFreeVector* theVec = fLogAtomicCrossSection[iZ];
  G4double logXS = theVec->Value(logene);
  cs = G4Exp(logXS);

  if (fVerboseLevel > 2)
    G4cout << "Gamma conversion cross section at " << energy/MeV << " MeV for Z=" << Z <<
      " = " << cs/barn << " barn" << G4endl;
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4PenelopeGammaConversionModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						  const G4MaterialCutsCouple* couple,
						  const G4DynamicParticle* aDynamicGamma,
						  G4double,
						  G4double)
{
  //
  // Penelope model v2008.
  // Final state is sampled according to the Bethe-Heitler model with Coulomb
  // corrections, according to the semi-empirical model of
  //  J. Baro' et al., Radiat. Phys. Chem. 44 (1994) 531.
  //
  // The model uses the high energy Coulomb correction from
  //  H. Davies et al., Phys. Rev. 93 (1954) 788
  // and atomic screening radii tabulated from
  //  J.H. Hubbel et al., J. Phys. Chem. Ref. Data 9 (1980) 1023
  // for Z= 1 to 92.
  //
  if (fVerboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopeGammaConversionModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();

  // Always kill primary
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->SetProposedKineticEnergy(0.);

  if (photonEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return ;
    }

  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();
  const G4Material* mat = couple->GetMaterial();

  //Either Initialize() was not called, or we are in a slave and InitializeLocal() was
  //not invoked
  if (!fEffectiveCharge)
    {
      //create a **thread-local** version of the table. Used only for G4EmCalculator and
      //Unit Tests
      fLocalTable = true;
      fEffectiveCharge = new std::map<const G4Material*,G4double>;
      fMaterialInvScreeningRadius = new std::map<const G4Material*,G4double>;
      fScreeningFunction = new std::map<const G4Material*,std::pair<G4double,G4double> >;
    }

  if (!fEffectiveCharge->count(mat))
    {
      //If we are here, it means that Initialize() was inkoved, but the MaterialTable was
      //not filled up. This can happen in a UnitTest or via G4EmCalculator
      if (fVerboseLevel > 0)
	{
	  //Issue a G4Exception (warning) only in verbose mode
	  G4ExceptionDescription ed;
	  ed << "Unable to allocate the EffectiveCharge data for " <<
	    mat->GetName() << G4endl;
	  ed << "This can happen only in Unit Tests" << G4endl;
	  G4Exception("G4PenelopeGammaConversionModel::SampleSecondaries()",
		      "em2019",JustWarning,ed);
	}
      //protect file reading via autolock
      G4AutoLock lock(&PenelopeGammaConversionModelMutex);
      InitializeScreeningFunctions(mat);
      lock.unlock();
    }

  // eps is the fraction of the photon energy assigned to e- (including rest mass)
  G4double eps = 0;
  G4double eki = electron_mass_c2/photonEnergy;

  //Do it fast for photon energy < 1.1 MeV (close to threshold)
  if (photonEnergy < fSmallEnergy)
    eps = eki + (1.0-2.0*eki)*G4UniformRand();
  else
    {
      //Complete calculation
      G4double effC = fEffectiveCharge->find(mat)->second;
      G4double alz = effC*fine_structure_const;
      G4double T = std::sqrt(2.0*eki);
      G4double F00=(-1.774-1.210e1*alz+1.118e1*alz*alz)*T
         +(8.523+7.326e1*alz-4.441e1*alz*alz)*T*T
         -(1.352e1+1.211e2*alz-9.641e1*alz*alz)*T*T*T
	+(8.946+6.205e1*alz-6.341e1*alz*alz)*T*T*T*T;

      G4double F0b = fScreeningFunction->find(mat)->second.second;
      G4double g0 = F0b + F00;
      G4double invRad = fMaterialInvScreeningRadius->find(mat)->second;
      G4double bmin = 4.0*eki/invRad;
      std::pair<G4double,G4double> scree =  GetScreeningFunctions(bmin);
      G4double g1 = scree.first;
      G4double g2 = scree.second;
      G4double g1min = g1+g0;
      G4double g2min = g2+g0;
      G4double xr = 0.5-eki;
      G4double a1 = 2.*g1min*xr*xr/3.;
      G4double p1 = a1/(a1+g2min);

      G4bool loopAgain = false;
      //Random sampling of eps
      do{
	loopAgain = false;
	if (G4UniformRand() <= p1)
	  {
	    G4double  ru2m1 = 2.0*G4UniformRand()-1.0;
	    if (ru2m1 < 0)
	      eps = 0.5-xr*std::pow(std::abs(ru2m1),1./3.);
	    else
	      eps = 0.5+xr*std::pow(ru2m1,1./3.);
	    G4double B = eki/(invRad*eps*(1.0-eps));
	    scree =  GetScreeningFunctions(B);
	    g1 = scree.first;
	    g1 = std::max(g1+g0,0.);
	    if (G4UniformRand()*g1min > g1)
	      loopAgain = true;
	  }
	else
	  {
	    eps = eki+2.0*xr*G4UniformRand();
	    G4double B = eki/(invRad*eps*(1.0-eps));
	    scree =  GetScreeningFunctions(B);
	    g2 = scree.second;
	    g2 = std::max(g2+g0,0.);
	    if (G4UniformRand()*g2min > g2)
	      loopAgain = true;
	  }
      }while(loopAgain);
    }
  if (fVerboseLevel > 4)
    G4cout << "Sampled eps = " << eps << G4endl;

  G4double electronTotEnergy = eps*photonEnergy;
  G4double positronTotEnergy = (1.0-eps)*photonEnergy;

  // Scattered electron (positron) angles. ( Z - axis along the parent photon)

  //electron kinematics
  G4double electronKineEnergy = std::max(0.,electronTotEnergy - electron_mass_c2) ;
  G4double costheta_el = G4UniformRand()*2.0-1.0;
  G4double kk = std::sqrt(electronKineEnergy*(electronKineEnergy+2.*electron_mass_c2));
  costheta_el = (costheta_el*electronTotEnergy+kk)/(electronTotEnergy+costheta_el*kk);
  G4double phi_el  = twopi * G4UniformRand() ;
  G4double dirX_el = std::sqrt(1.-costheta_el*costheta_el) * std::cos(phi_el);
  G4double dirY_el = std::sqrt(1.-costheta_el*costheta_el) * std::sin(phi_el);
  G4double dirZ_el = costheta_el;

  //positron kinematics
  G4double positronKineEnergy = std::max(0.,positronTotEnergy - electron_mass_c2) ;
  G4double costheta_po = G4UniformRand()*2.0-1.0;
  kk = std::sqrt(positronKineEnergy*(positronKineEnergy+2.*electron_mass_c2));
  costheta_po = (costheta_po*positronTotEnergy+kk)/(positronTotEnergy+costheta_po*kk);
  G4double phi_po  = twopi * G4UniformRand() ;
  G4double dirX_po = std::sqrt(1.-costheta_po*costheta_po) * std::cos(phi_po);
  G4double dirY_po = std::sqrt(1.-costheta_po*costheta_po) * std::sin(phi_po);
  G4double dirZ_po = costheta_po;

  // Kinematics of the created pair:
  // the electron and positron are assumed to have a symetric angular
  // distribution with respect to the Z axis along the parent photon
  G4double localEnergyDeposit = 0. ;

  if (electronKineEnergy > 0.0)
    {
      G4ThreeVector electronDirection ( dirX_el, dirY_el, dirZ_el);
      electronDirection.rotateUz(photonDirection);
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   electronDirection,
							   electronKineEnergy);
      fvect->push_back(electron);
    }
  else
    {
      localEnergyDeposit += electronKineEnergy;
      electronKineEnergy = 0;
    }

  //Generate the positron. Real particle in any case, because it will annihilate. If below
  //threshold, produce it at rest
  if (positronKineEnergy < 0.0)
    {
      localEnergyDeposit += positronKineEnergy;
      positronKineEnergy = 0; //produce it at rest
    }
  G4ThreeVector positronDirection(dirX_po,dirY_po,dirZ_po);
  positronDirection.rotateUz(photonDirection);
  G4DynamicParticle* positron = new G4DynamicParticle(G4Positron::Positron(),
						      positronDirection, positronKineEnergy);
  fvect->push_back(positron);

  //Add rest of energy to the local energy deposit
  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (fVerboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopeGammaConversion" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      if (electronKineEnergy)
	G4cout << "Electron (explicitly produced) " << electronKineEnergy/keV << " keV"
	       << G4endl;
      if (positronKineEnergy)
	G4cout << "Positron (not at rest) " << positronKineEnergy/keV << " keV" << G4endl;
      G4cout << "Rest masses of e+/- " << 2.0*electron_mass_c2/keV << " keV" << G4endl;
      if (localEnergyDeposit)
	G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (electronKineEnergy+positronKineEnergy+
					  localEnergyDeposit+2.0*electron_mass_c2)/keV <<
        " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
 if (fVerboseLevel > 0)
    {
      G4double energyDiff = std::fabs(electronKineEnergy+positronKineEnergy+
				      localEnergyDeposit+2.0*electron_mass_c2-photonEnergy);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopeGammaConversion: problem with energy conservation: "
	       << (electronKineEnergy+positronKineEnergy+
		   localEnergyDeposit+2.0*electron_mass_c2)/keV
	       << " keV (final) vs. " << photonEnergy/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeGammaConversionModel::ReadDataFile(const G4int Z)
{
  if (!IsMaster())
      //Should not be here!
    G4Exception("G4PenelopeGammaConversionModel::ReadDataFile()",
		"em0100",FatalException,"Worker thread in this method");

  if (fVerboseLevel > 2)
    {
      G4cout << "G4PenelopeGammaConversionModel::ReadDataFile()" << G4endl;
      G4cout << "Going to read Gamma Conversion data files for Z=" << Z << G4endl;
    }

    const char* path = G4FindDataDir("G4LEDATA");
    if(!path)
    {
      G4String excep =
	"G4PenelopeGammaConversionModel - G4LEDATA environment variable not set!";
      G4Exception("G4PenelopeGammaConversionModel::ReadDataFile()",
		  "em0006",FatalException,excep);
      return;
    }

  /*
    Read the cross section file
  */
  std::ostringstream ost;
  if (Z>9)
    ost << path << "/penelope/pairproduction/pdgpp" << Z << ".p08";
  else
    ost << path << "/penelope/pairproduction/pdgpp0" << Z << ".p08";
  std::ifstream file(ost.str().c_str());
  if (!file.is_open())
    {
      G4String excep = "G4PenelopeGammaConversionModel - data file " +
	G4String(ost.str()) + " not found!";
      G4Exception("G4PenelopeGammaConversionModel::ReadDataFile()",
		  "em0003",FatalException,excep);
    }

  //I have to know in advance how many points are in the data list
  //to initialize the G4PhysicsFreeVector()
  std::size_t ndata=0;
  G4String line;
  while( getline(file, line) )
    ndata++;
  ndata -= 1; //remove one header line

  file.clear();
  file.close();
  file.open(ost.str().c_str());
  G4int readZ =0;
  file >> readZ;

  if (fVerboseLevel > 3)
    G4cout << "Element Z=" << Z << G4endl;

  //check the right file is opened.
  if (readZ != Z)
    {
      G4ExceptionDescription ed;
      ed << "Corrupted data file for Z=" << Z << G4endl;
      G4Exception("G4PenelopeGammaConversionModel::ReadDataFile()",
		  "em0005",FatalException,ed);
    }

  fLogAtomicCrossSection[Z] = new G4PhysicsFreeVector(ndata);
  G4double ene=0,xs=0;
  for (std::size_t i=0;i<ndata;++i)
    {
      file >> ene >> xs;
      //dimensional quantities
      ene *= eV;
      xs *= barn;
      if (xs < 1e-40*cm2) //protection against log(0)
	xs = 1e-40*cm2;
      fLogAtomicCrossSection[Z]->PutValue(i,G4Log(ene),G4Log(xs));
    }
  file.close();

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeGammaConversionModel::InitializeScreeningFunctions(const G4Material* material)
{
  // This is subroutine GPPa0 of Penelope
  //
  // 1) calculate the effective Z for the purpose
  //
  G4double zeff = 0;
  G4int intZ = 0;
  G4int nElements = (G4int)material->GetNumberOfElements();
  const G4ElementVector* elementVector = material->GetElementVector();

  //avoid calculations if only one building element!
  if (nElements == 1)
    {
      zeff = (*elementVector)[0]->GetZ();
      intZ = (G4int) zeff;
    }
  else // many elements...let's do the calculation
    {
      const G4double* fractionVector = material->GetVecNbOfAtomsPerVolume();

      G4double atot = 0;
      for (G4int i=0;i<nElements;i++)
	{
	  G4double Zelement = (*elementVector)[i]->GetZ();
	  G4double Aelement = (*elementVector)[i]->GetAtomicMassAmu();
	  atot += Aelement*fractionVector[i];
	  zeff += Zelement*Aelement*fractionVector[i]; //average with the number of nuclei
	}
      atot /= material->GetTotNbOfAtomsPerVolume();
      zeff /= (material->GetTotNbOfAtomsPerVolume()*atot);

      intZ = (G4int) (zeff+0.25);
      if (intZ <= 0)
	intZ = 1;
      if (intZ > fMaxZ)
	intZ = fMaxZ;
    }

  if (fEffectiveCharge)
    fEffectiveCharge->insert(std::make_pair(material,zeff));

  //
  // 2) Calculate Coulomb Correction
  //
  G4double alz = fine_structure_const*zeff;
  G4double alzSquared = alz*alz;
  G4double fc =  alzSquared*(0.202059-alzSquared*
			     (0.03693-alzSquared*
			      (0.00835-alzSquared*(0.00201-alzSquared*
						   (0.00049-alzSquared*
						    (0.00012-alzSquared*0.00003)))))
			     +1.0/(alzSquared+1.0));
  //
  // 3) Screening functions and low-energy corrections
  //
  G4double matRadius = 2.0/ fAtomicScreeningRadius[intZ];
  if (fMaterialInvScreeningRadius)
    fMaterialInvScreeningRadius->insert(std::make_pair(material,matRadius));

  std::pair<G4double,G4double> myPair(0,0);
  G4double f0a = 4.0*G4Log(fAtomicScreeningRadius[intZ]);
  G4double f0b = f0a - 4.0*fc;
  myPair.first = f0a;
  myPair.second = f0b;

  if (fScreeningFunction)
    fScreeningFunction->insert(std::make_pair(material,myPair));

  if (fVerboseLevel > 2)
    {
      G4cout << "Average Z for material " << material->GetName() << " = " <<
	zeff << G4endl;
      G4cout << "Effective radius for material " << material->GetName() << " = " <<
	fAtomicScreeningRadius[intZ] << " m_e*c/hbar --> BCB = " <<
	matRadius << G4endl;
      G4cout << "Screening parameters F0 for material " << material->GetName() << " = " <<
	f0a << "," << f0b << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::pair<G4double,G4double>
G4PenelopeGammaConversionModel::GetScreeningFunctions(G4double B)
{
  // This is subroutine SCHIFF of Penelope
  //
  // Screening Functions F1(B) and F2(B) in the Bethe-Heitler differential cross
  // section for pair production
  //
  std::pair<G4double,G4double> result(0.,0.);
  G4double BSquared = B*B;
  G4double f1 = 2.0-2.0*G4Log(1.0+BSquared);
  G4double f2 = f1 - 6.66666666e-1; // (-2/3)
  if (B < 1.0e-10)
    f1 = f1-twopi*B;
  else
    {
      G4double a0 = 4.0*B*std::atan(1./B);
      f1 = f1 - a0;
      f2 += 2.0*BSquared*(4.0-a0-3.0*G4Log((1.0+BSquared)/BSquared));
    }
  G4double g1 = 0.5*(3.0*f1-f2);
  G4double g2 = 0.25*(3.0*f1+f2);

  result.first = g1;
  result.second = g2;

  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void G4PenelopeGammaConversionModel::SetParticle(const G4ParticleDefinition* p)
{
  if(!fParticle) {
    fParticle = p;
  }
}
