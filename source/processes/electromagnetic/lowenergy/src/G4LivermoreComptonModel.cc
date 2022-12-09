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
//
// Author: Alexander Bagulya
//         11 March 2012
//         on the base of G4LivermoreComptonModel
//
// History:
// --------
// 18 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - added protection against numerical problem in energy sampling 
//                  - use G4ElementSelector
// 26 Dec 2010   V Ivanchenko Load data tables only once to avoid memory leak

#include "G4LivermoreComptonModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4AutoLock.hh"
#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace { G4Mutex LivermoreComptonModelMutex = G4MUTEX_INITIALIZER; }
using namespace std;

G4PhysicsFreeVector* G4LivermoreComptonModel::data[] = {nullptr};
G4ShellData*       G4LivermoreComptonModel::shellData = nullptr;
G4DopplerProfile*  G4LivermoreComptonModel::profileData = nullptr;

static const G4double ln10 = G4Log(10.);

G4LivermoreComptonModel::G4LivermoreComptonModel(const G4ParticleDefinition*, 
						 const G4String& nam)
  : G4VEmModel(nam),maxZ(100),isInitialised(false)
{  
  verboseLevel=1 ;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(  verboseLevel>1 ) { 
    G4cout << "Livermore Compton model is constructed " << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  fParticleChange = 0;
  fAtomDeexcitation = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreComptonModel::~G4LivermoreComptonModel()
{
  if(IsMaster())
  {
    delete shellData;
    shellData = nullptr;
    delete profileData;
    profileData = nullptr;
    for(G4int i = 0; i <= maxZ; ++i)
    {
      if(data[i])
      {
        delete data[i];
        data[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModel::Initialise(const G4ParticleDefinition* particle,
					 const G4DataVector& cuts)
{
  if (verboseLevel > 1) {
    G4cout << "Calling G4LivermoreComptonModel::Initialise()" << G4endl;
  }

  // Initialise element selector  
  if(IsMaster()) {
    // Access to elements
    const char* path = G4FindDataDir("G4LEDATA");

    const G4ElementTable* elemTable = G4Element::GetElementTable();
    size_t numElems                 = (*elemTable).size();
    for(size_t ie = 0; ie < numElems; ++ie)
    {
      const G4Element* elem = (*elemTable)[ie];
      const G4int Z         = std::min(maxZ, elem->GetZasInt());

      if(data[Z] == nullptr)
      {
        ReadData(Z, path);
      }
    }

    // For Doppler broadening
    if(shellData == nullptr) {
      shellData = new G4ShellData(); 
      shellData->SetOccupancyData();
      G4String file = "/doppler/shell-doppler";
      shellData->LoadData(file);
    }
    if(profileData == nullptr) { profileData = new G4DopplerProfile(); }

    InitialiseElementSelectors(particle, cuts);
  }

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files" << G4endl;
  }
 
  if( verboseLevel>1 ) { 
    G4cout << "G4LivermoreComptonModel is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
  //  
  if(isInitialised) { return; }

  fParticleChange = GetParticleChangeForGamma();
  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModel::InitialiseLocal(const G4ParticleDefinition*,
					      G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreComptonModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
  {
    G4cout << "G4LivermoreComptonModel::ReadData()" 
	   << G4endl;
  }
  if(data[Z] != nullptr) { return; }
  const char* datadir = path;
  if(datadir == nullptr)
  {
    datadir = G4FindDataDir("G4LEDATA");
    if(datadir == nullptr)
    {
      G4Exception("G4LivermoreComptonModel::ReadData()",
		  "em0006",FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }
  
  data[Z] = new G4PhysicsFreeVector();
    
  std::ostringstream ost;
  if(G4EmParameters::Instance()->LivermoreDataDir() == "livermore"){
    ost << datadir << "/livermore/comp/ce-cs-" << Z <<".dat";
  }else{
    ost << datadir << "/epics2017/comp/ce-cs-" << Z <<".dat";
  }

  std::ifstream fin(ost.str().c_str());
  
  if( !fin.is_open()) 
    {
      G4ExceptionDescription ed;
      ed << "G4LivermoreComptonModel data file <" << ost.str().c_str()
	 << "> is not opened!" << G4endl;
    G4Exception("G4LivermoreComptonModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW8.0 or later");
      return;
    } else {
      if(verboseLevel > 3) {
	G4cout << "File " << ost.str() 
	       << " is opened by G4LivermoreComptonModel" << G4endl;
      }   
      data[Z]->Retrieve(fin, true);
      data[Z]->ScaleVector(MeV, MeV*barn);
    }   
  fin.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreComptonModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						    G4double GammaEnergy,
						    G4double Z, G4double,
						    G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "G4LivermoreComptonModel::ComputeCrossSectionPerAtom()" 
	   << G4endl;
  }
  G4double cs = 0.0; 

  if (GammaEnergy < LowEnergyLimit()) { return 0.0; }

  G4int intZ = G4lrint(Z);
  if(intZ < 1 || intZ > maxZ) { return cs; } 

  G4PhysicsFreeVector* pv = data[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(pv == nullptr)
    {
      InitialiseForElement(0, intZ);
      pv = data[intZ];
      if(pv == nullptr) { return cs; }
    }

  G4int n = G4int(pv->GetVectorLength() - 1);   
  G4double e1 = pv->Energy(0);
  G4double e2 = pv->Energy(n);

  if(GammaEnergy <= e1)      { cs = GammaEnergy/(e1*e1)*pv->Value(e1); }
  else if(GammaEnergy <= e2) { cs = pv->Value(GammaEnergy)/GammaEnergy; }
  else if(GammaEnergy > e2)  { cs = pv->Value(e2)/GammaEnergy; }
 
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void G4LivermoreComptonModel::SampleSecondaries(
                             std::vector<G4DynamicParticle*>* fvect,
			     const G4MaterialCutsCouple* couple,
			     const G4DynamicParticle* aDynamicGamma,
			     G4double, G4double)
{

  // The scattered gamma energy is sampled according to Klein - Nishina 
  // formula then accepted or rejected depending on the Scattering Function 
  // multiplied by factor from Klein - Nishina formula.
  // Expression of the angular distribution as Klein Nishina
  // angular and energy distribution and Scattering fuctions is taken from
  // D. E. Cullen "A simple model of photon transport" Nucl. Instr. Meth.
  // Phys. Res. B 101 (1995). Method of sampling with form factors is different
  // data are interpolated while in the article they are fitted.
  // Reference to the article is from J. Stepanek New Photon, Positron
  // and Electron Interaction Data for GEANT in Energy Range from 1 eV to 10
  // TeV (draft).
  // The random number techniques of Butcher & Messel are used
  // (Nucl Phys 20(1960),15).

  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  if (verboseLevel > 3) {
    G4cout << "G4LivermoreComptonModel::SampleSecondaries() E(MeV)= " 
	   << photonEnergy0/MeV << " in " << couple->GetMaterial()->GetName() 
	   << G4endl;
  }

  // do nothing below the threshold
  // should never get here because the XS is zero below the limit
  if (photonEnergy0 < LowEnergyLimit())     
    return ;    

  G4double e0m = photonEnergy0 / electron_mass_c2 ;
  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple,particle,photonEnergy0);

  G4int Z = elm->GetZasInt();

  G4double epsilon0Local = 1. / (1. + 2. * e0m);
  G4double epsilon0Sq = epsilon0Local * epsilon0Local;
  G4double alpha1 = -G4Log(epsilon0Local);
  G4double alpha2 = 0.5 * (1. - epsilon0Sq);

  G4double wlPhoton = h_Planck*c_light/photonEnergy0;

  // Sample the energy of the scattered photon
  G4double epsilon;
  G4double epsilonSq;
  G4double oneCosT;
  G4double sinT2;
  G4double gReject;

  if (verboseLevel > 3) {
    G4cout << "Started loop to sample gamma energy" << G4endl;
  }
  
  do {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand())
      {
        epsilon = G4Exp(-alpha1 * G4UniformRand());  
        epsilonSq = epsilon * epsilon;
      }
      else
      {
        epsilonSq = epsilon0Sq + (1. - epsilon0Sq) * G4UniformRand();
        epsilon = std::sqrt(epsilonSq);
      }

      oneCosT = (1. - epsilon) / ( epsilon * e0m);
      sinT2 = oneCosT * (2. - oneCosT);
      G4double x = std::sqrt(oneCosT/2.) * cm / wlPhoton;
      G4double scatteringFunction = ComputeScatteringFunction(x, Z);
      gReject = (1. - epsilon * sinT2 / (1. + epsilonSq)) * scatteringFunction;

  } while(gReject < G4UniformRand()*Z);

  G4double cosTheta = 1. - oneCosT;
  G4double sinTheta = std::sqrt (sinT2);
  G4double phi = twopi * G4UniformRand() ;
  G4double dirx = sinTheta * std::cos(phi);
  G4double diry = sinTheta * std::sin(phi);
  G4double dirz = cosTheta ;

  // Doppler broadening -  Method based on:
  // Y. Namito, S. Ban and H. Hirayama, 
  // "Implementation of the Doppler Broadening of a Compton-Scattered Photon 
  // into the EGS4 Code", NIM A 349, pp. 489-494, 1994
  
  // Maximum number of sampling iterations
  static G4int maxDopplerIterations = 1000;
  G4double bindingE = 0.;
  G4double photonEoriginal = epsilon * photonEnergy0;
  G4double photonE = -1.;
  G4int iteration = 0;
  G4double eMax = photonEnergy0;

  G4int shellIdx = 0;

  if (verboseLevel > 3) {
    G4cout << "Started loop to sample broading" << G4endl;
  }

  do
    {
      ++iteration;
      // Select shell based on shell occupancy
      shellIdx = shellData->SelectRandomShell(Z);
      bindingE = shellData->BindingEnergy(Z,shellIdx);

      if (verboseLevel > 3) {
	G4cout << "Shell ID= " << shellIdx 
	       << " Ebind(keV)= " << bindingE/keV << G4endl;
      }
      
      eMax = photonEnergy0 - bindingE;
     
      // Randomly sample bound electron momentum 
      // (memento: the data set is in Atomic Units)
      G4double pSample = profileData->RandomSelectMomentum(Z,shellIdx);
      if (verboseLevel > 3) {
	G4cout << "pSample= " << pSample << G4endl;
      }
      // Rescale from atomic units
      G4double pDoppler = pSample * fine_structure_const;
      G4double pDoppler2 = pDoppler * pDoppler;
      G4double var2 = 1. + oneCosT * e0m;
      G4double var3 = var2*var2 - pDoppler2;
      G4double var4 = var2 - pDoppler2 * cosTheta;
      G4double var = var4*var4 - var3 + pDoppler2 * var3;
      if (var > 0.)
	{
	  G4double varSqrt = std::sqrt(var);        
	  G4double scale = photonEnergy0 / var3;  
          // Random select either root
 	  if (G4UniformRand() < 0.5) { photonE = (var4 - varSqrt) * scale; }
	  else                       { photonE = (var4 + varSqrt) * scale; }
	} 
      else
	{
	  photonE = -1.;
	}
   } while (iteration <= maxDopplerIterations && photonE > eMax); 
	     
  // End of recalculation of photon energy with Doppler broadening
  // Revert to original if maximum number of iterations threshold 
  // has been reached
  if (iteration >= maxDopplerIterations)
    {
      photonE = photonEoriginal;
      bindingE = 0.;
    }

  // Update G4VParticleChange for the scattered photon
  G4ThreeVector photonDirection1(dirx,diry,dirz);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1) ;

  G4double photonEnergy1 = photonE;

  if (photonEnergy1 > 0.) {
    fParticleChange->SetProposedKineticEnergy(photonEnergy1) ;

  } else {
    // photon absorbed
    photonEnergy1 = 0.;
    fParticleChange->SetProposedKineticEnergy(0.) ;
    fParticleChange->ProposeTrackStatus(fStopAndKill);   
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0);
    return;
  }

  // Kinematics of the scattered electron
  G4double eKineticEnergy = photonEnergy0 - photonEnergy1 - bindingE;

  // protection against negative final energy: no e- is created
  if(eKineticEnergy < 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy0 - photonEnergy1);
    return;
  }

  G4double eTotalEnergy = eKineticEnergy + electron_mass_c2;

  G4double electronE = photonEnergy0 * (1. - epsilon) + electron_mass_c2; 
  G4double electronP2 = 
    electronE*electronE - electron_mass_c2*electron_mass_c2;
  G4double sinThetaE = -1.;
  G4double cosThetaE = 0.;
  if (electronP2 > 0.)
    {
      cosThetaE = (eTotalEnergy + photonEnergy1 )* 
	(1. - epsilon) / std::sqrt(electronP2);
      sinThetaE = -1. * sqrt(1. - cosThetaE * cosThetaE); 
    }
  
  G4double eDirX = sinThetaE * std::cos(phi);
  G4double eDirY = sinThetaE * std::sin(phi);
  G4double eDirZ = cosThetaE;

  G4ThreeVector eDirection(eDirX,eDirY,eDirZ);
  eDirection.rotateUz(photonDirection0);
  G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),
						 eDirection,eKineticEnergy) ;
  fvect->push_back(dp);

  // sample deexcitation
  if (verboseLevel > 3) {
    G4cout << "Started atomic de-excitation " << fAtomDeexcitation << G4endl;
  }

  if(nullptr != fAtomDeexcitation && iteration < maxDopplerIterations) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t i=nbefore; i<nafter; ++i) {
          //Check if there is enough residual energy 
          if (bindingE >= ((*fvect)[i])->GetKineticEnergy())
           {
             //Ok, this is a valid secondary: keep it
	     bindingE -= ((*fvect)[i])->GetKineticEnergy();
           }
          else
           {
 	     //Invalid secondary: not enough energy to create it!
 	     //Keep its energy in the local deposit
             delete (*fvect)[i]; 
             (*fvect)[i]=0;
           }
	} 
      }
    }
  }
  bindingE = std::max(bindingE, 0.0);
  fParticleChange->ProposeLocalEnergyDeposit(bindingE);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4LivermoreComptonModel::ComputeScatteringFunction(G4double x, G4int Z)
{
  G4double value = Z;
  if(x <= ScatFuncFitParam[Z][3])
  {
    G4double lgq = G4Log(x) / ln10;

    if(lgq < ScatFuncFitParam[Z][1])
    {
      value = ScatFuncFitParam[Z][4] + lgq * ScatFuncFitParam[Z][5];
    }
    else if(lgq >= ScatFuncFitParam[Z][1] && lgq < ScatFuncFitParam[Z][2])
    {
      value = ScatFuncFitParam[Z][6] +
              lgq * (ScatFuncFitParam[Z][7] +
                     lgq * (ScatFuncFitParam[Z][8] +
                            lgq * (ScatFuncFitParam[Z][9] +
                                   lgq * ScatFuncFitParam[Z][10])));
    }
    else
    {
      value = ScatFuncFitParam[Z][11] +
              lgq * (ScatFuncFitParam[Z][12] +
                     lgq * (ScatFuncFitParam[Z][13] +
                            lgq * (ScatFuncFitParam[Z][14] +
                                   lgq * ScatFuncFitParam[Z][15])));
    }
    value = G4Exp(value * ln10);
  }
  // G4cout << "    value= " << value << G4endl;
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4LivermoreComptonModel::InitialiseForElement(const G4ParticleDefinition*, 
					      G4int Z)
{
  G4AutoLock l(&LivermoreComptonModelMutex);
  if(data[Z] == nullptr) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Fitting data to compute scattering function based on EPICS2017
const G4double G4LivermoreComptonModel::ScatFuncFitParam[101][16] = {
  {0,              0.,             0.,              0.,              0.,              0.,               0.,              0.,               0.,              0.,               0.,              0.,               0.,              0.,               0.,              0.},
  {1, 6.000000000e+00, 7.087999300e+00, 1.499680000e+08, -1.435559123e+01, 2.000000043e+00, -3.925518125e+02, 2.434944521e+02, -5.784393623e+01, 6.160181204e+00, -2.461326602e-01, -1.649463594e+03, 8.121933215e+02, -1.498313316e+02, 1.227279742e+01, -3.765996345e-01},
  {2, 6.000000000e+00, 7.199000403e+00, 2.500350000e+08, -1.430103027e+01, 2.000000041e+00, 3.574019365e+02, -1.978574937e+02, 3.971327838e+01, -3.443224867e+00, 1.091825227e-01, -4.009960832e+02, 1.575831469e+02, -2.174763446e+01, 1.185163045e+00, -1.814503741e-02},
  {3, 6.000000000e+00, 7.301000136e+00, 3.999450000e+08, -1.357675458e+01, 2.000000074e+00, 7.051635443e+02, -4.223841786e+02, 9.318729225e+01, -9.002642767e+00, 3.220625771e-01, 1.524679907e+03, -7.851479582e+02, 1.509941052e+02, -1.285477984e+01, 4.089348830e-01},
  {4, 6.000000000e+00, 7.349500202e+00, 5.000350000e+08, -1.375202671e+01, 1.999999994e+00, -1.832909604e+02, 1.193997722e+02, -3.034328318e+01, 3.471545044e+00, -1.484222463e-01, 1.397476657e+03, -7.026416933e+02, 1.320720559e+02, -1.099824430e+01, 3.424610532e-01},
  {5, 6.000000000e+00, 7.388999972e+00, 5.997910000e+08, -1.380548571e+01, 2.000000004e+00, -2.334197545e+02, 1.467013466e+02, -3.574851109e+01, 3.925047955e+00, -1.616186492e-01, 6.784713308e+02, -3.419562074e+02, 6.433945831e+01, -5.354244209e+00, 1.663784966e-01},
  {6, 6.000000000e+00, 7.422500001e+00, 6.998420000e+08, -1.388639003e+01, 1.999999863e+00, -2.460254935e+02, 1.516613633e+02, -3.622024219e+01, 3.900099543e+00, -1.576557530e-01, -1.610185428e+02, 7.010907070e+01, -1.142375397e+01, 8.303365180e-01, -2.273786010e-02},
  {7, 6.000000000e+00, 7.451499931e+00, 7.998340000e+08, -1.388605429e+01, 1.999999612e+00, -3.054540719e+02, 1.877740247e+02, -4.440273010e+01, 4.718886370e+00, -1.881615004e-01, -2.263864349e+02, 1.017885461e+02, -1.716982752e+01, 1.292954622e+00, -3.668301946e-02},
  {8, 6.000000000e+00, 7.451499931e+00, 7.998340000e+08, -1.395860675e+01, 1.999999906e+00, -3.877174895e+02, 2.345831969e+02, -5.431822300e+01, 5.643262324e+00, -2.200840540e-01, -7.949384302e+02, 3.757293602e+02, -6.661741851e+01, 5.256265086e+00, -1.556986777e-01},
  {9, 6.000000000e+00, 7.451499931e+00, 7.998340000e+08, -1.400000063e+01, 2.000000106e+00, -2.939854827e+02, 1.784214589e+02, -4.168473845e+01, 4.377669850e+00, -1.724300716e-01, -1.169326170e+03, 5.545642014e+02, -9.863024948e+01, 7.801721240e+00, -2.315522357e-01},
  {10, 6.000000000e+00, 7.451499931e+00, 7.998340000e+08, -1.404575854e+01, 2.000000178e+00, -2.615701853e+02, 1.582596311e+02, -3.698114811e+01, 3.889093901e+00, -1.533613504e-01, -1.275287356e+03, 6.022076554e+02, -1.066410301e+02, 8.398773148e+00, -2.481899800e-01},
  {11, 6.000000000e+00, 7.500000000e+00, 1.000000000e+09, -1.344369665e+01, 1.999999860e+00, 1.112662501e+03, -6.807056448e+02, 1.545837472e+02, -1.548462180e+01, 5.785425068e-01, -1.007702307e+03, 4.699937040e+02, -8.220352105e+01, 6.396099420e+00, -1.867816054e-01},
  {12, 6.000000000e+00, 7.500000000e+00, 1.000000000e+09, -1.339794047e+01, 2.000000080e+00, 9.895649717e+02, -5.983228286e+02, 1.340681576e+02, -1.323046651e+01, 4.863434994e-01, -5.790532602e+02, 2.626052403e+02, -4.463548055e+01, 3.376239891e+00, -9.588786915e-02},
  {13, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.340893585e+01, 2.000000078e+00, 7.335256091e+02, -4.405291562e+02, 9.770954287e+01, -9.519317788e+00, 3.448067237e-01, -5.328832253e+02, 2.398514938e+02, -4.044557740e+01, 3.034597500e+00, -8.547410419e-02},
  {14, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.345593195e+01, 2.000000000e+00, 3.978691889e+02, -2.370975001e+02, 5.158692183e+01, -4.884868277e+00, 1.707270518e-01, -2.340256277e+02, 9.813362251e+01, -1.527892110e+01, 1.051070768e+00, -2.692716945e-02},
  {15, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.349485049e+01, 2.000000083e+00, 2.569833671e+02, -1.513623448e+02, 3.210087153e+01, -2.925756803e+00, 9.724379436e-02, -1.345727293e+01, -6.291081167e+00, 3.235960888e+00, -4.059236666e-01, 1.601245178e-02},
  {16, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.353760159e+01, 1.999999937e+00, 1.015293074e+02, -5.721639224e+01, 1.078607152e+01, -7.890593144e-01, 1.726056327e-02, 1.854818165e+02, -1.000803879e+02, 1.979815884e+01, -1.704221744e+00, 5.413372375e-02},
  {17, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.358502705e+01, 2.000000066e+00, -4.294163461e+01, 2.862162412e+01, -8.285972104e+00, 1.087745268e+00, -5.172153610e-02, 1.676674074e+02, -8.976414784e+01, 1.763329621e+01, -1.507161653e+00, 4.753277254e-02},
  {18, 6.000000000e+00, 7.587999300e+00, 1.499680000e+09, -1.361978902e+01, 2.000000042e+00, -3.573422746e+01, 2.403066369e+01, -7.173617800e+00, 9.657608431e-01, -4.662317662e-02, 1.811925229e+02, -9.574636323e+01, 1.861940167e+01, -1.578810247e+00, 4.946799877e-02},
  {19, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.320760816e+01, 1.999999979e+00, 1.263152069e+02, -8.738932892e+01, 2.109042182e+01, -2.166733566e+00, 8.146018979e-02, 9.183312428e+01, -5.232836676e+01, 1.072450810e+01, -9.419512971e-01, 3.023884410e-02},
  {20, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.314266674e+01, 1.999999876e+00, 6.620218058e+02, -4.057504297e+02, 9.180787767e+01, -9.124184449e+00, 3.372518137e-01, 7.034138711e+01, -4.198325416e+01, 8.861351614e+00, -7.930506530e-01, 2.578454342e-02},
  {21, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.317392498e+01, 1.999999966e+00, 6.766093786e+02, -4.129087029e+02, 9.305090790e+01, -9.212128925e+00, 3.392408033e-01, 1.916559096e+01, -1.807294109e+01, 4.677205921e+00, -4.679350245e-01, 1.632115420e-02},
  {22, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.320065945e+01, 1.999999999e+00, 6.969823082e+02, -4.236620289e+02, 9.513714106e+01, -9.388294642e+00, 3.446942719e-01, -6.501317146e+01, 2.138553133e+01, -2.250998891e+00, 7.219326079e-02, 5.467529893e-04},
  {23, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.322914744e+01, 1.999999909e+00, 6.889749928e+02, -4.181421624e+02, 9.373529727e+01, -9.233142268e+00, 3.383772151e-01, -1.382770534e+02, 5.540647456e+01, -8.170017489e+00, 5.295569200e-01, -1.269556386e-02},
  {24, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.333724128e+01, 1.999999854e+00, 4.365566411e+02, -2.672774427e+02, 6.001631369e+01, -5.895458454e+00, 2.149710735e-01, -2.393534124e+02, 1.020845165e+02, -1.624744211e+01, 1.150387566e+00, -3.057723021e-02},
  {25, 6.000000000e+00, 7.650499797e+00, 1.999860000e+09, -1.328399669e+01, 2.000000008e+00, 6.461381990e+02, -3.918546518e+02, 8.769548644e+01, -8.618784385e+00, 3.150660827e-01, -2.597409979e+02, 1.113332866e+02, -1.782124571e+01, 1.269519197e+00, -3.396126698e-02},
  {26, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.330103000e+01, 1.999999998e+00, 4.261007401e+02, -2.588846763e+02, 5.764613910e+01, -5.609660122e+00, 2.024165636e-01, -1.982896712e+02, 8.274273985e+01, -1.284074215e+01, 8.845687432e-01, -2.282143299e-02},
  {27, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.332790165e+01, 1.999999922e+00, 4.006816638e+02, -2.439311564e+02, 5.435031497e+01, -5.287693457e+00, 1.906696163e-01, -2.205075564e+02, 9.262919772e+01, -1.448909443e+01, 1.006686819e+00, -2.621294059e-02},
  {28, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.334678710e+01, 1.999999939e+00, 3.967750019e+02, -2.411866801e+02, 5.364872608e+01, -5.210295834e+00, 1.875525119e-01, -2.516823030e+02, 1.065117131e+02, -1.680533335e+01, 1.178363534e+00, -3.098194406e-02},
  {29, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.344369664e+01, 1.999999853e+00, 2.437671888e+02, -1.499592208e+02, 3.332221026e+01, -3.206587185e+00, 1.138639692e-01, -2.874130637e+02, 1.223381969e+02, -1.943178054e+01, 1.371979484e+00, -3.633119448e-02},
  {30, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.338721562e+01, 1.999999911e+00, 3.914867984e+02, -2.378147085e+02, 5.284517777e+01, -5.126420186e+00, 1.843322562e-01, -3.235063319e+02, 1.384252948e+02, -2.211844479e+01, 1.571300198e+00, -4.187323186e-02},
  {31, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.335654643e+01, 1.999999847e+00, 4.325820127e+02, -2.614587597e+02, 5.793273998e+01, -5.611190206e+00, 2.015836827e-01, -3.359152840e+02, 1.437507638e+02, -2.297457475e+01, 1.632470701e+00, -4.351215346e-02},
  {32, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.337675047e+01, 1.999999960e+00, 4.388195965e+02, -2.642662297e+02, 5.834159168e+01, -5.629419790e+00, 2.014339673e-01, -3.430730654e+02, 1.467102631e+02, -2.343160019e+01, 1.663765504e+00, -4.431369286e-02},
  {33, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.339794046e+01, 2.000000074e+00, 3.931399547e+02, -2.363700718e+02, 5.197696913e+01, -4.987097655e+00, 1.772567576e-01, -3.501570134e+02, 1.497141578e+02, -2.390888062e+01, 1.697503580e+00, -4.520887478e-02},
  {34, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.342021680e+01, 2.000000064e+00, 3.772588127e+02, -2.256347960e+02, 4.929790851e+01, -4.694628847e+00, 1.654667382e-01, -3.481053019e+02, 1.486490112e+02, -2.370745096e+01, 1.680991482e+00, -4.471064364e-02},
  {35, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.344369666e+01, 1.999999864e+00, 3.344685842e+02, -1.994816236e+02, 4.332267376e+01, -4.090542180e+00, 1.426839031e-01, -3.227660675e+02, 1.370301996e+02, -2.171543883e+01, 1.529681552e+00, -4.041331983e-02},
  {36, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.345593194e+01, 1.999999999e+00, 3.004054446e+02, -1.781334135e+02, 3.834850324e+01, -3.580074471e+00, 1.232168921e-01, -2.980827664e+02, 1.257508661e+02, -1.978792154e+01, 1.383723149e+00, -3.628014907e-02},
  {37, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.310790583e+01, 2.000000075e+00, -3.687188343e+01, 1.054409719e+01, -8.516586814e-01, 9.339751003e-03, 8.809383936e-04, -2.699384784e+02, 1.129635316e+02, -1.761447452e+01, 1.219971043e+00, -3.166503704e-02},
  {38, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.304095795e+01, 1.999999892e+00, 1.969969064e+02, -1.286503864e+02, 3.008431767e+01, -3.031946980e+00, 1.124456346e-01, -2.331258613e+02, 9.627987243e+01, -1.478515961e+01, 1.007215642e+00, -2.567873120e-02},
  {39, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.306048023e+01, 1.999999916e+00, 2.891710763e+02, -1.819536752e+02, 4.158265841e+01, -4.128940218e+00, 1.515168697e-01, -1.997404800e+02, 8.119476676e+01, -1.223426670e+01, 8.159269666e-01, -2.031079820e-02},
  {40, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.308092198e+01, 2.000000013e+00, 3.393782172e+02, -2.103908454e+02, 4.758278737e+01, -4.688308235e+00, 1.709723418e-01, -1.549247582e+02, 6.091403935e+01, -8.799307373e+00, 5.578963961e-01, -1.305663921e-02},
  {41, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.316749062e+01, 1.999999920e+00, 2.748604341e+02, -1.706429616e+02, 3.843757441e+01, -3.759045290e+00, 1.358263430e-01, -1.163607425e+02, 4.350905533e+01, -5.859305970e+00, 3.376426246e-01, -6.881281652e-03},
  {42, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.318708720e+01, 2.000000093e+00, 3.203285955e+02, -1.966282865e+02, 4.398204769e+01, -4.283031482e+00, 1.543480828e-01, -9.364181222e+01, 3.329814493e+01, -4.141689265e+00, 2.095170962e-01, -3.304665813e-03},
  {43, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.314266674e+01, 1.999999876e+00, 4.184977165e+02, -2.552902161e+02, 5.707764818e+01, -5.576436872e+00, 2.020184726e-01, -8.395646154e+01, 2.898228589e+01, -3.422356654e+00, 1.564059753e-01, -1.838508896e-03},
  {44, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.322914744e+01, 1.999999909e+00, 3.243555305e+02, -1.978255470e+02, 4.397580841e+01, -4.256142657e+00, 1.524431452e-01, -5.506292375e+01, 1.599310639e+01, -1.237152904e+00, -6.611574411e-03, 2.712232383e-03},
  {45, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.325181249e+01, 2.000000089e+00, 3.037823599e+02, -1.856628295e+02, 4.128167884e+01, -3.991656133e+00, 1.427469878e-01, -5.014186072e+01, 1.386962969e+01, -8.950806420e-01, -3.095321225e-02, 3.357984426e-03},
  {46, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.340893584e+01, 2.000000073e+00, 3.529797051e+02, -2.101512262e+02, 4.563946029e+01, -4.315279704e+00, 1.509248358e-01, -4.815922691e+01, 1.301508788e+01, -7.580854951e-01, -4.059091985e-02, 3.608993811e-03},
  {47, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.328399669e+01, 2.000000008e+00, 3.074953924e+02, -1.872462583e+02, 4.149827252e+01, -4.000811852e+00, 1.426973118e-01, -4.897188379e+01, 1.335300002e+01, -8.110051997e-01, -3.684788190e-02, 3.508156457e-03},
  {48, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.322914743e+01, 1.999999904e+00, 4.059717166e+02, -2.462737702e+02, 5.472040126e+01, -5.311320062e+00, 1.911670149e-01, -5.901534554e+01, 1.791385249e+01, -1.587065943e+00, 2.182673278e-02, 1.845559896e-03},
  {49, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.320760815e+01, 1.999999973e+00, 4.369774251e+02, -2.639721687e+02, 5.849617557e+01, -5.667842049e+00, 2.037342202e-01, -7.399698219e+01, 2.469785523e+01, -2.737881327e+00, 1.085351830e-01, -6.022720695e-04},
  {50, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.322184869e+01, 1.999999993e+00, 4.289361021e+02, -2.585593024e+02, 5.714058683e+01, -5.518600115e+00, 1.976499817e-01, -9.269047286e+01, 3.314422349e+01, -4.167341855e+00, 2.159629039e-01, -3.626802503e-03},
  {51, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.323657166e+01, 1.999999946e+00, 3.866985836e+02, -2.328379698e+02, 5.128884878e+01, -4.929614910e+00, 1.755331333e-01, -1.067869310e+02, 3.950715983e+01, -5.243321447e+00, 2.967791238e-01, -5.901223876e-03},
  {52, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.325181248e+01, 2.000000083e+00, 3.947511198e+02, -2.363799049e+02, 5.179393756e+01, -4.951603918e+00, 1.753404387e-01, -1.069681982e+02, 3.995521754e+01, -5.382071424e+00, 3.120248901e-01, -6.467957474e-03},
  {53, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.326760745e+01, 2.000000205e+00, 3.694394448e+02, -2.204699428e+02, 4.806381052e+01, -4.565474883e+00, 1.604614344e-01, -1.180749905e+02, 4.460080701e+01, -6.105217447e+00, 3.616537171e-01, -7.733059623e-03},
  {54, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.328399667e+01, 2.000000001e+00, 3.423943987e+02, -2.041330669e+02, 4.437639784e+01, -4.197363553e+00, 1.467594367e-01, -1.288973984e+02, 4.985324046e+01, -7.056041375e+00, 4.378018318e-01, -1.000965926e-02},
  {55, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.297881025e+01, 1.999999927e+00, -7.663422017e+01, 3.462700567e+01, -6.273553579e+00, 5.487612834e-01, -1.912897528e-02, -1.318428276e+02, 5.081036112e+01, -7.154907590e+00, 4.405355674e-01, -9.955685075e-03},
  {56, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.290657751e+01, 1.999999869e+00, 1.084179205e+02, -7.602229206e+01, 1.843754298e+01, -1.892451591e+00, 7.085434176e-02, -1.346311376e+02, 5.207427468e+01, -7.369834199e+00, 4.568138610e-01, -1.041859875e-02},
  {57, 6.000000000e+00, 7.725500002e+00, 2.824880000e+09, -1.292445241e+01, 1.999999898e+00, 2.995898890e+02, -1.889477671e+02, 4.336642429e+01, -4.330424108e+00, 1.599942758e-01, 5.503972208e+00, -1.227641064e+01, 3.699182312e+00, -3.884476060e-01, 1.375966896e-02},
  {58, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.293554133e+01, 1.999999890e+00, 1.709135500e+02, -1.120124681e+02, 2.615893820e+01, -2.624416758e+00, 9.674223967e-02, -1.375860132e+02, 5.337811974e+01, -7.586786386e+00, 4.730023198e-01, -1.087482303e-02},
  {59, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.293554133e+01, 1.999999890e+00, 1.214691988e+02, -8.336119630e+01, 1.996468944e+01, -2.032283439e+00, 7.562254632e-02, -1.631005912e+02, 6.472051894e+01, -9.476098737e+00, 6.127875286e-01, -1.475060958e-02},
  {60, 6.000000000e+00, 7.849500202e+00, 5.000350000e+09, -1.294309494e+01, 1.999999967e+00, 1.302719596e+02, -8.835087414e+01, 2.101971144e+01, -2.131084478e+00, 7.908549730e-02, -1.692901279e+02, 6.742727614e+01, -9.920661139e+00, 6.453186854e-01, -1.564524492e-02},
  {61, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.295078139e+01, 1.999999905e+00, 1.127680235e+02, -7.782238836e+01, 1.865126163e+01, -1.895116816e+00, 7.030502833e-02, -2.059821608e+02, 8.384774285e+01, -1.267344799e+01, 8.502354115e-01, -2.135994609e-02},
  {62, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.295860692e+01, 1.999999936e+00, 1.203145109e+02, -8.212556537e+01, 1.956606386e+01, -1.981212240e+00, 7.333626288e-02, -2.158058793e+02, 8.810144391e+01, -1.336380022e+01, 9.000362964e-01, -2.270715579e-02},
  {63, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.296657573e+01, 1.999999918e+00, 1.212159597e+02, -8.256559477e+01, 1.964122173e+01, -1.986442056e+00, 7.345564343e-02, -2.278531434e+02, 9.336519465e+01, -1.422588608e+01, 9.627883381e-01, -2.441986614e-02},
  {64, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.298296617e+01, 1.999999921e+00, 1.689382403e+02, -1.099987696e+02, 2.551961464e+01, -2.543234152e+00, 9.313568005e-02, -2.282716670e+02, 9.348611199e+01, -1.423588448e+01, 9.628551072e-01, -2.440492772e-02},
  {65, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.299139910e+01, 1.999999880e+00, 1.724155378e+02, -1.120798437e+02, 2.598264738e+01, -2.588807295e+00, 9.481417896e-02, -2.322687147e+02, 9.517466656e+01, -1.450332749e+01, 9.817069914e-01, -2.490386807e-02},
  {66, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.298716240e+01, 1.999999941e+00, 1.286079419e+02, -8.646296410e+01, 2.039801258e+01, -2.050839207e+00, 7.549033493e-02, -2.420048480e+02, 9.935663043e+01, -1.517653800e+01, 1.029875015e+00, -2.619626869e-02},
  {67, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.299567846e+01, 1.999999971e+00, 1.182799697e+02, -8.043389241e+01, 1.908027783e+01, -1.923209794e+00, 7.087268462e-02, -2.464462609e+02, 1.012059056e+02, -1.546468270e+01, 1.049814070e+00, -2.671320158e-02},
  {68, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.300436459e+01, 1.999999966e+00, 1.150510247e+02, -7.859576077e+01, 1.868688175e+01, -1.885844183e+00, 6.954765052e-02, -2.457555063e+02, 1.007538481e+02, -1.536692833e+01, 1.041070997e+00, -2.643279207e-02},
  {69, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.300877391e+01, 2.000000000e+00, 1.266280406e+02, -8.514491730e+01, 2.007089332e+01, -2.015475088e+00, 7.409191965e-02, -2.492442707e+02, 1.021615320e+02, -1.557878384e+01, 1.055183253e+00, -2.678362279e-02},
  {70, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.301772826e+01, 1.999999912e+00, 1.224253568e+02, -8.281395858e+01, 1.958609738e+01, -1.970785167e+00, 7.255458061e-02, -2.488808342e+02, 1.018569466e+02, -1.550601866e+01, 1.048325396e+00, -2.655661748e-02},
  {71, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.303151733e+01, 2.000000051e+00, 1.862181262e+02, -1.199038630e+02, 2.763107534e+01, -2.742586837e+00, 1.001956495e-01, -2.403102476e+02, 9.796272016e+01, -1.484525920e+01, 9.987147871e-01, -2.516533876e-02},
  {72, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.304575796e+01, 2.000000081e+00, 2.297759959e+02, -1.448485621e+02, 3.295877082e+01, -3.245850428e+00, 1.179456377e-01, -2.282155654e+02, 9.249921555e+01, -1.392266984e+01, 9.297052139e-01, -2.323558576e-02},
  {73, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.306048022e+01, 1.999999910e+00, 2.646909006e+02, -1.647716545e+02, 3.719903613e+01, -3.645113853e+00, 1.319890617e-01, -2.165150972e+02, 8.722660467e+01, -1.303415548e+01, 8.633600348e-01, -2.138300143e-02},
  {74, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.308092196e+01, 2.000000008e+00, 2.251239174e+02, -1.414731209e+02, 3.206048507e+01, -3.142433101e+00, 1.135971917e-01, -2.070173544e+02, 8.296725365e+01, -1.231986936e+01, 8.102887128e-01, -1.990853407e-02},
  {75, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.309151488e+01, 1.999999984e+00, 2.627532736e+02, -1.629008146e+02, 3.661592385e+01, -3.571257833e+00, 1.286871297e-01, -1.945762063e+02, 7.740995255e+01, -1.139129234e+01, 7.415172466e-01, -1.800335280e-02},
  {76, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.310790581e+01, 2.000000068e+00, 2.644549626e+02, -1.637369900e+02, 3.675734857e+01, -3.580665992e+00, 1.288721975e-01, -1.725967865e+02, 6.755389456e+01, -9.737633351e+00, 6.184954292e-01, -1.457897448e-02},
  {77, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.311918599e+01, 1.999999933e+00, 2.677629012e+02, -1.650589135e+02, 3.690999414e+01, -3.582378706e+00, 1.284763849e-01, -1.584140848e+02, 6.122430396e+01, -8.680876005e+00, 5.402879020e-01, -1.241386995e-02},
  {78, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.319382006e+01, 2.000000009e+00, 2.420702029e+02, -1.484461630e+02, 3.292288306e+01, -3.162757529e+00, 1.121487556e-01, -1.319886050e+02, 4.940494114e+01, -6.702740089e+00, 3.934770465e-01, -8.336673895e-03},
  {79, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.320760814e+01, 1.999999969e+00, 2.346714957e+02, -1.439356552e+02, 3.189416251e+01, -3.059071523e+00, 1.082595858e-01, -1.130109430e+02, 4.093029258e+01, -5.286747014e+00, 2.885753389e-01, -5.428939868e-03},
  {80, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.316115147e+01, 2.000000093e+00, 2.747370538e+02, -1.689673404e+02, 3.771696324e+01, -3.655841153e+00, 1.309852214e-01, -9.001823908e+01, 3.066094857e+01, -3.570459523e+00, 1.613797666e-01, -1.901561361e-03},
  {81, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.313667715e+01, 2.000000002e+00, 3.142563781e+02, -1.916613838e+02, 4.259167223e+01, -4.119713271e+00, 1.474792530e-01, -7.642731867e+01, 2.462410146e+01, -2.566977318e+00, 8.741068396e-02, 1.388590928e-04},
  {82, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.314266674e+01, 1.999999876e+00, 3.509258060e+02, -2.125470710e+02, 4.702461797e+01, -4.535380912e+00, 1.620138781e-01, -5.173355302e+01, 1.362015056e+01, -7.321282362e-01, -4.826261322e-02, 3.892879264e-03},
  {83, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.315490164e+01, 1.999999944e+00, 3.399729483e+02, -2.056319770e+02, 4.539614689e+01, -4.366195994e+00, 1.554792165e-01, -4.131443229e+01, 8.986236911e+00, 3.924628986e-02, -1.052060828e-01, 5.466043586e-03},
  {84, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.316749062e+01, 1.999999920e+00, 3.640602841e+02, -2.190164327e+02, 4.815603439e+01, -4.616573783e+00, 1.639147626e-01, -3.256862965e+01, 5.115606198e+00, 6.800853161e-01, -1.522315744e-01, 6.756786448e-03},
  {85, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.318045630e+01, 2.000000044e+00, 3.766488275e+02, -2.257321142e+02, 4.947300991e+01, -4.728919006e+00, 1.674240471e-01, -2.300947210e+01, 8.615223509e-01, 1.388425307e+00, -2.045157608e-01, 8.200511055e-03},
  {86, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.319382005e+01, 2.000000006e+00, 3.443622947e+02, -2.064342780e+02, 4.516044966e+01, -4.302253084e+00, 1.516667044e-01, -5.399039282e+00, -7.002814559e+00, 2.702516748e+00, -3.018766003e-01, 1.089953798e-02},
  {87, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.291364147e+01, 2.000000217e+00, -3.706791591e+01, 1.118013187e+01, -1.057728859e+00, 3.312859839e-02, -3.138341244e-06, -3.451314336e+00, -7.779254134e+00, 2.816269849e+00, -3.090776388e-01, 1.106424389e-02},
  {88, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.284163724e+01, 1.999999954e+00, 6.125934670e+01, -4.855548659e+01, 1.248551381e+01, -1.323304763e+00, 5.060744172e-02, -6.021643455e+00, -6.580234329e+00, 2.607440108e+00, -2.929625239e-01, 1.059951856e-02},
  {89, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.285387248e+01, 2.000000090e+00, 1.350863292e+02, -9.126618691e+01, 2.169932948e+01, -2.201947573e+00, 8.186860720e-02, 1.937135880e+01, -1.787129621e+01, 4.485878662e+00, -4.315325969e-01, 1.442445798e-02},
  {90, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.286966604e+01, 1.999999970e+00, 1.784388998e+02, -1.161623817e+02, 2.702376618e+01, -2.704797298e+00, 9.957279361e-02, 2.216057166e+01, -1.904990091e+01, 4.671627339e+00, -4.444534802e-01, 1.475921763e-02},
  {91, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.287289489e+01, 1.999999951e+00, 1.368355213e+02, -9.179790820e+01, 2.169910915e+01, -2.190249857e+00, 8.102241740e-02, 4.516580666e+00, -1.118102949e+01, 3.357662550e+00, -3.470694353e-01, 1.205639951e-02},
  {92, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.287942629e+01, 2.000000032e+00, 1.427130850e+02, -9.499714618e+01, 2.234475916e+01, -2.247599931e+00, 8.291713193e-02, 1.341991149e+01, -1.518503354e+01, 4.030838171e+00, -3.972060658e-01, 1.345248084e-02},
  {93, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.288605524e+01, 1.999999761e+00, 2.341801100e+01, -2.506119713e+01, 7.023029272e+00, -7.610742531e-01, 2.903245750e-02, -3.575331738e+01, 7.276302226e+00, 1.906771859e-01, -1.059475755e-01, 5.184029625e-03},
  {94, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.288272835e+01, 1.999999941e+00, 1.287618322e+02, -8.721780968e+01, 2.073255323e+01, -2.100572716e+00, 7.794295578e-02, -2.307262580e+01, 1.113132278e+00, 1.305250601e+00, -1.948949139e-01, 7.829116438e-03},
  {95, 6.000000000e+00, 7.951499931e+00, 7.998340000e+09, -1.288940956e+01, 1.999999880e+00, 1.334821220e+02, -8.985337775e+01, 2.127928526e+01, -2.150628571e+00, 7.965294640e-02, -3.518662723e+01, 6.514543434e+00, 4.030862442e-01, -1.279850170e-01, 5.970168353e-03},
  {96, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.290553004e+01, 2.000000198e+00, 4.545581472e+01, -3.771304300e+01, 9.729129321e+00, -1.017037014e+00, 3.807733199e-02, -4.973805034e+01, 1.342335334e+01, -8.221139917e-01, -3.176841835e-02, 3.146810827e-03},
  {97, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.291150963e+01, 2.000000019e+00, 4.689042092e+01, -3.843347264e+01, 9.859294531e+00, -1.027014690e+00, 3.834833665e-02, -4.657434145e+01, 1.204637835e+01, -5.982449163e-01, -4.786919243e-02, 3.579251285e-03},
  {98, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.290833198e+01, 1.999999824e+00, 1.337584189e+01, -1.907284620e+01, 5.691614909e+00, -6.307838734e-01, 2.430868142e-02, -5.573362773e+01, 1.615667599e+01, -1.288960621e+00, 3.655033732e-03, 2.140047522e-03},
  {99, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.291435263e+01, 1.999999988e+00, 1.376201293e+01, -1.919251815e+01, 5.693799461e+00, -6.287500644e-01, 2.416045199e-02, -4.914211254e+01, 1.314247998e+01, -7.739336035e-01, -3.530513333e-02, 3.241293077e-03},
  {100, 6.000000000e+00, 8.000000000e+00, 1.000000000e+10, -1.292045700e+01, 2.000000004e+00, 1.277081775e+01, -1.854047224e+01, 5.534680382e+00, -6.118054153e-01, 2.349768815e-02, -5.074293980e+01, 1.383260974e+01, -8.858904786e-01, -2.718885953e-02, 3.019620454e-03}
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
