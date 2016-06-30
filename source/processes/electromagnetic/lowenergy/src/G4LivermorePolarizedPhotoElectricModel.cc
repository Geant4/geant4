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
// $Id: G4LivermorePhotoElectricModel.cc 90578 2015-06-04 09:55:07Z gcosmo $
//
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPhotoElectric developed by A.Forti and M.G.Pia
//
// 22 Oct 2012   A & V Ivanchenko Migration data structure to G4PhysicsVector
// 

#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LossTableManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4CrossSectionHandler.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4SauterGavrilaAngularDistribution.hh"
#include "G4AtomicShell.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LPhysicsFreeVector*  G4LivermorePolarizedPhotoElectricModel::fCrossSection[] = {nullptr};
G4LPhysicsFreeVector*  G4LivermorePolarizedPhotoElectricModel::fCrossSectionLE[] = {nullptr};
std::vector<G4double>* G4LivermorePolarizedPhotoElectricModel::fParam[] = {nullptr};
G4int                  G4LivermorePolarizedPhotoElectricModel::fNShells[] = {0};
G4int                  G4LivermorePolarizedPhotoElectricModel::fNShellsUsed[] = {0};
G4ElementData*         G4LivermorePolarizedPhotoElectricModel::fShellCrossSection = nullptr;
G4Material*            G4LivermorePolarizedPhotoElectricModel::fWater = nullptr;
G4double               G4LivermorePolarizedPhotoElectricModel::fWaterEnergyLimit = 0.0;

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::G4LivermorePolarizedPhotoElectricModel(const G4String& nam)
  : G4VEmModel(nam),fParticleChange(nullptr),maxZ(99),
    nShellLimit(100),fDeexcitationActive(false),isInitialised(false),
    fAtomDeexcitation(nullptr)
{ 
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();

  // default generator
  //  SetAngularDistribution(new G4SauterGavrilaAngularDistribution());
  SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());

  if(verboseLevel>0) {
    G4cout << "Livermore PhotoElectric is constructed " 
	   << " nShellLimit= " << nShellLimit << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fSandiaCof.resize(4,0.0);
  fCurrSection = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::~G4LivermorePolarizedPhotoElectricModel()
{  
  if(IsMaster()) {
    delete fShellCrossSection;
    for(G4int i=0; i<maxZ; ++i) { 
      delete fParam[i];
      fParam[i] = 0;
      delete fCrossSection[i];
      fCrossSection[i] = 0;
      delete fCrossSectionLE[i];
      fCrossSectionLE[i] = 0;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedPhotoElectricModel::Initialise(const G4ParticleDefinition*,
					  const G4DataVector&)
{
  if (verboseLevel > 2) {
    G4cout << "Calling G4LivermorePolarizedPhotoElectricModel::Initialise()" << G4endl;
  }

  if(IsMaster()) {

    if(!fWater) { 
      fWater = G4Material::GetMaterial("G4_WATER", false); 
      if(fWater) { fWaterEnergyLimit = 13.6*eV; }
    }

    if(!fShellCrossSection) { fShellCrossSection = new G4ElementData(); }

    char* path = getenv("G4LEDATA");

    G4ProductionCutsTable* theCoupleTable =
      G4ProductionCutsTable::GetProductionCutsTable();
    G4int numOfCouples = theCoupleTable->GetTableSize();
  
    for(G4int i=0; i<numOfCouples; ++i) {
      const G4MaterialCutsCouple* couple = 
	theCoupleTable->GetMaterialCutsCouple(i);
      const G4Material* material = couple->GetMaterial();
      const G4ElementVector* theElementVector = material->GetElementVector();
      G4int nelm = material->GetNumberOfElements();
    
      for (G4int j=0; j<nelm; ++j) {        
	G4int Z = (G4int)(*theElementVector)[j]->GetZ();
	if(Z < 1)          { Z = 1; }
	else if(Z > maxZ)  { Z = maxZ; }
	if(!fCrossSection[Z]) { ReadData(Z, path); }
      }
    }
  }  
  //  
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for LivermorePolarizedPhotoElectric model" 
	   << G4endl;
  }
  if(!isInitialised) {
    isInitialised = true;
    fParticleChange = GetParticleChangeForGamma();

    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  }
  fDeexcitationActive = false;
  if(fAtomDeexcitation) { 
    fDeexcitationActive = fAtomDeexcitation->IsFluoActive(); 
  }

  if (verboseLevel > 0) { 
    G4cout << "LivermorePolarizedPhotoElectric model is initialized " << G4endl
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedPhotoElectricModel::CrossSectionPerVolume(
                                       const G4Material* material,
				       const G4ParticleDefinition* p,
				       G4double energy,
				       G4double, G4double)
{
  fCurrSection = 0.0;
  if(fWater && (material == fWater || 
		material->GetBaseMaterial() == fWater)) {
    if(energy <= fWaterEnergyLimit) { 
      fWater->GetSandiaTable()->GetSandiaCofWater(energy, fSandiaCof);

      G4double energy2 = energy*energy;
      G4double energy3 = energy*energy2;
      G4double energy4 = energy2*energy2;

      fCurrSection = material->GetDensity()*
	(fSandiaCof[0]/energy  + fSandiaCof[1]/energy2 +
	 fSandiaCof[2]/energy3 + fSandiaCof[3]/energy4);
    } 
  }
  if(0.0 == fCurrSection) {
    fCurrSection = G4VEmModel::CrossSectionPerVolume(material, p, energy);
  }
  return fCurrSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedPhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double ZZ, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "G4LivermorePolarizedPhotoElectricModel::ComputeCrossSectionPerAtom():" 
	   << " Z= " << ZZ << "  R(keV)= " << energy/keV << G4endl;
  }
  G4double cs = 0.0;
  G4int Z = G4lrint(ZZ);
  if(Z < 1 || Z >= maxZ) { return cs; }

  // if element was not initialised
  // do initialisation safely for MT mode
  if(!fCrossSection[Z]) {
    InitialiseForElement(0, Z);
    if(!fCrossSection[Z]) { return cs; }
  }

  G4int idx = fNShells[Z]*6 - 4;
  if (energy < (*(fParam[Z]))[idx-1]) { energy = (*(fParam[Z]))[idx-1]; }
  
  G4double x1 = 1.0/energy;
  G4double x2 = x1*x1;
  G4double x3 = x2*x1;

  // parameterisation
  if(energy >= (*(fParam[Z]))[0]) {
    G4double x4 = x2*x2;
    cs = x1*((*(fParam[Z]))[idx] + x1*(*(fParam[Z]))[idx+1]
	     + x2*(*(fParam[Z]))[idx+2] + x3*(*(fParam[Z]))[idx+3] 
	     + x4*(*(fParam[Z]))[idx+4]);
    // high energy part
  } else if(energy >= (*(fParam[Z]))[1]) {
    cs = x3*(fCrossSection[Z])->Value(energy);

    // low energy part
  } else {
    cs = x3*(fCrossSectionLE[Z])->Value(energy);
  }
  if (verboseLevel > 1) { 
    G4cout << "LivermorePolarizedPhotoElectricModel: E(keV)= " << energy/keV
	   << " Z= " << Z << " cross(barn)= " << cs/barn << G4endl;
  }
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedPhotoElectricModel::SampleSecondaries(
                              std::vector<G4DynamicParticle*>* fvect,
			      const G4MaterialCutsCouple* couple,
			      const G4DynamicParticle* aDynamicGamma,
			      G4double,
			      G4double)
{
  G4double gammaEnergy = aDynamicGamma->GetKineticEnergy();
  if (verboseLevel > 3) {
    G4cout << "G4LivermorePolarizedPhotoElectricModel::SampleSecondaries() Egamma(keV)= "
	   << gammaEnergy/keV << G4endl;
  }
  
  // kill incident photon
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
  fParticleChange->SetProposedKineticEnergy(0.);

  // low-energy photo-effect in water - full absorption
  const G4Material* material = couple->GetMaterial();
  if(fWater && (material == fWater || 
		material->GetBaseMaterial() == fWater)) {
    if(gammaEnergy <= fWaterEnergyLimit) { 
      fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);      
      return;
    }
  }
 
  // Returns the normalized direction of the momentum
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection(); 

  // Select randomly one element in the current material
  //G4cout << "Select random atom Egamma(keV)= " << gammaEnergy/keV << G4endl;
  const G4Element* elm = SelectRandomAtom(material, theGamma, gammaEnergy);
  G4int Z = G4lrint(elm->GetZ());

  // Select the ionised shell in the current atom according to shell 
  //   cross sections
  // G4cout << "Select random shell Z= " << Z << G4endl;

  if(Z >= maxZ) { Z = maxZ-1; }

  // element was not initialised gamma should be absorbed
  if(!fCrossSection[Z]) {
    fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);
    return;
  }
  
  // shell index
  size_t shellIdx = 0;
  size_t nn = fNShellsUsed[Z];

  if(nn > 1) {
    if(gammaEnergy >= (*(fParam[Z]))[0]) {
      G4double x1 = 1.0/gammaEnergy;
      G4double x2 = x1*x1;
      G4double x3 = x2*x1;
      G4double x4 = x3*x1;
      G4int idx   = nn*6 - 4;
      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs0 = G4UniformRand()*((*(fParam[Z]))[idx] 
				      + x1*(*(fParam[Z]))[idx+1]
				      + x2*(*(fParam[Z]))[idx+2] 
				      + x3*(*(fParam[Z]))[idx+3] 
				      + x4*(*(fParam[Z]))[idx+4]);
      for(shellIdx=0; shellIdx<nn; ++shellIdx) {
	idx = shellIdx*6 + 2;
	if(gammaEnergy > (*(fParam[Z]))[idx-1]) {
	  G4double cs = (*(fParam[Z]))[idx] + x1*(*(fParam[Z]))[idx+1] 
	    + x2*(*(fParam[Z]))[idx+2] + x3*(*(fParam[Z]))[idx+3] 
	    + x4*(*(fParam[Z]))[idx+4];
	  if(cs >= cs0) { break; }
	}
      }
      if(shellIdx >= nn) { shellIdx = nn-1; }

    } else {

      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs = G4UniformRand();

      if(gammaEnergy >= (*(fParam[Z]))[1]) {
	cs *= (fCrossSection[Z])->Value(gammaEnergy);
      } else {
	cs *= (fCrossSectionLE[Z])->Value(gammaEnergy);
      }

      for(size_t j=0; j<nn; ++j) {
	shellIdx = (size_t)fShellCrossSection->GetComponentID(Z, j);
	if(gammaEnergy > (*(fParam[Z]))[6*shellIdx+1]) {
	  cs -= fShellCrossSection->GetValueForComponent(Z, j, gammaEnergy);
	}
	if(cs <= 0.0 || j+1 == nn) { break; }
      }
    }
  }

  G4double bindingEnergy = (*(fParam[Z]))[shellIdx*6 + 1];
  //G4cout << "Z= " << Z << " shellIdx= " << shellIdx 
  //       << " nShells= " << fNShells[Z] 
  //       << " Ebind(keV)= " << bindingEnergy/keV 
  //       << " Egamma(keV)= " << gammaEnergy/keV << G4endl;

  const G4AtomicShell* shell = 0;

  // no de-excitation from the last shell
  if(fDeexcitationActive && shellIdx + 1 < nn) {
    G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIdx);
    shell = fAtomDeexcitation->GetAtomicShell(Z, as);
  }

  // If binding energy of the selected shell is larger than photon energy
  //    do not generate secondaries
  if(gammaEnergy < bindingEnergy) {
    fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);
    return;
  }

  // Primary outcoming electron
  G4double eKineticEnergy = gammaEnergy - bindingEnergy;
  G4double edep = bindingEnergy;

  // Calculate direction of the photoelectron
  G4ThreeVector electronDirection = 
    GetAngularDistribution()->SampleDirection(aDynamicGamma, 
					      eKineticEnergy,
					      shellIdx, 
					      couple->GetMaterial());

  // The electron is created 
  G4DynamicParticle* electron = new G4DynamicParticle (theElectron,
						       electronDirection,
						       eKineticEnergy);
  fvect->push_back(electron);

  // Sample deexcitation
  if(shell) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      G4int nbefore = fvect->size();

      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      G4int nafter = fvect->size();
      if(nafter > nbefore) {
	G4double esec = 0.0;
	for (G4int j=nbefore; j<nafter; ++j) {

	  G4double e = ((*fvect)[j])->GetKineticEnergy();
	  if(esec + e > edep) {
	    // correct energy in order to have energy balance
	    e = edep - esec;
	    ((*fvect)[j])->SetKineticEnergy(e);
	    esec += e;
	    // delete the rest of secondaries (should not happens)
	    for (G4int jj=nafter-1; jj>j; --jj) { 
	      delete (*fvect)[jj]; 
	      fvect->pop_back(); 
	    }
	    break;	      
	  }
	  esec += e; 
	}
        edep -= esec;
      }
    }
  }
  // energy balance - excitation energy left
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedPhotoElectricModel::ReadData(G4int Z, const char* path)
{
  if (verboseLevel > 1) 
  {
    G4cout << "Calling ReadData() of G4LivermoreGammaConversionModel" 
           << G4endl;
  }

  if(fCrossSection[Z]) { return; }
  
  const char* datadir = path;

  if(!datadir) 
  {
    datadir = getenv("G4LEDATA");
    if(!datadir) 
    {
      G4Exception("G4LivermorePolarizedPhotoElectricModel::ReadData()",
                  "em0006",FatalException,
                  "Environment variable G4LEDATA not defined");
      return;
    }
  }

  // spline for photoeffect total x-section above K-shell
  fCrossSection[Z] = new G4LPhysicsFreeVector();
  fCrossSection[Z]->SetSpline(true);

  std::ostringstream ost;
  ost << datadir << "/livermore/phot/pe-cs-" << Z <<".dat";
  std::ifstream fin(ost.str().c_str());
  if( !fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePolarizedPhotoElectricModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePolarizedPhotoElectricModel::ReadData()",
                "em0003",FatalException,
                ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
    return;
  } else {
    if(verboseLevel > 3) { G4cout << "File " << ost.str().c_str() 
             << " is opened by G4LivermorePolarizedPhotoElectricModel" << G4endl;}
    fCrossSection[Z]->Retrieve(fin, true);
    fCrossSection[Z]->ScaleVector(MeV, barn);
    fin.close();
  }

  fParam[Z] = new std::vector<G4double>;

  // read fit parameters
  G4int n1 = 0;
  G4int n2 = 0;
  G4double x;
  std::ostringstream ost1;
  ost1 << datadir << "/livermore/phot/pe-" << Z <<".dat";
  std::ifstream fin1(ost1.str().c_str());
  if( !fin1.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePolarizedPhotoElectricModel data file <" << ost1.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePolarizedPhotoElectricModel::ReadData()",
                "em0003",FatalException,
                ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
    return;
  } else {
    if(verboseLevel > 3) { 
      G4cout << "File " << ost1.str().c_str()
             << " is opened by G4LivermorePolarizedPhotoElectricModel" << G4endl;
    }
    fin1 >> n1;
    if(fin1.fail()) { return; }
    if(0 > n1 || n1 >= INT_MAX) { n1 = 0; }

    fin1 >> n2;
    if(fin1.fail()) { return; }
    if(0 > n2 || n2 >= INT_MAX) { n2 = 0; }

    fin1 >> x;
    if(fin1.fail()) { return; }

    fNShells[Z] = n1;
    fParam[Z]->reserve(6*n1+1);
    fParam[Z]->push_back(x*MeV);
    for(G4int i=0; i<n1; ++i) {
      for(G4int j=0; j<6; ++j) {
	fin1 >> x;
        if(0 == j) { x *= MeV; }
        else       { x *= barn; }
	fParam[Z]->push_back(x);
      }
    }
    fin1.close();
  }
  // there is a possibility to used only main shells
  if(nShellLimit < n2) { n2 = nShellLimit; }
  fShellCrossSection->InitialiseForComponent(Z, n2);
  fNShellsUsed[Z] = n2;

  if(1 < n2) {
    std::ostringstream ost2;
    ost2 << datadir << "/livermore/phot/pe-ss-cs-" << Z <<".dat";
    std::ifstream fin2(ost2.str().c_str());
    if( !fin2.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedPhotoElectricModel data file <" << ost2.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4LivermorePolarizedPhotoElectricModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
      return;
    } else {
      if(verboseLevel > 3) { 
	G4cout << "File " << ost2.str().c_str()
	       << " is opened by G4LivermorePolarizedPhotoElectricModel" << G4endl;
      }

      G4int n3, n4;
      G4double y;
      for(G4int i=0; i<n2; ++i) {
	fin2 >> x >> y >> n3 >> n4;
	G4LPhysicsFreeVector* v = new G4LPhysicsFreeVector(n3, x, y);
	for(G4int j=0; j<n3; ++j) {
	  fin2 >> x >> y;
	  v->PutValues(j, x*MeV, y*barn);
	}
	fShellCrossSection->AddComponent(Z, n4, v);
      }
      fin2.close();
    }
  }

  // no spline for photoeffect total x-section below K-shell
  if(1 < fNShells[Z]) {
    fCrossSectionLE[Z] = new G4LPhysicsFreeVector();

    std::ostringstream ost3;
    ost3 << datadir << "/livermore/phot/pe-le-cs-" << Z <<".dat";
    std::ifstream fin3(ost3.str().c_str());
    if( !fin3.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePolarizedPhotoElectricModel data file <" << ost3.str().c_str()
	 << "> is not opened!" << G4endl;
      G4Exception("G4LivermorePolarizedPhotoElectricModel::ReadData()",
		  "em0003",FatalException,
		  ed,"G4LEDATA version should be G4EMLOW6.32 or later.");
      return;
    } else {
      if(verboseLevel > 3) { 
	G4cout << "File " << ost3.str().c_str() 
	       << " is opened by G4LivermorePolarizedPhotoElectricModel" << G4endl;
      }
      fCrossSectionLE[Z]->Retrieve(fin3, true);
      fCrossSectionLE[Z]->ScaleVector(MeV, barn);
      fin3.close();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoLock.hh"
namespace { G4Mutex LivermorePolarizedPhotoElectricModelMutex = G4MUTEX_INITIALIZER; }

void G4LivermorePolarizedPhotoElectricModel::InitialiseForElement(
                                    const G4ParticleDefinition*, G4int Z)
{
  G4AutoLock l(&LivermorePolarizedPhotoElectricModelMutex);
  //  G4cout << "G4LivermorePolarizedPhotoElectricModel::InitialiseForElement Z= " 
  //   << Z << G4endl;
  if(!fCrossSection[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
