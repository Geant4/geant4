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
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPhotoElectric developed by A.Forti and M.G.Pia
//
// 22 Oct 2012   A & V Ivanchenko Migration data structure to G4PhysicsVector
// 1 June 2017   M Bandieramonte: New model based on livermore/epics2014
//               evaluated data - parameterization fits in two ranges

#include "G4LivermorePhotoElectricModel.hh"

#include "G4AtomicShell.hh"
#include "G4AutoLock.hh"
#include "G4CrossSectionHandler.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4SauterGavrilaAngularDistribution.hh"
#include "G4SystemOfUnits.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4EmParameters.hh"
#include <thread>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ElementData* G4LivermorePhotoElectricModel::fCrossSection = nullptr;
G4ElementData* G4LivermorePhotoElectricModel::fCrossSectionLE = nullptr;
std::vector<G4double>* G4LivermorePhotoElectricModel::fParamHigh[] = {nullptr};
std::vector<G4double>* G4LivermorePhotoElectricModel::fParamLow[] = {nullptr};
G4int G4LivermorePhotoElectricModel::fNShells[] = {0};
G4int G4LivermorePhotoElectricModel::fNShellsUsed[] = {0};
G4Material* G4LivermorePhotoElectricModel::fWater = nullptr;
G4double G4LivermorePhotoElectricModel::fWaterEnergyLimit = 0.0;
G4String G4LivermorePhotoElectricModel::fDataDirectory = "";

static std::once_flag applyOnce;

namespace
{
  G4Mutex livPhotoeffMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::G4LivermorePhotoElectricModel(const G4String& nam) : G4VEmModel(nam)
{
  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();

  // default generator
  SetAngularDistribution(new G4SauterGavrilaAngularDistribution());

  if (verboseLevel > 0) {
    G4cout << "Livermore PhotoElectric is constructed "
           << " nShellLimit= " << nShellLimit << G4endl;
  }

  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);

  // For water
  fSandiaCof.resize(4, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::~G4LivermorePhotoElectricModel()
{
  if (isInitializer) {
    for (G4int i = 0; i < ZMAXPE; ++i) {
      if (fParamHigh[i]) {
        delete fParamHigh[i];
        fParamHigh[i] = nullptr;
      }
      if (fParamLow[i]) {
        delete fParamLow[i];
        fParamLow[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::Initialise(const G4ParticleDefinition*,
					       const G4DataVector&)
{
  if (verboseLevel > 1) {
    G4cout << "Calling G4LivermorePhotoElectricModel::Initialise() " << G4endl;
  }

  // initialise static tables only once
  std::call_once(applyOnce, [this]() { isInitializer = true; });

  if (isInitializer) {
    G4AutoLock l(&livPhotoeffMutex);
    FindDirectoryPath();
    if (fWater == nullptr) {
      fWater = G4Material::GetMaterial("G4_WATER", false);
      if (fWater == nullptr) {
        fWater = G4Material::GetMaterial("Water", false);
      }
      if (fWater != nullptr) {
        fWaterEnergyLimit = 13.6 * CLHEP::eV;
      }
    }

    if (fCrossSection == nullptr) {
      fCrossSection = new G4ElementData(ZMAXPE);
      fCrossSection->SetName("PhotoEffXS");
      fCrossSectionLE = new G4ElementData(ZMAXPE);
      fCrossSectionLE->SetName("PhotoEffLowXS");
    }

    const G4ElementTable* elemTable = G4Element::GetElementTable();
    std::size_t numElems = (*elemTable).size();
    for (std::size_t ie = 0; ie < numElems; ++ie) {
      const G4Element* elem = (*elemTable)[ie];
      G4int Z = elem->GetZasInt();
      if (Z < ZMAXPE) {
	if (fCrossSection->GetElementData(Z) == nullptr) {
	  ReadData(Z);
	}
      }
    }
    l.unlock();
  }

  if (verboseLevel > 1) {
    G4cout << "Loaded cross section files for new LivermorePhotoElectric model" << G4endl;
  }
  if (nullptr == fParticleChange) {
    fParticleChange = GetParticleChangeForGamma();
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  }

  fDeexcitationActive = false;
  if (nullptr != fAtomDeexcitation) {
    fDeexcitationActive = fAtomDeexcitation->IsFluoActive();
  }

  if (verboseLevel > 1) {
    G4cout << "LivermorePhotoElectric model is initialized " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4LivermorePhotoElectricModel::CrossSectionPerVolume(const G4Material* material,
                                                     const G4ParticleDefinition* p,
                                                     G4double energy,
						     G4double, G4double)
{
  fCurrSection = 0.0;
  if (fWater && (material == fWater || material->GetBaseMaterial() == fWater)) {
    if (energy <= fWaterEnergyLimit) {
      fWater->GetSandiaTable()->GetSandiaCofWater(energy, fSandiaCof);

      G4double energy2 = energy * energy;
      G4double energy3 = energy * energy2;
      G4double energy4 = energy2 * energy2;

      fCurrSection = material->GetDensity()
                     * (fSandiaCof[0] / energy + fSandiaCof[1] / energy2 +
			fSandiaCof[2] / energy3 + fSandiaCof[3] / energy4);
    }
  }
  if (0.0 == fCurrSection) {
    fCurrSection = G4VEmModel::CrossSectionPerVolume(material, p, energy);
  }
  return fCurrSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4LivermorePhotoElectricModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                          G4double energy,
							  G4double ZZ, G4double,
							  G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "\n G4LivermorePhotoElectricModel::ComputeCrossSectionPerAtom():"
           << " Z= " << ZZ << "  R(keV)= " << energy / keV << G4endl;
  }
  G4double cs = 0.0;
  G4int Z = G4lrint(ZZ);
  if (Z >= ZMAXPE || Z <= 0) {
    return cs;
  }
  // if element was not initialised

  // do initialisation safely for MT mode
  if (fCrossSection->GetElementData(Z) == nullptr) {
    InitialiseOnFly(Z);
    if (fCrossSection->GetElementData(Z) == nullptr) { return cs; }
  }

  // 7: rows in the parameterization file; 5: number of parameters
  G4int idx = fNShells[Z] * 7 - 5;

  energy = std::max(energy, (*(fParamHigh[Z]))[idx - 1]);

  G4double x1 = 1.0 / energy;
  G4double x2 = x1 * x1;
  G4double x3 = x2 * x1;

  // high energy parameterisation
  if (energy >= (*(fParamHigh[Z]))[0]) {
    G4double x4 = x2 * x2;
    G4double x5 = x4 * x1;

    cs = x1 *
      ((*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1]
       + x2 * (*(fParamHigh[Z]))[idx + 2] + x3 * (*(fParamHigh[Z]))[idx + 3]
       + x4 * (*(fParamHigh[Z]))[idx + 4] + x5 * (*(fParamHigh[Z]))[idx + 5]);
  }
  // low energy parameterisation
  else if (energy >= (*(fParamLow[Z]))[0]) {
    G4double x4 = x2 * x2;
    G4double x5 = x4 * x1;  // this variable usage can probably be optimized
    cs = x1 *
      ((*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1]
       + x2 * (*(fParamLow[Z]))[idx + 2] + x3 * (*(fParamLow[Z]))[idx + 3]
       + x4 * (*(fParamLow[Z]))[idx + 4] + x5 * (*(fParamLow[Z]))[idx + 5]);
  }
  // Tabulated values above k-shell ionization energy
  else if (energy >= (*(fParamHigh[Z]))[1]) {
    cs = x3 * (fCrossSection->GetElementData(Z))->Value(energy);
  }
  // Tabulated values below k-shell ionization energy
  else {
    cs = x3 * (fCrossSectionLE->GetElementData(Z))->Value(energy);
  }
  if (verboseLevel > 1) {
    G4cout << "G4LivermorePhotoElectricModel: E(keV)= " << energy / keV << " Z= " << Z
           << " cross(barn)= " << cs / barn << G4endl;
  }
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                      const G4MaterialCutsCouple* couple,
                                                      const G4DynamicParticle* aDynamicGamma,
                                                      G4double, G4double)
{
  G4double gammaEnergy = aDynamicGamma->GetKineticEnergy();
  if (verboseLevel > 3) {
    G4cout << "G4LivermorePhotoElectricModel::SampleSecondaries() Egamma(keV)= "
           << gammaEnergy / keV << G4endl;
  }

  // kill incident photon
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->SetProposedKineticEnergy(0.);

  // low-energy photo-effect in water - full absorption
  const G4Material* material = couple->GetMaterial();
  if (fWater && (material == fWater || material->GetBaseMaterial() == fWater)) {
    if (gammaEnergy <= fWaterEnergyLimit) {
      fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);
      return;
    }
  }

  // Returns the normalized direction of the momentum
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  const G4Element* elm = SelectRandomAtom(material, theGamma, gammaEnergy);
  G4int Z = elm->GetZasInt();

  // Select the ionised shell in the current atom according to shell
  // cross sections

  // If element was not initialised gamma should be absorbed
  if (Z >= ZMAXPE || fCrossSection->GetElementData(Z) == nullptr) {
    fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);
    return;
  }

  // SAMPLING OF THE SHELL INDEX
  std::size_t shellIdx = 0;
  std::size_t nn = fNShellsUsed[Z];
  if (nn > 1) {
    if (gammaEnergy >= (*(fParamHigh[Z]))[0]) {
      G4double x1 = 1.0 / gammaEnergy;
      G4double x2 = x1 * x1;
      G4double x3 = x2 * x1;
      G4double x4 = x3 * x1;
      G4double x5 = x4 * x1;
      std::size_t idx = nn * 7 - 5;
      // when do sampling common factors are not taken into account
      // so cross section is not real

      G4double rand = G4UniformRand();
      G4double cs0 = rand *
        ((*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1]
        + x2 * (*(fParamHigh[Z]))[idx + 2] + x3 * (*(fParamHigh[Z]))[idx + 3]
        + x4 * (*(fParamHigh[Z]))[idx + 4] + x5 * (*(fParamHigh[Z]))[idx + 5]);

      for (shellIdx = 0; shellIdx < nn; ++shellIdx) {
        idx = shellIdx * 7 + 2;
        if (gammaEnergy > (*(fParamHigh[Z]))[idx - 1]) {
          G4double cs = (*(fParamHigh[Z]))[idx] + x1 * (*(fParamHigh[Z]))[idx + 1]
            + x2 * (*(fParamHigh[Z]))[idx + 2] + x3 * (*(fParamHigh[Z]))[idx + 3]
            + x4 * (*(fParamHigh[Z]))[idx + 4] + x5 * (*(fParamHigh[Z]))[idx + 5];

          if (cs >= cs0) {
            break;
          }
        }
      }
      if (shellIdx >= nn) {
        shellIdx = nn - 1;
      }
    }
    else if (gammaEnergy >= (*(fParamLow[Z]))[0]) {
      G4double x1 = 1.0 / gammaEnergy;
      G4double x2 = x1 * x1;
      G4double x3 = x2 * x1;
      G4double x4 = x3 * x1;
      G4double x5 = x4 * x1;
      std::size_t idx = nn * 7 - 5;
      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs0 = G4UniformRand() *
        ((*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1]
        + x2 * (*(fParamLow[Z]))[idx + 2] + x3 * (*(fParamLow[Z]))[idx + 3]
        + x4 * (*(fParamLow[Z]))[idx + 4] + x5 * (*(fParamLow[Z]))[idx + 5]);
      for (shellIdx = 0; shellIdx < nn; ++shellIdx) {
        idx = shellIdx * 7 + 2;
        if (gammaEnergy > (*(fParamLow[Z]))[idx - 1]) {
          G4double cs = (*(fParamLow[Z]))[idx] + x1 * (*(fParamLow[Z]))[idx + 1]
            + x2 * (*(fParamLow[Z]))[idx + 2] + x3 * (*(fParamLow[Z]))[idx + 3]
            + x4 * (*(fParamLow[Z]))[idx + 4] + x5 * (*(fParamLow[Z]))[idx + 5];
          if (cs >= cs0) {
            break;
          }
        }
      }
      if (shellIdx >= nn) {
        shellIdx = nn - 1;
      }
    }
    else {
      // when do sampling common factors are not taken into account
      // so cross section is not real
      G4double cs = G4UniformRand();

      if (gammaEnergy >= (*(fParamHigh[Z]))[1]) {
        // above K-shell binding energy
        cs *= fCrossSection->GetElementData(Z)->Value(gammaEnergy);
      }
      else {
        // below K-shell binding energy
        cs *= fCrossSectionLE->GetElementData(Z)->Value(gammaEnergy);
      }

      for (G4int j = 0; j < (G4int)nn; ++j) {
        shellIdx = (std::size_t)fCrossSection->GetComponentID(Z, j);
        if (gammaEnergy > (*(fParamLow[Z]))[7 * shellIdx + 1]) {
          cs -= fCrossSection->GetValueForComponent(Z, j, gammaEnergy);
        }
        if (cs <= 0.0 || j + 1 == (G4int)nn) {
          break;
        }
      }
    }
  }
  // END: SAMPLING OF THE SHELL

  G4double bindingEnergy = (*(fParamHigh[Z]))[shellIdx * 7 + 1];
  const G4AtomicShell* shell = nullptr;

  // no de-excitation from the last shell
  if (fDeexcitationActive && shellIdx + 1 < nn) {
    auto as = G4AtomicShellEnumerator(shellIdx);
    shell = fAtomDeexcitation->GetAtomicShell(Z, as);
  }

  // If binding energy of the selected shell is larger than photon energy
  //    do not generate secondaries
  if (gammaEnergy < bindingEnergy) {
    fParticleChange->ProposeLocalEnergyDeposit(gammaEnergy);
    return;
  }

  // Primary outcoming electron
  G4double eKineticEnergy = gammaEnergy - bindingEnergy;
  G4double edep = bindingEnergy;

  // Calculate direction of the photoelectron
  G4ThreeVector electronDirection = GetAngularDistribution()->SampleDirection(
    aDynamicGamma, eKineticEnergy, (G4int)shellIdx, couple->GetMaterial());

  // The electron is created
  auto electron =
    new G4DynamicParticle(theElectron, electronDirection, eKineticEnergy);
  fvect->push_back(electron);

  // Sample deexcitation
  if (nullptr != shell) {
    G4int index = couple->GetIndex();
    if (fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      std::size_t nbefore = fvect->size();

      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      std::size_t nafter = fvect->size();
      if (nafter > nbefore) {
        G4double esec = 0.0;
        for (std::size_t j = nbefore; j < nafter; ++j) {
          G4double e = ((*fvect)[j])->GetKineticEnergy();
          if (esec + e > edep) {
            // correct energy in order to have energy balance
            e = edep - esec;
            ((*fvect)[j])->SetKineticEnergy(e);
            esec += e;
            // delete the rest of secondaries (should not happens)
            for (std::size_t jj = nafter - 1; jj > j; --jj) {
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
  if (edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& G4LivermorePhotoElectricModel::FindDirectoryPath()
{
  // no check in this method - environment variable is check by utility
  if (fDataDirectory.empty()) {
    auto param = G4EmParameters::Instance();
    std::ostringstream ost;
    if (param->LivermoreDataDir() == "livermore") {
      ost << param->GetDirLEDATA() << "/livermore/phot_epics2014/";
    }
    else {
      ost << param->GetDirLEDATA() << "/epics2017/phot/";
    }
    fDataDirectory = ost.str();
  }
  return fDataDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::ReadData(G4int Z)
{
  if (Z <= 0 || Z >= ZMAXPE) {
    G4cout << "G4LivermorePhotoElectricModel::ReadData() Warning: Z="
	   << Z << " is out of range - request ignored" << G4endl;
    return;
  }
  if (verboseLevel > 1) {
    G4cout << "G4LivermorePhotoElectricModel::ReadData() for Z=" << Z << G4endl;
  }

  if (fCrossSection->GetElementData(Z) != nullptr) {
    return;
  }

  // spline for photoeffect total x-section above K-shell when using EPDL97
  // but below the parameterized ones

  G4bool spline = (G4EmParameters::Instance()->LivermoreDataDir() == "livermore");
  G4int number = G4EmParameters::Instance()->NumberForFreeVector();
  auto pv = new G4PhysicsFreeVector(spline);

  // fDataDirectory will be defined after these lines
  std::ostringstream ost;
  ost << FindDirectoryPath() << "pe-cs-" << Z << ".dat";
  std::ifstream fin(ost.str().c_str());
  if (!fin.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePhotoElectricModel data file <"
       << ost.str().c_str() << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePhotoElectricModel::ReadData()", "em0003", FatalException, ed,
                "G4LEDATA version should be G4EMLOW8.0 or later.");
    return;
  }
  if (verboseLevel > 3) {
    G4cout << "File " << ost.str().c_str() << " is opened by G4LivermorePhotoElectricModel"
           << G4endl;
  }
  pv->Retrieve(fin, true);
  pv->ScaleVector(MeV, barn);
  pv->FillSecondDerivatives();
  pv->EnableLogBinSearch(number);
  fCrossSection->InitialiseForElement(Z, pv);
  fin.close();

  // read high-energy fit parameters
  fParamHigh[Z] = new std::vector<G4double>;
  G4int n1 = 0;
  G4int n2 = 0;
  G4double x;
  std::ostringstream ost1;
  ost1 << fDataDirectory << "pe-high-" << Z << ".dat";
  std::ifstream fin1(ost1.str().c_str());
  if (!fin1.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePhotoElectricModel data file <"
       << ost1.str().c_str() << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePhotoElectricModel::ReadData()", "em0003", FatalException, ed,
                "G4LEDATA version should be G4EMLOW7.2 or later.");
    return;
  }
  if (verboseLevel > 3) {
    G4cout << "File " << ost1.str().c_str() << " is opened by G4LivermorePhotoElectricModel"
           << G4endl;
  }
  fin1 >> n1;
  if (fin1.fail()) {
    return;
  }
  if (0 > n1 || n1 >= INT_MAX) {
    n1 = 0;
  }

  fin1 >> n2;
  if (fin1.fail()) {
    return;
  }
  if (0 > n2 || n2 >= INT_MAX) {
    n2 = 0;
  }

  fin1 >> x;
  if (fin1.fail()) {
    return;
  }

  fNShells[Z] = n1;
  fParamHigh[Z]->reserve(7 * n1 + 1);
  fParamHigh[Z]->push_back(x * MeV);
  for (G4int i = 0; i < n1; ++i) {
    for (G4int j = 0; j < 7; ++j) {
      fin1 >> x;
      if (0 == j) {
        x *= MeV;
      }
      else {
        x *= barn;
      }
      fParamHigh[Z]->push_back(x);
    }
  }
  fin1.close();

  // read low-energy fit parameters
  fParamLow[Z] = new std::vector<G4double>;
  G4int n1_low = 0;
  G4int n2_low = 0;
  G4double x_low;
  std::ostringstream ost1_low;
  ost1_low << fDataDirectory << "pe-low-" << Z << ".dat";
  std::ifstream fin1_low(ost1_low.str().c_str());
  if (!fin1_low.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermorePhotoElectricModel data file <" << ost1_low.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermorePhotoElectricModel::ReadData()", "em0003", FatalException, ed,
                "G4LEDATA version should be G4EMLOW8.0 or later.");
    return;
  }
  if (verboseLevel > 3) {
    G4cout << "File " << ost1_low.str().c_str() << " is opened by G4LivermorePhotoElectricModel"
           << G4endl;
  }
  fin1_low >> n1_low;
  if (fin1_low.fail()) {
    return;
  }
  if (0 > n1_low || n1_low >= INT_MAX) {
    n1_low = 0;
  }

  fin1_low >> n2_low;
  if (fin1_low.fail()) {
    return;
  }
  if (0 > n2_low || n2_low >= INT_MAX) {
    n2_low = 0;
  }

  fin1_low >> x_low;
  if (fin1_low.fail()) {
    return;
  }

  fNShells[Z] = n1_low;
  fParamLow[Z]->reserve(7 * n1_low + 1);
  fParamLow[Z]->push_back(x_low * MeV);
  for (G4int i = 0; i < n1_low; ++i) {
    for (G4int j = 0; j < 7; ++j) {
      fin1_low >> x_low;
      if (0 == j) {
        x_low *= MeV;
      }
      else {
        x_low *= barn;
      }
      fParamLow[Z]->push_back(x_low);
    }
  }
  fin1_low.close();

  // there is a possibility to use only main shells
  if (nShellLimit < n2) {
    n2 = nShellLimit;
  }
  fCrossSection->InitialiseForComponent(Z, n2);  // number of shells
  fNShellsUsed[Z] = n2;

  if (1 < n2) {
    std::ostringstream ost2;
    ost2 << fDataDirectory << "pe-ss-cs-" << Z << ".dat";
    std::ifstream fin2(ost2.str().c_str());
    if (!fin2.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePhotoElectricModel data file <" << ost2.str().c_str() << "> is not opened!"
         << G4endl;
      G4Exception("G4LivermorePhotoElectricModel::ReadData()", "em0003", FatalException, ed,
                  "G4LEDATA version should be G4EMLOW7.2 or later.");
      return;
    }
    if (verboseLevel > 3) {
      G4cout << "File " << ost2.str().c_str() << " is opened by G4LivermorePhotoElectricModel"
             << G4endl;
    }

    G4int n3, n4;
    G4double y;
    for (G4int i = 0; i < n2; ++i) {
      fin2 >> x >> y >> n3 >> n4;
      auto v = new G4PhysicsFreeVector(n3, x, y);
      for (G4int j = 0; j < n3; ++j) {
        fin2 >> x >> y;
        v->PutValues(j, x * CLHEP::MeV, y * CLHEP::barn);
      }
      v->EnableLogBinSearch(number);
      fCrossSection->AddComponent(Z, n4, v);
    }
    fin2.close();
  }

  // no spline for photoeffect total x-section below K-shell
  if (1 < fNShells[Z]) {
    auto pv1 = new G4PhysicsFreeVector(false);
    std::ostringstream ost3;
    ost3 << fDataDirectory << "pe-le-cs-" << Z << ".dat";
    std::ifstream fin3(ost3.str().c_str());
    if (!fin3.is_open()) {
      G4ExceptionDescription ed;
      ed << "G4LivermorePhotoElectricModel data file <" << ost3.str().c_str() << "> is not opened!"
         << G4endl;
      G4Exception("G4LivermorePhotoElectricModel::ReadData()", "em0003", FatalException, ed,
                  "G4LEDATA version should be G4EMLOW8.0 or later.");
      return;
    }
    if (verboseLevel > 3) {
      G4cout << "File " << ost3.str().c_str() << " is opened by G4LivermorePhotoElectricModel"
             << G4endl;
    }
    pv1->Retrieve(fin3, true);
    pv1->ScaleVector(CLHEP::MeV, CLHEP::barn);
    pv1->EnableLogBinSearch(number);
    fCrossSectionLE->InitialiseForElement(Z, pv1);
    fin3.close();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePhotoElectricModel::GetBindingEnergy(G4int Z, G4int shell)
{
  if (Z < 1 || Z >= ZMAXPE) {
    return -1;
  }  // If Z is out of the supported return 0

  InitialiseOnFly(Z);
  if (fCrossSection->GetElementData(Z) == nullptr ||
      shell < 0 || shell >= fNShellsUsed[Z]) {
    return -1;
  }

  if (Z > 2) {
    return fCrossSection->GetComponentDataByIndex(Z, shell)->Energy(0);
  }
  else {
    return fCrossSection->GetElementData(Z)->Energy(0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4LivermorePhotoElectricModel::InitialiseForElement(const G4ParticleDefinition*, G4int Z)
{
  if (fCrossSection == nullptr) {
    fCrossSection = new G4ElementData(ZMAXPE);
    fCrossSection->SetName("PhotoEffXS");
    fCrossSectionLE = new G4ElementData(ZMAXPE);
    fCrossSectionLE->SetName("PhotoEffLowXS");
  }
  ReadData(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LivermorePhotoElectricModel::InitialiseOnFly(G4int Z)
{
  if (fCrossSection->GetElementData(Z) == nullptr && Z > 0 && Z < ZMAXPE) {
    G4AutoLock l(&livPhotoeffMutex);
    if (fCrossSection->GetElementData(Z) == nullptr) {
      ReadData(Z);
    }
    l.unlock();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
