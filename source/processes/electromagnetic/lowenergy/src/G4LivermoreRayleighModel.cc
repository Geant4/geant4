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
// Author: Sebastien Incerti
//         31 March 2012
//         on base of G4LivermoreRayleighModel
//

#include "G4LivermoreRayleighModel.hh"

#include "G4AutoLock.hh"
#include "G4EmParameters.hh"
#include "G4RayleighAngularGenerator.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace
{
G4Mutex LivermoreRayleighModelMutex = G4MUTEX_INITIALIZER;
}

G4PhysicsFreeVector* G4LivermoreRayleighModel::dataCS[] = {nullptr};
G4String G4LivermoreRayleighModel::gDataDirectory = "";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreRayleighModel::G4LivermoreRayleighModel() : G4VEmModel("LivermoreRayleigh")
{
  fParticleChange = nullptr;
  lowEnergyLimit = 10 * CLHEP::eV;

  SetAngularDistribution(new G4RayleighAngularGenerator());

  verboseLevel = 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods

  if (verboseLevel > 0) {
    G4cout << "G4LivermoreRayleighModel is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreRayleighModel::~G4LivermoreRayleighModel()
{
  if (IsMaster()) {
    for (G4int i = 0; i <= maxZ; ++i) {
      if (nullptr != dataCS[i]) {
        delete dataCS[i];
        dataCS[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& cuts)
{
  if (verboseLevel > 1) {
    G4cout << "Calling Initialise() of G4LivermoreRayleighModel." << G4endl
           << "Energy range: " << LowEnergyLimit() / eV << " eV - " << HighEnergyLimit() / GeV
           << " GeV" << G4endl;
  }

  if (IsMaster()) {
    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);

    // Access to elements
    const G4ElementTable* elemTable = G4Element::GetElementTable();
    std::size_t numElems = (*elemTable).size();
    for (std::size_t ie = 0; ie < numElems; ++ie) {
      const G4Element* elem = (*elemTable)[ie];
      const G4int Z = std::min(maxZ, elem->GetZasInt());
      if (dataCS[Z] == nullptr) {
        ReadData(Z);
      }
    }
  }
  if (isInitialised) {
    return;
  }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::InitialiseLocal(const G4ParticleDefinition*, G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4String& G4LivermoreRayleighModel::FindDirectoryPath()
{
  // no check in this method - environment variable is check by utility
  if (gDataDirectory.empty()) {
    auto param = G4EmParameters::Instance();
    std::ostringstream ost;
    if (param->LivermoreDataDir() == "livermore") {
      ost << param->GetDirLEDATA() << "/livermore/rayl/";
    }
    else {
      ost << param->GetDirLEDATA() << "/epics2017/rayl/";
    }
    gDataDirectory = ost.str();
  }
  return gDataDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::ReadData(const G4int ZZ)
{
  if (verboseLevel > 1) {
    G4cout << "Calling ReadData() of G4LivermoreRayleighModel for Z=" << ZZ << G4endl;
  }
  const G4int Z = std::min(ZZ, maxZ);

  if (nullptr != dataCS[Z]) {
    return;
  }

  dataCS[Z] = new G4PhysicsFreeVector();

  std::ostringstream ostCS;
  ostCS << FindDirectoryPath() << "re-cs-" << Z << ".dat";

  std::ifstream finCS(ostCS.str().c_str());

  if (!finCS.is_open()) {
    G4ExceptionDescription ed;
    ed << "G4LivermoreRayleighModel data file <" << ostCS.str().c_str() << "> is not opened!"
       << G4endl;
    G4Exception("G4LivermoreRayleighModel::ReadData()", "em0003", FatalException, ed,
                "G4LEDATA version should be G4EMLOW8.0 or later.");
    return;
  }
  else {
    if (verboseLevel > 3) {
      G4cout << "File " << ostCS.str() << " is opened by G4LivermoreRayleighModel" << G4endl;
    }
    dataCS[Z]->Retrieve(finCS, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreRayleighModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                              G4double GammaEnergy, G4double Z,
                                                              G4double, G4double, G4double)
{
  if (verboseLevel > 1) {
    G4cout << "G4LivermoreRayleighModel::ComputeCrossSectionPerAtom()" << G4endl;
  }

  if (GammaEnergy < lowEnergyLimit) {
    return 0.0;
  }

  G4double xs = 0.0;
  G4int intZ = G4lrint(Z);
  if (intZ < 1 || intZ > maxZ) {
    return xs;
  }

  G4PhysicsFreeVector* pv = dataCS[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if (nullptr == pv) {
    InitialiseForElement(nullptr, intZ);
    pv = dataCS[intZ];
    if (nullptr == pv) {
      return xs;
    }
  }

  auto n = G4int(pv->GetVectorLength() - 1);
  G4double e = GammaEnergy / MeV;
  if (e >= pv->Energy(n)) {
    xs = (*pv)[n] / (e * e);
  }
  else if (e >= pv->Energy(0)) {
    xs = pv->Value(e) / (e * e);
  }

  if (verboseLevel > 1) {
    G4cout << "****** DEBUG: tcs value for Z=" << Z << " at energy (MeV)=" << e << G4endl;
    G4cout << "  cs (Geant4 internal unit)=" << xs << G4endl;
    G4cout << "    -> first E*E*cs value in CS data file (iu) =" << (*pv)[0] << G4endl;
    G4cout << "    -> last  E*E*cs value in CS data file (iu) =" << (*pv)[n] << G4endl;
    G4cout << "*********************************************************" << G4endl;
  }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                                 const G4MaterialCutsCouple* couple,
                                                 const G4DynamicParticle* aDynamicGamma, G4double,
                                                 G4double)
{
  if (verboseLevel > 1) {
    G4cout << "Calling SampleSecondaries() of G4LivermoreRayleighModel" << G4endl;
  }
  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  // Select randomly one element in the current material
  const G4ParticleDefinition* particle = aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple, particle, photonEnergy0);
  G4int Z = elm->GetZasInt();

  // Sample the angle of the scattered photon
  G4ThreeVector photonDirection = GetAngularDistribution()->SampleDirection(
    aDynamicGamma, photonEnergy0, Z, couple->GetMaterial());
  fParticleChange->ProposeMomentumDirection(photonDirection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LivermoreRayleighModel::InitialiseForElement(const G4ParticleDefinition*, G4int Z)
{
  if (nullptr != dataCS[Z]) {
    return;
  }
  G4AutoLock l(&LivermoreRayleighModelMutex);
  if (nullptr == dataCS[Z]) {
    ReadData(Z);
  }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
