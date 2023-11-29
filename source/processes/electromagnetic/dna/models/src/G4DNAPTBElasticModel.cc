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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Models come from
// M. Bug et al, Rad. Phys and Chem. 130, 459-479 (2017)
//

#include "G4DNAPTBElasticModel.hh"

#include "G4DNAChampionElasticModel.hh"
#include "G4DNAMaterialManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"

G4DNAPTBElasticModel::G4DNAPTBElasticModel(const G4String& applyToMaterial,
                                           const G4ParticleDefinition*, const G4String& nam)
  : G4VDNAModel(nam, applyToMaterial)
{
  if (verboseLevel > 0) {
    G4cout << "PTB Elastic model is constructed : " << G4endl;
  }
  fpTHF = G4Material::GetMaterial("THF", false);
  fpPY = G4Material::GetMaterial("PY", false);
  fpPU = G4Material::GetMaterial("PU", false);
  fpTMP = G4Material::GetMaterial("TMP", false);
  fpG4_WATER = G4Material::GetMaterial("G4_WATER", false);
  fpBackbone_THF = G4Material::GetMaterial("backbone_THF", false);
  fpCytosine_PY = G4Material::GetMaterial("cytosine_PY", false);
  fpThymine_PY = G4Material::GetMaterial("thymine_PY", false);
  fpAdenine_PU = G4Material::GetMaterial("adenine_PU", false);
  fpBackbone_TMP = G4Material::GetMaterial("backbone_TMP", false);
  fpGuanine_PU = G4Material::GetMaterial("guanine_PU", false);
  fpN2 = G4Material::GetMaterial("N2", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::Initialise(const G4ParticleDefinition* particle,
                                      const G4DataVector& /*cuts*/)
{
  if (isInitialised) {
    return;
  }
  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNAPTBElasticModel::Initialise()" << G4endl;
  }
  if (particle != G4Electron::ElectronDefinition()) {
    std::ostringstream oss;
    oss << " Model is not applied for this particle " << particle->GetParticleName();
    G4Exception("G4DNAPTBElasticModel::G4DNAPTBElasticModel", "PTB001", FatalException,
                oss.str().c_str());
  }
  G4double scaleFactor = 1e-16 * cm * cm;
  //*******************************************************
  // Cross section data
  //*******************************************************

  std::size_t index;
  // MPietrzak, adding paths for N2
  if (fpN2 != nullptr) {
    index = fpN2->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_N2",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_N2", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1.02 * MeV);
  }
  // MPietrzak

  if (fpTHF != nullptr) {
    index = fpTHF->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_THF",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_THF", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpPY != nullptr) {
    index = fpPY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PY",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PY", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpPU != nullptr) {
    index = fpPU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PU",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PU", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpTMP != nullptr) {
    index = fpTMP->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_TMP",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_TMP", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }
  //????
  if (fpG4_WATER != nullptr) {
    index = fpG4_WATER->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e_champion",
                        "dna/sigmadiff_cumulated_elastic_e_champion", scaleFactor);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }
  // DNA materials
  //
  if (fpBackbone_THF != nullptr) {
    index = fpBackbone_THF->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_THF",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_THF", scaleFactor * 33. / 30);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpCytosine_PY != nullptr) {
    index = fpCytosine_PY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PY",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PY", scaleFactor * 42. / 30);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpThymine_PY != nullptr) {
    index = fpThymine_PY->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PY",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PY", scaleFactor * 48. / 30);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpAdenine_PU != nullptr) {
    index = fpAdenine_PU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PU",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PU", scaleFactor * 50. / 44);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }
  if (fpGuanine_PU != nullptr) {
    index = fpGuanine_PU->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_PU",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_PU", scaleFactor * 56. / 44);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (fpBackbone_TMP != nullptr) {
    index = fpBackbone_TMP->GetIndex();
    AddCrossSectionData(index, particle, "dna/sigma_elastic_e-_PTB_TMP",
                        "dna/sigmadiff_cumulated_elastic_e-_PTB_TMP", scaleFactor * 33. / 50);
    SetLowELimit(index, particle, 10 * eV);
    SetHighELimit(index, particle, 1 * keV);
  }

  if (!G4DNAMaterialManager::Instance()->IsLocked()) {
    // Load the data
    LoadCrossSectionData(particle);
    G4DNAMaterialManager::Instance()->SetMasterDataModel(DNAModelType::fDNAElastics, this);
    fpModelData = this;
  }
  else {
    auto dataModel = dynamic_cast<G4DNAPTBElasticModel*>(
      G4DNAMaterialManager::Instance()->GetModel(DNAModelType::fDNAElastics));
    if (dataModel == nullptr) {
      G4cout << "G4DNAPTBElasticModel::Initialise:: not good modelData" << G4endl;
      G4Exception("G4DNAPTBElasticModel::Initialise", "PTB0006", FatalException,
                  "not good modelData");
    }
    else {
      fpModelData = dataModel;
    }
  }

  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for PTB Elastic model" << G4endl;
  }

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::ReadDiffCSFile(const std::size_t& materialName,
                                          const G4ParticleDefinition* particleName,
                                          const G4String& file, const G4double&)
{
  // Method to read and save the information contained within the differential cross section files.
  // This method is not yet standard.

  // get the path of the G4LEDATA data folder
  const char* path = G4FindDataDir("G4LEDATA");
  // if it is not found then quit and print error message
  if (!path) {
    G4Exception("G4DNAPTBElasticModel::ReadAllDiffCSFiles", "em0006", FatalException,
                "G4LEDATA environment variable not set.");
    return;
  }

  // build the fullFileName path of the data file
  std::ostringstream fullFileName;
  fullFileName << path << "/" << file << ".dat";

  // open the data file
  std::ifstream diffCrossSection(fullFileName.str().c_str());
  // error if file is not there
  std::stringstream endPath;
  if (!diffCrossSection) {
    endPath << "Missing data file: " << file;
    G4Exception("G4DNAPTBElasticModel::Initialise", "em0003", FatalException,
                endPath.str().c_str());
  }

  tValuesVec[materialName][particleName].push_back(0.);

  G4String line;

  // read the file line by line until we reach the end of file point
  while (std::getline(diffCrossSection, line)) {
    // check if the line is comment or empty
    //
    std::istringstream testIss(line);
    G4String test;
    testIss >> test;
    // check first caracter to determine if following information is data or comments
    if (test == "#") {
      // skip the line by beginning a new while loop.
      continue;
    }
    // check if line is empty
    else if (line.empty()) {
      // skip the line by beginning a new while loop.
      continue;
    }
    //
    // end of the check

    // transform the line into a iss
    std::istringstream iss(line);

    // Variables to be filled by the input file
    G4double tDummy;
    G4double eDummy;

    // fill the variables with the content of the line
    iss >> tDummy >> eDummy;

    // SI : mandatory Vecm initialization

    // Fill two vectors contained in maps of types:
    // [materialName][particleName]=vector
    // [materialName][particleName][T]=vector
    // to list all the incident energies (tValues) and all the output energies (eValues) within the
    // file
    //
    // Check if we already have the current T value in the vector.
    // If not then add it
    if (tDummy != tValuesVec[materialName][particleName].back()) {
      // Add the current T value
      tValuesVec[materialName][particleName].push_back(tDummy);
      // Make it correspond to a default zero E value
      eValuesVect[materialName][particleName][tDummy].push_back(0.);
    }

    // Put the differential cross section value of the input file within the diffCrossSectionData
    // map
    iss >> diffCrossSectionData[materialName][particleName][tDummy][eDummy];

    // If the current E value (eDummy) is different from the one already registered in the eVector
    // then add it to the vector
    if (eDummy != eValuesVect[materialName][particleName][tDummy].back()) {
      eValuesVect[materialName][particleName][tDummy].push_back(eDummy);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::CrossSectionPerVolume(const G4Material* pMaterial,
                                                     const G4ParticleDefinition* p, G4double ekin,
                                                     G4double /*emin*/, G4double /*emax*/)
{
  if (verboseLevel > 3){
    G4cout << "Calling CrossSectionPerVolume() of G4DNAPTBElasticModel" << G4endl;
  }

  // Get the name of the current particle
  const G4String& particleName = p->GetParticleName();
  const std::size_t& materialID = pMaterial->GetIndex();

  // set killBelowEnergy value for current material
  fKillBelowEnergy = fpModelData->GetLowELimit(materialID, p);
  // initialise the return value (cross section) to zero
  G4double sigma = 0.;

  // check if we are below the high energy limit
  if (ekin < fpModelData->GetHighELimit(materialID, p)) {
    // This is used to kill the particle if its kinetic energy is below fKillBelowEnergy.
    // If the energy is lower then we return a maximum cross section and thus the SampleSecondaries
    // method will be called for sure. SampleSecondaries will remove the particle from the
    // simulation.
    //
    // SI : XS must not be zero otherwise sampling of secondaries method ignored
    if (ekin < fKillBelowEnergy) {
      return DBL_MAX;
    }

    // Get the tables with the cross section data
    auto tableData = fpModelData->GetData();
    if ((*tableData)[materialID][p] == nullptr) {
      G4Exception("G4DNAPTBElasticModel::CrossSectionPerVolume", "em00236", FatalException,
                  "No model is registered");
    }
    // Retrieve the cross section value
    sigma = (*tableData)[materialID][p]->FindValue(ekin);
  }

  if (verboseLevel > 2) {
    G4cout << "__________________________________" << G4endl;
    G4cout << "°°° G4DNAPTBElasticModel - XS INFO START" << G4endl;
    G4cout << "°°° Kinetic energy(eV)=" << ekin / eV << " particle : " << particleName << G4endl;
    G4cout << "°°° Cross section per molecule (cm^2)=" << sigma / cm / cm << G4endl;
    G4cout << "°°° G4DNAPTBElasticModel - XS INFO END" << G4endl;
  }

  // Return the cross section
  auto MolDensity =
    (*G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(pMaterial))[materialID];
  return sigma * MolDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAPTBElasticModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                             const G4MaterialCutsCouple* couple,
                                             const G4DynamicParticle* aDynamicElectron,
                                             G4double /*tmin*/, G4double /*tmax*/)
{
  if (verboseLevel > 3) {
    G4cout << "Calling SampleSecondaries() of G4DNAPTBElasticModel" << G4endl;
  }

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  const std::size_t& materialID = couple->GetIndex();
  auto p = aDynamicElectron->GetParticleDefinition();

  // set killBelowEnergy value for material
  fKillBelowEnergy = fpModelData->GetLowELimit(materialID, p);

  // If the particle (electron here) energy is below the kill limit then we remove it from the
  // simulation
  if (electronEnergy0 < fKillBelowEnergy) {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
  }
  // If we are above the kill limite and below the high limit then we proceed
  else if (electronEnergy0 >= fKillBelowEnergy && electronEnergy0 < GetHighELimit(materialID, p)) {
    // Random sampling of the cosTheta
    G4double cosTheta = fpModelData->RandomizeCosTheta(electronEnergy0, materialID);

    // Random sampling of phi
    G4double phi = 2. * CLHEP::pi * G4UniformRand();

    auto zVers = aDynamicElectron->GetMomentumDirection();
    auto xVers = zVers.orthogonal();
    auto yVers = zVers.cross(xVers);

    G4double xDir = std::sqrt(1. - cosTheta * cosTheta);
    G4double yDir = xDir;
    xDir *= std::cos(phi);
    yDir *= std::sin(phi);

    // Particle direction after ModelInterface
    G4ThreeVector zPrikeVers((xDir * xVers + yDir * yVers + cosTheta * zVers));

    // Give the new direction
    fParticleChangeForGamma->ProposeMomentumDirection(zPrikeVers.unit());

    // Update the energy which does not change here
    fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::Theta(const G4ParticleDefinition* p, G4double k, G4double integrDiff,
                                     const std::size_t& materialID)
{
  G4double theta = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;
  G4double xs12 = 0;
  G4double xs21 = 0;
  G4double xs22 = 0;
  if (p == G4Electron::ElectronDefinition()) {
    auto t2 =
      std::upper_bound(tValuesVec[materialID][p].begin(), tValuesVec[materialID][p].end(), k);
    auto t1 = t2 - 1;

    auto e12 = std::upper_bound(eValuesVect[materialID][p][(*t1)].begin(),
                                eValuesVect[materialID][p][(*t1)].end(), integrDiff);
    auto e11 = e12 - 1;

    auto e22 = std::upper_bound(eValuesVect[materialID][p][(*t2)].begin(),
                                eValuesVect[materialID][p][(*t2)].end(), integrDiff);
    auto e21 = e22 - 1;

    valueT1 = *t1;
    valueT2 = *t2;
    valueE21 = *e21;
    valueE22 = *e22;
    valueE12 = *e12;
    valueE11 = *e11;

    xs11 = diffCrossSectionData[materialID][p][valueT1][valueE11];
    xs12 = diffCrossSectionData[materialID][p][valueT1][valueE12];
    xs21 = diffCrossSectionData[materialID][p][valueT2][valueE21];
    xs22 = diffCrossSectionData[materialID][p][valueT2][valueE22];
  }

  if (xs11 == 0 && xs12 == 0 && xs21 == 0 && xs22 == 0) {
    return (0.);
  }

  theta = QuadInterpolator(valueE11, valueE12, valueE21, valueE22, xs11, xs12, xs21, xs22, valueT1,
                           valueT2, k, integrDiff);

  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1,
                                                 G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = std::exp(d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1,
                                                 G4double xs2)
{
  G4double d1 = xs1;
  G4double d2 = xs2;
  G4double value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1,
                                                 G4double xs2)
{
  G4double a = (std::log10(xs2) - std::log10(xs1)) / (std::log10(e2) - std::log10(e1));
  G4double b = std::log10(xs2) - a * std::log10(e2);
  G4double sigma = a * std::log10(e) + b;
  G4double value = (std::pow(10., sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::QuadInterpolator(G4double e11, G4double e12, G4double e21,
                                                G4double e22, G4double xs11, G4double xs12,
                                                G4double xs21, G4double xs22, G4double t1,
                                                G4double t2, G4double t, G4double e)
{
  // Log-Log
  /*
   G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
   G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
   G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);


   // Lin-Log
   G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
   G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
   G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
 */

  // Lin-Lin
  G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAPTBElasticModel::RandomizeCosTheta(const G4double& k, const std::size_t& materialID)
{
  G4double integrdiff = 0;
  G4double uniformRand = G4UniformRand();
  integrdiff = uniformRand;

  G4double theta = 0.;
  G4double cosTheta = 0.;
  theta = Theta(G4Electron::ElectronDefinition(), k / eV, integrdiff, materialID);

  cosTheta = std::cos(theta * CLHEP::pi / 180);

  return cosTheta;
}
