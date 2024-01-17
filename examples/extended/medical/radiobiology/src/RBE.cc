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
/// \file radiobiology/src/RBE.cc
/// \brief Implementation of the RadioBio::RBE class

#include "RBE.hh"

#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "RBEAccumulable.hh"
#include "RBEMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

// Note that dose is needed in order to fully calculate RBE using
// this class. Therefore, here, the correct dependencies must be
// added.
#include "Dose.hh"
#include "Manager.hh"
#include <G4NistManager.hh>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

#define width 15L

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RBE::RBE() : VRadiobiologicalQuantity()
{
  fPath = "RadioBio";

  // CreateMessenger
  fMessenger = new RBEMessenger(this);

  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RBE::~RBE()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::Initialize()
{
  fLnS.resize(VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber());
  fDoseX.resize(VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber());

  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::Store()
{
  StoreAlphaAndBeta();
  StoreRBE();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::PrintParameters()
{
  G4cout << "*******************************************" << G4endl
         << "******* Parameters of the class RBE *******" << G4endl
         << "*******************************************" << G4endl;
  PrintVirtualParameters();
  G4cout << "** RBE Cell line: " << fActiveCellLine << G4endl;
  G4cout << "** RBE Dose threshold value: " << fDoseCut / gray << " gray" << G4endl;
  G4cout << "** RBE Alpha X value: " << fAlphaX * gray << " 1/gray" << G4endl;
  G4cout << "** RBE Beta X value: " << fBetaX * std::pow(gray, 2.0) << " 1/gray2" << G4endl;
  G4cout << "*******************************************" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @short Split string into parts using a delimiter.
 *
 * @param line a string to be splitted
 * @param delimiter a string to be looked for
 *
 * Usage: Help function for reading CSV
 * Also remove spaces and quotes around.
 * Note: If delimiter is inside a string, the function fails!
 */
std::vector<G4String> split(const G4String& line, const G4String& delimiter)
{
  std::vector<G4String> result;

  size_t current = 0;
  size_t pos = 0;

  while (pos != G4String::npos) {
    pos = line.find(delimiter, current);
    G4String token = line.substr(current, pos - current);

    G4StrUtil::strip(token);
    G4StrUtil::strip(token, '\"');

    result.push_back(token);
    current = pos + delimiter.size();
  }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::LoadLEMTable(G4String path)
{
  std::ifstream in(path);
  if (!in) {
    std::stringstream ss;
    ss << "Cannot open LEM table input file: " << path;
    G4Exception("RBE::LoadData", "WrongTable", FatalException, ss.str().c_str());
  }

  // Start with the first line
  G4String line;
  if (!getline(in, line)) {
    std::stringstream ss;
    ss << "Cannot read header from the LEM table file: " << path;
    G4Exception("RBE::LoadLEMTable", "CannotReadHeader", FatalException, ss.str().c_str());
  }
  std::vector<G4String> header = split(line, ",");

  // Find the order of columns
  std::vector<G4String> columns = {"alpha_x", "beta_x", "D_t",  "specific_energy",
                                   "alpha",   "beta",   "cell", "particle"};
  std::map<G4String, int> columnIndices;
  for (auto columnName : columns) {
    auto pos = find(header.begin(), header.end(), columnName);
    if (pos == header.end()) {
      std::stringstream ss;
      ss << "Column " << columnName << " not present in the LEM table file.";
      G4Exception("RBE::LoadLEMTable", "ColumnNotPresent", FatalException, ss.str().c_str());
    }
    else {
      columnIndices[columnName] = distance(header.begin(), pos);
    }
  }

  // Continue line by line
  while (getline(in, line)) {
    std::vector<G4String> lineParts = split(line, ",");
    G4String cellLine = lineParts[columnIndices["cell"]];
    // G4int element = elements.at(lineParts[columnIndices["particle"]]);
    G4NistManager* man = G4NistManager::Instance();
    G4int element = man->FindOrBuildElement(lineParts[columnIndices["particle"]])->GetZasInt();

    fTablesEnergy[cellLine][element].push_back(
      std::stod(lineParts[columnIndices["specific_energy"]]) * MeV);
    fTablesAlpha[cellLine][element].push_back(stod(lineParts[columnIndices["alpha"]]));
    /* fTablesLet[cellLine][element].push_back(
        stod(lineParts[columnIndices["let"]])
    ); */
    fTablesBeta[cellLine][element].push_back(stod(lineParts[columnIndices["beta"]]));

    fTablesAlphaX[cellLine] = stod(lineParts[columnIndices["alpha_x"]]) / gray;
    fTablesBetaX[cellLine] = stod(lineParts[columnIndices["beta_x"]]) / (gray * gray);
    fTablesDoseCut[cellLine] = stod(lineParts[columnIndices["D_t"]]) * gray;
  }

  // Sort the tables by energy
  // (https://stackoverflow.com/a/12399290/2692780)
  for (auto aPair : fTablesEnergy) {
    for (auto ePair : aPair.second) {
      std::vector<G4double>& tableEnergy = fTablesEnergy[aPair.first][ePair.first];
      std::vector<G4double>& tableAlpha = fTablesAlpha[aPair.first][ePair.first];
      std::vector<G4double>& tableBeta = fTablesBeta[aPair.first][ePair.first];

      std::vector<size_t> idx(tableEnergy.size());
      iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(),
           [&tableEnergy](size_t i1, size_t i2) { return tableEnergy[i1] < tableEnergy[i2]; });

      std::vector<std::vector<G4double>*> tables = {&tableEnergy, &tableAlpha, &tableBeta};
      for (std::vector<G4double>* table : tables) {
        std::vector<G4double> copy = *table;
        for (size_t i = 0; i < copy.size(); ++i) {
          (*table)[i] = copy[idx[i]];
        }
        // G4cout << (*table)[0];
        // reorder(*table, idx);
        // G4cout << (*table)[0] << G4endl;
      }
    }
  }

  if (fVerboseLevel > 0) {
    G4cout << "RBE: read LEM data for the following cell lines and elements [number of points]:"
           << G4endl;
    for (auto aPair : fTablesEnergy) {
      G4cout << "- " << aPair.first << ": ";
      for (auto ePair : aPair.second) {
        G4cout << ePair.first << "[" << ePair.second.size() << "] ";
      }
      G4cout << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::SetCellLine(G4String name)
{
  G4cout << "*************************" << G4endl << "*******SetCellLine*******" << G4endl
         << "*************************" << G4endl;
  if (fTablesEnergy.size() == 0) {
    G4Exception("RBE::SetCellLine", "NoCellLine", FatalException,
                "Cannot select cell line, probably LEM table not loaded yet?");
  }
  if (fTablesEnergy.find(name) == fTablesEnergy.end()) {
    std::stringstream str;
    str << "Cell line " << name << " not found in the LEM table.";
    G4Exception("RBE::SetCellLine", "", FatalException, str.str().c_str());
  }
  else {
    fAlphaX = fTablesAlphaX[name];
    fBetaX = fTablesBetaX[name];
    fDoseCut = fTablesDoseCut[name];

    fActiveTableEnergy = &fTablesEnergy[name];
    fActiveTableAlpha = &fTablesAlpha[name];
    fActiveTableBeta = &fTablesBeta[name];

    fMinZ = 0;
    fMaxZ = 0;
    fMinEnergies.clear();
    fMaxEnergies.clear();

    for (auto energyPair : *fActiveTableEnergy) {
      if (!fMinZ || (energyPair.first < fMinZ)) fMinZ = energyPair.first;
      if (energyPair.first > fMaxZ) fMaxZ = energyPair.first;

      fMinEnergies[energyPair.first] = energyPair.second[0];
      fMaxEnergies[energyPair.first] = energyPair.second[energyPair.second.size() - 1];
    }
  }

  fActiveCellLine = name;

  if (fVerboseLevel > 0) {
    G4cout << "RBE: cell line " << name << " selected." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::tuple<G4double, G4double> RBE::GetHitAlphaAndBeta(G4double E, G4int Z)
{
  if (!fActiveTableEnergy) {
    G4Exception("RBE::GetHitAlphaAndBeta", "NoCellLine", FatalException,
                "No cell line selected. Please, do it using the /rbe/cellLine command.");
  }

  // Check we are in the tables
  if ((Z < fMinZ) || (Z > fMaxZ)) {
    if (fVerboseLevel > 2) {
      std::stringstream str;
      str << "Alpha & beta calculation supported only for ";
      str << fMinZ << " <= Z <= " << fMaxZ << ", but " << Z << " requested.";
      G4Exception("RBE::GetHitAlphaAndBeta", "", JustWarning, str.str().c_str());
    }
    return std::make_tuple<G4double, G4double>(0.0, 0.0);  // out of table!
  }
  if ((E < fMinEnergies[Z]) || (E >= fMaxEnergies[Z])) {
    if (fVerboseLevel > 2) {
      G4cout << "RBE hit: Z=" << Z << ", E=" << E << " => out of LEM table";
      if (fVerboseLevel > 3) {
        G4cout << " (" << fMinEnergies[Z] << " to " << fMaxEnergies[Z] << " MeV)";
      }
      G4cout << G4endl;
    }
    return std::make_tuple<G4double, G4double>(0.0, 0.0);  // out of table!
  }

  std::vector<G4double>& vecEnergy = (*fActiveTableEnergy)[Z];
  std::vector<G4double>& vecAlpha = (*fActiveTableAlpha)[Z];
  std::vector<G4double>& vecBeta = (*fActiveTableBeta)[Z];

  // Find the row in energy tables
  const auto eLarger = upper_bound(begin(vecEnergy), end(vecEnergy), E);
  const G4int lower = distance(begin(vecEnergy), eLarger) - 1;
  const G4int upper = lower + 1;

  // Interpolation
  const G4double energyLower = vecEnergy[lower];
  const G4double energyUpper = vecEnergy[upper];
  const G4double energyFraction = (E - energyLower) / (energyUpper - energyLower);

  // linear interpolation along E
  const G4double alpha =
    ((1 - energyFraction) * vecAlpha[lower] + energyFraction * vecAlpha[upper]);
  const G4double beta = ((1 - energyFraction) * vecBeta[lower] + energyFraction * vecBeta[upper]);
  if (fVerboseLevel > 2) {
    G4cout << "RBE hit: Z=" << Z << ", E=" << E << " => alpha=" << alpha << ", beta=" << beta
           << G4endl;
  }

  return std::make_tuple(alpha, beta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Zaider & Rossi alpha & Beta mean
void RBE::ComputeAlphaAndBeta()
{
    // Skip RBE computation if calculation not enabled.
    if (!fCalculationEnabled) {
        if (fVerboseLevel > 0) {
            G4cout << "RBE::ComputeAlphaAndBeta() called but skipped as calculation not enabled"
                   << G4endl;
        }
        return;
    }

    if (fVerboseLevel > 0) {
        G4cout << "RBE: Computing alpha and beta..." << G4endl;
    }

    // Re-initialize the number of voxels
    fAlpha.resize(fAlphaNumerator.size()); // Initialize with the same number of elements
    fBeta.resize(fBetaNumerator.size());   // Initialize with the same number of elements

    for (size_t ii = 0; ii < fDenominator.size(); ii++) {
        if (fDenominator[ii] > 0) {
            fAlpha[ii] = fAlphaNumerator[ii] / (fDenominator[ii] * gray);
            fBeta[ii] = std::pow(fBetaNumerator[ii] / (fDenominator[ii] * gray), 2.0);
        } 
        
        else {
            fAlpha[ii] = 0.;
            fBeta[ii] = 0.;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::ComputeRBE()
{
  // Skip RBE computation if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "RBE::ComputeRBE() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fVerboseLevel > 0) {
    G4cout << "RBE: Computing survival and RBE..." << G4endl;
  }
  G4double smax = fAlphaX + 2 * fBetaX * fDoseCut;

  for (G4int i = 0; i < VoxelizedSensitiveDetector::GetInstance()->GetTotalVoxelNumber(); i++) {
    if (std::isnan(fAlpha[i]) || std::isnan(fBeta[i])) {
      fLnS[i] = 0.0;
      fDoseX[i] = 0.0;
    }
    else if (fDose[i] <= fDoseCut) {
      fLnS[i] = -(fAlpha[i] * fDose[i]) - (fBeta[i] * (std::pow(fDose[i], 2.0)));
      fDoseX[i] =
        std::sqrt((-fLnS[i] / fBetaX) + std::pow((fAlphaX / (2 * fBetaX)), 2.0)) - (fAlphaX / (2 * fBetaX));
    }
    else {
      G4double ln_Scut = -(fAlpha[i] * fDoseCut) - (fBeta[i] * (std::pow((fDoseCut), 2.0)));
      fLnS[i] = ln_Scut - ((fDose[i] - fDoseCut) * smax);

      fDoseX[i] = ((-fLnS[i] + ln_Scut) / smax) + fDoseCut;
    }
  }
  fRBE = fDoseX / fDose;
  fSurvival = std::exp(fLnS);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::Compute()
{
  // Skip RBE computation if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "RBE::Compute() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  if (fCalculated == true) return;

  GetDose();

  ComputeAlphaAndBeta();
  ComputeRBE();

  // If this method reached the bottom, set calculated to true
  fCalculated = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::GetDose()
{
  // Get the pointer to dose. If it does not exist, launch an exception
  const Dose* dose = dynamic_cast<const Dose*>(Manager::GetInstance()->GetQuantity("Dose"));
  if (dose == nullptr)
    G4Exception("RBE::Compute", "RBEMissingDose", FatalException,
                "Trying to compute RBE without knowing the dose. Please add a valid dose or "
                "disable RBE calculation");

  // Check whether dose has been calculated.
  // If not, give a warning
  if (!dose->HasBeenCalculated())
    G4Exception("RBE::Compute", "RBEDoseNotCalculated", JustWarning,
                "Dose has not been calculated yet while computing RBE, that will be wrong."
                " Please, first calculate dose");
  // Copy the proper vector from Dose class to RBE class
  fDose = dose->GetDose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::SetDenominator(const RBE::array_type denom)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Setting denominator..." << G4endl;
  }
  fDenominator = denom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::AddDenominator(const RBE::array_type denom)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Adding denominator...";
  }
  if (fDenominator.size()) {
    fDenominator += denom;
  }
  else {
    if (fVerboseLevel > 1) {
      G4cout << " (created empty array)";
    }
    fDenominator = denom;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::SetAlphaNumerator(const RBE::array_type alpha)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Setting alpha numerator..." << G4endl;
  }
  fAlphaNumerator = alpha;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::SetBetaNumerator(const RBE::array_type beta)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Setting beta numerator..." << G4endl;
  }
  fBetaNumerator = beta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::AddAlphaNumerator(const RBE::array_type alpha)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Adding alpha numerator...";
  }
  if (fAlphaNumerator.size()) {
    fAlphaNumerator += alpha;
  }
  else {
    if (fVerboseLevel > 1) {
      G4cout << " (created empty array)";
    }
    fAlphaNumerator = alpha;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::AddBetaNumerator(const RBE::array_type beta)
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Adding beta numerator...";
  }
  if (fBetaNumerator.size()) {
    fBetaNumerator += beta;
  }
  else {
    if (fVerboseLevel > 1) {
      G4cout << " (created empty array)";
    }
    fBetaNumerator = beta;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::StoreAlphaAndBeta()
{
  // Skip RBE storing if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "RBE::StoreAlphaAndBeta() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  G4String AlphaBetaPath = fPath + "_AlphaAndBeta.out";
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Writing alpha and beta..." << G4endl;
  }
  std::ofstream ofs(AlphaBetaPath);

  Compute();

  if (ofs.is_open()) {
    ofs << std::left << std::setw(width) << "i" << std::setw(width) << "j" << std::setw(width)
        << "k" << std::setw(width) << "alpha" << std::setw(width) << "beta " << G4endl;

    auto voxSensDet = VoxelizedSensitiveDetector::GetInstance();

    // Alpha and beta are written only if valid. If value is -nan, 0 is written
    // on the text file
    for (G4int i = 0; i < voxSensDet->GetVoxelNumberAlongX(); i++)
      for (G4int j = 0; j < voxSensDet->GetVoxelNumberAlongY(); j++)
        for (G4int k = 0; k < voxSensDet->GetVoxelNumberAlongZ(); k++) {
          G4int v = voxSensDet->GetThisVoxelNumber(i, j, k);

          ofs << std::left << std::setw(width) << i << std::setw(width) << j << std::setw(width)
              << k << std::setw(width) << (std::isnan(fAlpha[v]) ? 0 : (fAlpha[v] * gray))
              << std::setw(width) << (std::isnan(fBeta[v]) ? 0 : (fBeta[v] * std::pow(gray, 2.0)))
              << G4endl;
        }
  }
  if (fVerboseLevel > 0) {
    G4cout << "RBE: Alpha and beta written to " << AlphaBetaPath << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::StoreRBE()
{
  // Skip RBE storing if calculation not enabled.
  if (!fCalculationEnabled) {
    if (fVerboseLevel > 0) {
      G4cout << "RBE::StoreRBE() called but skipped as calculation not enabled" << G4endl;
    }
    return;
  }

  G4String RBEPath = fPath + "_RBE.out";
  if (fSaved == true)
    G4Exception("RBE::StoreRBE", "RBEOverwrite", JustWarning,
                "Overwriting RBE file. For multiple runs, change filename.");
  std::ofstream ofs(RBEPath);

  Compute();

  if (ofs.is_open()) {
    ofs << std::left << std::setw(width) << "i" << std::setw(width) << "j" << std::setw(width)
        << "k" << std::setw(width) << "Dose(Gy)" << std::setw(width) << "ln(S) " << std::setw(width)
        << "Survival" << std::setw(width) << "DoseB(Gy)" << std::setw(width) << "RBE" << G4endl;

    auto voxSensDet = VoxelizedSensitiveDetector::GetInstance();

    for (G4int i = 0; i < voxSensDet->GetVoxelNumberAlongX(); i++)
      for (G4int j = 0; j < voxSensDet->GetVoxelNumberAlongY(); j++)
        for (G4int k = 0; k < voxSensDet->GetVoxelNumberAlongZ(); k++) {
          G4int v = voxSensDet->GetThisVoxelNumber(i, j, k);

          ofs << std::left << std::setw(width) << i << std::setw(width) << j << std::setw(width)
              << k << std::setw(width) << (fDose[v] / gray) << std::setw(width) << fLnS[v]
              << std::setw(width) << fSurvival[v] << std::setw(width) << fDoseX[v] / gray
              << std::setw(width) << (std::isnan(fRBE[v]) ? 0. : fRBE[v]) << G4endl;
        }
  }
  if (fVerboseLevel > 0) {
    G4cout << "RBE: RBE written to " << RBEPath << G4endl;
  }

  fSaved = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::Reset()
{
  if (fVerboseLevel > 1) {
    G4cout << "RBE: Reset(): ";
  }
  fAlphaNumerator = 0.0;
  fBetaNumerator = 0.0;
  fDenominator = 0.0;
  fDose = 0.0;
  fCalculated = false;
  if (fVerboseLevel > 1) {
    G4cout << fAlphaNumerator.size() << " points." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBE::AddFromAccumulable(G4VAccumulable* GenAcc)
{
  RBEAccumulable* acc = (RBEAccumulable*)GenAcc;
  AddAlphaNumerator(acc->GetAlphaNumerator());
  AddBetaNumerator(acc->GetBetaNumerator());
  AddDenominator(acc->GetDenominator());

  fCalculated = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
