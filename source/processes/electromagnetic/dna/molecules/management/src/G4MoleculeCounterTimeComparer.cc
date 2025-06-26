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
// Author: Christian Velten (2025)

#include "G4MoleculeCounterTimeComparer.hh"

#include "G4Exception.hh"

G4MoleculeCounterTimeComparer::G4MoleculeCounterTimeComparer()
  : fType(TimeComparerType::FixedPrecision)
{}

G4MoleculeCounterTimeComparer::G4MoleculeCounterTimeComparer(
  const G4MoleculeCounterTimeComparer& other)
{
  *this = other;
}

G4MoleculeCounterTimeComparer&
G4MoleculeCounterTimeComparer::operator=(const G4MoleculeCounterTimeComparer& other)
{
  switch (other.fType) {
    case TimeComparerType::VariablePrecision:
      SetVariablePrecision(other.fVariablePrecision);
      break;
    case TimeComparerType::FixedPrecision:
    default:
      SetFixedPrecision(other.fPrecision);
      break;
  }
  return *this;
}

void G4MoleculeCounterTimeComparer::SetFixedPrecision(G4double precision)
{
  fType = TimeComparerType::FixedPrecision;
  fPrecision = precision;
}

void G4MoleculeCounterTimeComparer::SetVariablePrecision(
  const std::vector<G4double>& globalTimeHighEdges, const std::vector<G4double>& resolutions)
{
  fType = TimeComparerType::VariablePrecision;
  fVariablePrecision.clear();
  for (auto key = globalTimeHighEdges.cbegin(), val = resolutions.cbegin();
       key != globalTimeHighEdges.cend() && val != resolutions.cend(); ++key, ++val)
  {
    fVariablePrecision.emplace(*key, *val);
  }

  if (fVariablePrecision.empty()) {
    G4Exception("G4MoleculeCounterTimeComparer::SetVariablePrecision()",
                "G4MoleculeCounterTimeComparer", FatalException, "Precision map cannot be empty!");
  }
}

void G4MoleculeCounterTimeComparer::SetVariablePrecision(
  const std::map<G4double, G4double>& resolutionMap)
{
  fType = TimeComparerType::VariablePrecision;
  fVariablePrecision = resolutionMap;

  if (fVariablePrecision.empty()) {
    G4Exception("G4MoleculeCounterTimeComparer::SetVariablePrecision()",
                "G4MoleculeCounterTimeComparer", FatalException, "Precision map cannot be empty!");
  }
}

G4double G4MoleculeCounterTimeComparer::GetPrecisionAtTime(G4double time) const
{
  if (fType == G4MoleculeCounterTimeComparer::FixedPrecision)
    return fPrecision;
  else
    return fVariablePrecision.upper_bound(time)->second;
}

G4bool G4MoleculeCounterTimeComparer::operator()(const G4double& a, const G4double& b) const
{
  if (fType == G4MoleculeCounterTimeComparer::FixedPrecision) {
    if (std::fabs(a - b) < fPrecision) {
      return false;
    }
  }
  else if (fType == G4MoleculeCounterTimeComparer::VariablePrecision) {
    const auto it_a = fVariablePrecision.upper_bound(a);
    const auto it_b = fVariablePrecision.upper_bound(b);
    auto precision = fVariablePrecision.crbegin()->second;
    if (it_a != fVariablePrecision.cend() && it_b != fVariablePrecision.cend())
      precision = std::min(it_a->second, it_b->second);
    if (std::fabs(a - b) < precision) {
      return false;
    }
  }
  else {
    G4Exception("G4MoleculeCounterTimeComparer::operator()", "G4MoleculeCounterTimeComparer",
                FatalException, "Unknown comparison type");
  }
  return a < b;
}

// ---------------------------------------------------------------------

G4MoleculeCounterTimeComparer
G4MoleculeCounterTimeComparer::CreateWithFixedPrecision(G4double precision)
{
  auto obj = G4MoleculeCounterTimeComparer();
  obj.SetFixedPrecision(precision);
  return obj;
}

G4MoleculeCounterTimeComparer
G4MoleculeCounterTimeComparer::CreateWithVariablePrecision(const std::map<G4double, G4double>& map)
{
  auto obj = G4MoleculeCounterTimeComparer();
  obj.SetVariablePrecision(map);
  return obj;
}
