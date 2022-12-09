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

// Author: Ivana Hrivnacova, 10/08/2022  (ivana@ipno.in2p3.fr)

#include "G4THnToolsManager.hh"

#include "tools/histo/h1d"

using namespace G4Analysis;

#include <fstream>

// Specialization for H1 type


//_____________________________________________________________________________
template <>
tools::histo::h1d* G4THnToolsManager<1, tools::histo::h1d>::CreateToolsHT(
  const G4String& title,
  const std::array<G4HnDimension, 1>& bins,
  const std::array<G4HnDimensionInformation, 1>& hnInfo)
{
  // Apply hn information to bins
  auto newBins(bins[kX]);
  Update(newBins, hnInfo[kX]);

  if (hnInfo[kX].fBinScheme == G4BinScheme::kLinear) {
    return new tools::histo::h1d(title, newBins.fNBins, newBins.fMinValue, newBins.fMaxValue);
  }

  return new tools::histo::h1d(title, newBins.fEdges);
}

//_____________________________________________________________________________
template <>
void G4THnToolsManager<1, tools::histo::h1d>::ConfigureToolsHT(
  tools::histo::h1d* ht,
  const std::array<G4HnDimension, 1>& bins,
  const std::array<G4HnDimensionInformation, 1>& hnInfo)
{
  // Apply hn information to bins
  auto newBins(bins[kX]);
  Update(newBins, hnInfo[kX]);

  if (hnInfo[kX].fBinScheme == G4BinScheme::kLinear) {
    ht->configure(newBins.fNBins, newBins.fMinValue, newBins.fMaxValue);
    return;
  }

  ht->configure(newBins.fEdges);
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<1, tools::histo::h1d>::FillHT(
  tools::histo::h1d* ht, const G4HnInformation& hnInformation, 
  std::array<G4double, 1>& value, G4double weight)
{
  auto xInfo = hnInformation.GetHnDimensionInformation(kX);

  // Apply hn information to value
  Update(value[kX], xInfo);

  // Fill updated value
  ht->fill(value[kX], weight);

  return true;
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<1, tools::histo::h1d>::WriteOnAscii(
  std::ofstream& output)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.

  // Do nothing if no histograms are selected
  if ( ! GetHnManager()->IsAscii() ) return true;

  // Write h1 histograms
  auto id = GetHnManager()->GetFirstId();
  for (const auto& [h1, info] : *GetTHnVector()) {

    if ( ! info->GetAscii() ) {
      // skip writing if activation is enabled and H1 is inactivated
      id++;
      continue;
    }

    Message(kVL3, "write on ascii", "h1d", info->GetName());

    output << "\n  1D histogram " << id++ << ": " << h1->title()
           << "\n \n \t     X \t\t Bin Height" << G4endl;

    for (G4int j=0; j< G4int(h1->axis().bins()); ++j) {
       output << "  " << j << "\t"
              << h1->axis().bin_center(j) << "\t"
              << h1->bin_height(j) << G4endl;
    }
  }

  return output.good();
}
