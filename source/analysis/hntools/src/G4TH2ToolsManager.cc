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

#include "tools/histo/h2d"

using namespace G4Analysis;

#include <fstream>

// Specialization for H2 type


//_____________________________________________________________________________
template <>
tools::histo::h2d* G4THnToolsManager<2, tools::histo::h2d>::CreateToolsHT(
  const G4String& title,
  const std::array<G4HnDimension, 2>& bins,
  const std::array<G4HnDimensionInformation, 2>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  Update(newYBins, hnInfo[kY]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear)) {
    return new tools::histo::h2d(title,
                                 newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
                                 newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue);
  }

  return new tools::histo::h2d(title, newXBins.fEdges, newYBins.fEdges);
}

//_____________________________________________________________________________
template <>
void G4THnToolsManager<2, tools::histo::h2d>::ConfigureToolsHT(
  tools::histo::h2d* ht,
  const std::array<G4HnDimension, 2>& bins,
  const std::array<G4HnDimensionInformation, 2>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  Update(newYBins, hnInfo[kY]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear)) {
    ht->configure(
          newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
          newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue);
    return;
  }

  ht->configure(newXBins.fEdges, newYBins.fEdges);
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<2, tools::histo::h2d>::FillHT(
  tools::histo::h2d* ht, const G4HnInformation& hnInformation, 
  std::array<G4double, 2>& value, G4double weight)
{
  auto xInfo = hnInformation.GetHnDimensionInformation(kX);
  auto yInfo = hnInformation.GetHnDimensionInformation(kY);

  // Apply hn information to value
  Update(value[kX], xInfo);
  Update(value[kY], yInfo);

  // Fill updated value
  ht->fill(value[kX], value[kY], weight);

  return true;
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<2, tools::histo::h2d>::WriteOnAscii(
  std::ofstream& output)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.

  // Do nothing if no histograms are selected
  if ( ! GetHnManager()->IsAscii() ) return true;

  // Write h2 histograms
  auto id = GetHnManager()->GetFirstId();
  for (const auto& [h2, info] : *GetTHnVector()) {

    if ( ! info->GetAscii() ) {
      // skip writing if activation is enabled and H1 is inactivated
      id++;
      continue;
    }

    Message(kVL3, "write on ascii", "h2d", info->GetName());

    output << "\n  2D histogram " << id++ << ": " << h2->title()
           << "\n \n \t \t     X \t\t     Y \t\t Bin Height" << G4endl;

    for (G4int j=0; j< G4int(h2->axis_x().bins()); ++j) {
      for (G4int k=0; k< G4int(h2->axis_y().bins()); ++k) {
        output << "  " << j << "\t" << k << "\t"
               << h2->axis_x().bin_center(j) << "\t"
               << h2->axis_y().bin_center(k) << "\t"
               << h2->bin_height(j, k) << G4endl;
      }
    }
  }

  return output.good();
}
