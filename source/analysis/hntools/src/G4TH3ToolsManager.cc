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

#include "tools/histo/h3d"

using namespace G4Analysis;

#include <fstream>

// Specialization for H3 type


//_____________________________________________________________________________
template <>
tools::histo::h3d* G4THnToolsManager<3, tools::histo::h3d>::CreateToolsHT(
  const G4String& title,
  const std::array<G4HnDimension, 3>& bins,
  const std::array<G4HnDimensionInformation, 3>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  Update(newYBins, hnInfo[kY]);
  auto newZBins(bins[kZ]);
  Update(newZBins, hnInfo[kZ]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kZ].fBinScheme == G4BinScheme::kLinear)) {
    return new tools::histo::h3d(title,
                 newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
                 newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue,
                 newZBins.fNBins, newZBins.fMinValue, newZBins.fMaxValue);
  }

  return new tools::histo::h3d(title,
               newXBins.fEdges, newYBins.fEdges, newZBins.fEdges);
}

//_____________________________________________________________________________
template <>
void G4THnToolsManager<3, tools::histo::h3d>::ConfigureToolsHT(
  tools::histo::h3d* ht,
  const std::array<G4HnDimension, 3>& bins,
  const std::array<G4HnDimensionInformation, 3>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  Update(newYBins, hnInfo[kY]);
  auto newZBins(bins[kZ]);
  Update(newZBins, hnInfo[kZ]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kZ].fBinScheme == G4BinScheme::kLinear)) {
    ht->configure(
          newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
          newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue,
          newZBins.fNBins, newZBins.fMinValue, newZBins.fMaxValue);
    return;
  }

  ht->configure(newXBins.fEdges, newYBins.fEdges, newZBins.fEdges);
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<3, tools::histo::h3d>::FillHT(
  tools::histo::h3d* ht, const G4HnInformation& hnInformation, 
  std::array<G4double, 3>& value, G4double weight)
{
  auto xInfo = hnInformation.GetHnDimensionInformation(kX);
  auto yInfo = hnInformation.GetHnDimensionInformation(kY);
  auto zInfo = hnInformation.GetHnDimensionInformation(kZ);

  // Apply hn information to value
  Update(value[kX], xInfo);
  Update(value[kY], yInfo);
  Update(value[kZ], zInfo);

  // Fill updated value
  ht->fill(value[kX], value[kY], value[kZ], weight);

  return true;
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<3, tools::histo::h3d>::WriteOnAscii(
  std::ofstream& output)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.

  // Do nothing if no histograms are selected
  if ( ! GetHnManager()->IsAscii() ) return true;

  // Write h3 histograms
  auto id = GetHnManager()->GetFirstId();
  for (const auto& [h3, info] : *GetTHnVector()) {

    if ( ! info->GetAscii() ) {
      // skip writing if activation is enabled and H1 is inactivated
      id++;
      continue;
    }

    Message(kVL3, "write on ascii", "h3d", info->GetName());

    output << "\n  3D histogram " << id++ << ": " << h3->title()
           << "\n \n \t \t \t     X \t\t     Y \t\t     Z \t\t Bin Height" << G4endl;

    for (G4int j=0; j< G4int(h3->axis_x().bins()); ++j) {
      for (G4int k=0; k< G4int(h3->axis_y().bins()); ++k) {
        for (G4int l=0; l< G4int(h3->axis_y().bins()); ++l) {
          output << "  " << j << "\t" << k << "\t" << l << "\t"
                 << h3->axis_x().bin_center(j) << "\t"
                 << h3->axis_y().bin_center(k) << "\t"
                 << h3->axis_y().bin_center(l) << "\t"
                 << h3->bin_height(j, k, l) << G4endl;
        }
      }
    }
  }

  return output.good();
}
