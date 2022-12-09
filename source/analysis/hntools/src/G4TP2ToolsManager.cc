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

#include "tools/histo/p2d"

using namespace G4Analysis;

#include <fstream>

// Specialization for P2 type


//_____________________________________________________________________________
template <>
tools::histo::p2d* G4THnToolsManager<3, tools::histo::p2d>::CreateToolsHT(
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
  UpdateValues(newZBins, hnInfo[kZ]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear)) {
    if ( newZBins.fMinValue == 0. && newZBins.fMaxValue == 0.) {
      return new tools::histo::p2d(
        title, newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
        newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue);
    }
    return new tools::histo::p2d(
      title, newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
      newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue,
      newZBins.fMinValue, newZBins.fMaxValue);
  }

  if ( newZBins.fMinValue == 0. && newZBins.fMaxValue == 0.) {
    return new tools::histo::p2d(title, newXBins.fEdges, newYBins.fEdges);
  }
  return new tools::histo::p2d(title, newXBins.fEdges, newYBins.fEdges,
    newZBins.fMinValue, newZBins.fMaxValue);
}

//_____________________________________________________________________________
template <>
void G4THnToolsManager<3, tools::histo::p2d>::ConfigureToolsHT(
  tools::histo::p2d* ht,
  const std::array<G4HnDimension, 3>& bins,
  const std::array<G4HnDimensionInformation, 3>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  Update(newYBins, hnInfo[kY]);
  auto newZBins(bins[kZ]);
  UpdateValues(newZBins, hnInfo[kZ]);

  if ((hnInfo[kX].fBinScheme == G4BinScheme::kLinear) &&
      (hnInfo[kY].fBinScheme == G4BinScheme::kLinear)) {
    if ( newZBins.fMinValue == 0. && newZBins.fMaxValue == 0.) {
      ht->configure(
            newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
            newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue);
      return;
    }
    ht->configure(
          newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
          newYBins.fNBins, newYBins.fMinValue, newYBins.fMaxValue,
          newZBins.fMinValue, newZBins.fMaxValue);
    return;
  }

  if ( newZBins.fMinValue == 0. && newZBins.fMaxValue == 0.) {
    ht->configure(newXBins.fEdges, newYBins.fEdges);
    return;
  }

  ht->configure(newXBins.fEdges, newYBins.fEdges,
        newZBins.fMinValue, newZBins.fMaxValue);
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<3, tools::histo::p2d>::FillHT(
  tools::histo::p2d* ht, const G4HnInformation& hnInformation, 
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
G4bool G4THnToolsManager<3, tools::histo::p2d>::WriteOnAscii(
  std::ofstream& output)
{
// Write selected objects on ASCII file

  // Do nothing if no histograms are selected
  if ( ! GetHnManager()->IsAscii() ) return true;

  // Write p2 histograms
  auto id = GetHnManager()->GetFirstId();
  for (const auto& [p2, info] : *GetTHnVector()) {

    if ( ! info->GetAscii() ) {
      // skip writing if activation is enabled and H1 is inactivated
      id++;
      continue;
    }

    Message(kVL3, "write on ascii", "p2d", info->GetName());

    output << "\n  2D profile " << id++ << ": " << p2->title()
           << "\n \n \t \t     X \t\t     Y \t\t    MeanZ" << G4endl;

    for (G4int j=0; j< G4int(p2->axis_x().bins()); ++j) {
      for (G4int k=0; k< G4int(p2->axis_y().bins()); ++k) {
        auto sw = p2->bin_Sw(j, k);
        auto svw = p2->bin_Svw(j, k);
        auto mean = ( sw != 0. ) ?  (svw / sw) : 0.;
        output << "  " << j << "\t" << k << "\t"
               << p2->axis_x().bin_center(j) << "\t"
               << p2->axis_y().bin_center(k) << "\t"
               << mean << G4endl;
      }
    }
  }

  return output.good();
}
