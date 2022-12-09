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

#include "tools/histo/p1d"

using namespace G4Analysis;

#include <fstream>

// Specialization for P1 type


//_____________________________________________________________________________
template <>
tools::histo::p1d* G4THnToolsManager<2, tools::histo::p1d>::CreateToolsHT(
  const G4String& title,
  const std::array<G4HnDimension, 2>& bins,
  const std::array<G4HnDimensionInformation, 2>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  UpdateValues(newYBins, hnInfo[kY]);

  if (hnInfo[kX].fBinScheme == G4BinScheme::kLinear) {
    if ( newYBins.fMinValue == 0. && newYBins.fMaxValue == 0.) {
      return new tools::histo::p1d(
        title, newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue);
    }
    return new tools::histo::p1d(
      title, newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
      newYBins.fMinValue, newYBins.fMaxValue);
  }

  if ( newYBins.fMinValue == 0. && newYBins.fMaxValue == 0.) {
    return new tools::histo::p1d(title, newXBins.fEdges);
  }

  return new tools::histo::p1d(title, newXBins.fEdges,
    newYBins.fMinValue, newYBins.fMaxValue);
}

//_____________________________________________________________________________
template <>
void G4THnToolsManager<2, tools::histo::p1d>::ConfigureToolsHT(
  tools::histo::p1d* ht,
  const std::array<G4HnDimension, 2>& bins,
  const std::array<G4HnDimensionInformation, 2>& hnInfo)
{
  // Apply hn information to bins
  auto newXBins(bins[kX]);
  Update(newXBins, hnInfo[kX]);
  auto newYBins(bins[kY]);
  UpdateValues(newYBins, hnInfo[kY]);

  if (hnInfo[kX].fBinScheme == G4BinScheme::kLinear) {
    if ( newYBins.fMinValue == 0. && newYBins.fMaxValue == 0.) {
      ht->configure(
            newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue);
      return;
    }
    ht->configure(
          newXBins.fNBins, newXBins.fMinValue, newXBins.fMaxValue,
          newYBins.fMinValue, newYBins.fMaxValue);
    return;
  }

  if ( newYBins.fMinValue == 0. && newYBins.fMaxValue == 0.) {
    ht->configure(newXBins.fEdges);
  return;
  }
  ht->configure(newXBins.fEdges, newYBins.fMinValue, newYBins.fMaxValue);
  return;
}

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<2, tools::histo::p1d>::FillHT(
  tools::histo::p1d* ht, const G4HnInformation& hnInformation, 
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
G4bool G4THnToolsManager<2, tools::histo::p1d>::WriteOnAscii(
  std::ofstream& output)
{
// Write selected objects on ASCII file

  // Do nothing if no histograms are selected
  if ( ! GetHnManager()->IsAscii() ) return true;

  // Write p1 histograms
  auto id = GetHnManager()->GetFirstId();
  for (const auto& [p1, info] : *GetTHnVector()) {

    if ( ! info->GetAscii() ) {
      // skip writing if activation is enabled and H1 is inactivated
      id++;
      continue;
    }

    Message(kVL3, "write on ascii", "p1d", info->GetName());

    output << "\n  1D profile " << id++ << ": " << p1->title()
           << "\n \n \t \t     X \t\t MeanY" << G4endl;

    for (G4int j=0; j< G4int(p1->axis().bins()); ++j) {
       auto sw = p1->bin_Sw(j);
       auto svw = p1->bin_Svw(j);
       auto mean = ( sw != 0. ) ?  (svw / sw) : 0.;
       output << "  " << j << "\t"
              << p1->axis().bin_center(j) << "\t"
              << mean << G4endl;
    }
  }

  return output.good();
}
