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

// Author: Ivana Hrivnacova, 22/08/2013  (ivana@ipno.in2p3.fr)

#include "G4BinScheme.hh"
#include "G4AnalysisUtilities.hh"

namespace G4Analysis
{

//_____________________________________________________________________________
G4BinScheme GetBinScheme(const G4String& binSchemeName)
{
  if (binSchemeName == "linear") return G4BinScheme::kLinear;
  if (binSchemeName == "log")    return G4BinScheme::kLog;
  if (binSchemeName == "user")   return G4BinScheme::kUser;

  // No other name is supported
  Warn("\"" + binSchemeName + "\" binning scheme is not supported.\n"
       "Linear binning will be applied.", kNamespaceName, "GetBinScheme");

  return G4BinScheme::kLinear;
}

//_____________________________________________________________________________
void ComputeEdges(G4int nbins, G4double xmin, G4double xmax,
                  G4double unit, G4Fcn fcn, G4BinScheme binScheme,
                  std::vector<G4double>& edges)
{
// Compute edges from parameters

  if ( binScheme == G4BinScheme::kUser ) {
    // This call should never happen for user bin scheme
    Warn("There is no need to compute edges for G4BinScheme::kUser\n"
         "Call is ignored.", kNamespaceName, "GetBinScheme");
    return;
  }

  if (unit == 0.) {
    // Should never happen
    Warn("Illegal unit value (0), 1. will be used instead",
      kNamespaceName, "ComputeEdges");
    unit = 1.;
  }

  if (nbins == 0) {
    // Should never happen
    Warn("Illegal number of nbins value (0), call will be ignored",
      kNamespaceName, "ComputeEdges");
    return;
  }

  // Apply units
  auto xumin = xmin/unit;
  auto xumax = xmax/unit;

  if ( binScheme == G4BinScheme::kLinear ) {
    auto dx = (fcn(xumax) - fcn(xumin) ) / nbins;
    auto binValue = fcn(xumin);
    while ( G4int(edges.size()) <= nbins ) {  // Loop checking, 23.06.2015, I. Hrivnacova
      edges.push_back(binValue);
      binValue += dx;
    }
    return;
  }

  if ( binScheme == G4BinScheme::kLog ) {
    // do not apply fcn
    auto dlog = (std::log10(xumax) - std::log10(xumin))/ nbins;
    auto dx = std::pow(10, dlog);
    auto binValue = xumin;
    while ( G4int(edges.size()) <= nbins ) { // Loop checking, 23.06.2015, I. Hrivnacova
      edges.push_back(binValue);
      binValue *= dx;
    }
    return;
  }
}

//_____________________________________________________________________________
void ComputeEdges(const std::vector<G4double>& edges,
                  G4double unit, G4Fcn fcn,
                  std::vector<G4double>& newBins)
{
// Apply function & unit to defined edges

  if (unit == 0.) {
    // Should never happen
    Warn("Illegal unit value (0), 1. will be used instead",
      kNamespaceName, "ComputeEdges");
    unit = 1.;
  }

  for (auto element : edges) {
    newBins.push_back(fcn(element/unit));
  }
}

}
