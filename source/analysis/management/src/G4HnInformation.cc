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

// Author: Ivana Hrivnacova,  10/08/2022  (ivana@ipno.in2p3.fr)

#include "G4HnInformation.hh"
#include "G4AnalysisUtilities.hh"


//_____________________________________________________________________________
void G4HnDimension::Print() const
{
  G4cout
    << "NBins: " << fNBins << " minValue: " << fMinValue << " maxValue: "
    << fMaxValue << ";" << " edges: ";
  for ( auto value : fEdges ) {
    G4cout << value << ", ";
  }
  G4cout << G4endl;
}

//_____________________________________________________________________________
void G4HnDimensionInformation::Print() const
{
  G4cout
    << "Unit name: " << fUnitName << " Fcn Name: " << fFcnName << " BinSchemeName: "
    << fBinSchemeName << " Unit: " << fUnit << " BinScheme: " << static_cast<int>(fBinScheme)
    << G4endl;
}

namespace G4Analysis
{

//_____________________________________________________________________________
void Update(G4double& value, const G4HnDimensionInformation& hnInfo)
{
// Apply hnInfo to a value

  auto unit = hnInfo.fUnit;
  auto fcn = hnInfo.fFcn;

  if (unit == 0.) {
    // Should never happen
    Warn("Illegal unit value (0), 1. will be used instead",
      kNamespaceName, "UpdateBins");
    unit = 1.;
  }
  value = fcn(value/unit);
}

//_____________________________________________________________________________
void UpdateValues(
  G4HnDimension& bins, const G4HnDimensionInformation& hnInfo)
{
// Apply hnInfo to bins min and max value

  auto unit = hnInfo.fUnit;
  auto fcn = hnInfo.fFcn;

  if (unit == 0.) {
    // Should never happen
    Warn("Illegal unit value (0), 1. will be used instead",
      kNamespaceName, "UpdateBins");
    unit = 1.;
  }
  // Update min/max values
  bins.fMinValue = fcn(bins.fMinValue/unit);
  bins.fMaxValue = fcn(bins.fMaxValue/unit);
}

//_____________________________________________________________________________
void Update(
  G4HnDimension& bins, const G4HnDimensionInformation& hnInfo)
{
// Apply hnInfo to bins, compute edges

  auto unit = hnInfo.fUnit;
  auto fcn = hnInfo.fFcn;
  auto binScheme = hnInfo.fBinScheme;

  if (binScheme == G4BinScheme::kLinear) {
    // Compute edges, as they may be needed in the context of 2D or 3D histograms
    // with log binning in other dimension
    G4Analysis::ComputeEdges(
      bins.fNBins, bins.fMinValue, bins.fMaxValue, unit, fcn, binScheme, bins.fEdges);

    // Update min/max Values
    UpdateValues(bins, hnInfo);

    return;
  }

  if (binScheme == G4BinScheme::kLog) {
    // Logarithmic bin scheme
    // compute edges from parameters
    G4Analysis::ComputeEdges(
      bins.fNBins, bins.fMinValue, bins.fMaxValue, unit, fcn, binScheme, bins.fEdges);
  }

  if (binScheme == G4BinScheme::kUser) {
    G4Analysis::ComputeEdges(bins.fEdges, unit, fcn, bins.fEdges);
  }
}

//_____________________________________________________________________________
void UpdateTitle(G4String& title, const G4HnDimensionInformation& hnInfo)
{
  if ( hnInfo.fFcnName != "none" )  { title += " "; title +=  hnInfo.fFcnName; title += "("; }
  if ( hnInfo.fUnitName != "none" ) { title += " ["; title += hnInfo.fUnitName; title += "]";}
  if ( hnInfo.fFcnName != "none" )  { title += ")"; }
}

//_____________________________________________________________________________
G4bool CheckMinMax(G4double minValue, G4double maxValue)
{
  auto result = true;

  // Do not check default values
  if ( minValue == 0. && maxValue == 0. ) return result;

  if ( maxValue <= minValue ) {
    Warn("Illegal value of (minValue >= maxMaxValue)", kNamespaceName, "CheckMinMax");
    result = false;
  }

  return result;
}

//_____________________________________________________________________________
G4bool CheckDimension(unsigned int idim,
         const G4HnDimension& dimension, const G4HnDimensionInformation& info)
{
  auto result = true;
  G4String xyz {"xyz"};

  // Check nbins
  if ( (dimension.fNBins <= 0) && (info.fBinScheme != G4BinScheme::kUser) ) {
    Warn("Illegal value of number of " + xyz.substr(idim,1) + " bins: nbins <= 0.",
      kNamespaceName, "CheckDimension");
    result = false;
  }

  // Check edges
  if ( dimension.fEdges.empty() && (info.fBinScheme == G4BinScheme::kUser) ) {
    Warn("Illegal value of number of " + xyz.substr(idim,1) + " edges vector size",
      kNamespaceName, "CheckDimension");
    result = false;
  }

  // Check min/max
  if ( dimension.fMaxValue <= dimension.fMinValue ) {
    Warn("Illegal value of " + xyz.substr(idim,1) + " (min >= max)",
      kNamespaceName, "CheckDimension");
    result = false;
  }

  // Check function
  if ( ( info.fFcnName != "none" ) && ( info.fBinScheme != G4BinScheme::kLinear ) ) {
    Warn("Combining  " + xyz.substr(idim,1) + " Function and Binning scheme is not supported.",
      kNamespaceName, "CheckDimension");
    result = false;
  }

  // Check minValue if log binning or log function
  if ( ( info.fBinScheme == G4BinScheme::kLog ||
         info.fFcnName == "log" || info.fFcnName == "log10" ) && ( dimension.fMinValue == 0 ) ) {
    Warn("Illegal value of " + xyz.substr(idim,1) + " (min = 0) with logarithmic function or binning",
      kNamespaceName, "CheckDimension");
    result = false;
  }

  return result;
}

}
