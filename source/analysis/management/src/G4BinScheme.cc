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
// $Id:$

// Author: Ivana Hrivnacova, 22/08/2013  (ivana@ipno.in2p3.fr)

#include "G4BinScheme.hh"

namespace G4Analysis
{

//_____________________________________________________________________________
G4BinScheme GetBinScheme(const G4String& binSchemeName)
{
  G4BinScheme binScheme = kLinearBinScheme;
  if ( binSchemeName != "linear" ) {
    if  ( binSchemeName == "log" )  
      binScheme = kLogBinScheme;
    else {
      // There is no name associated with kUserBinScheme
      G4ExceptionDescription description;
      description 
        << "    \"" << binScheme << "\" binning scheme is not supported." << G4endl
        << "    " << "Linear binning will be applied.";
      G4Exception("G4Analysis::GetBinScheme",
                "Analysis_W013", JustWarning, description);
    }              
  }
  return binScheme;            
}

//_____________________________________________________________________________
void ComputeEdges(G4int nbins, G4double xmin, G4double xmax, 
                  G4double unit, G4Fcn fcn, G4BinScheme binScheme,
                  std::vector<G4double>& edges)
{
// Compute edges from parameters

  // Apply units
  G4double xumin = xmin/unit;
  G4double xumax = xmax/unit;

  if ( binScheme == kLinearBinScheme ) {
    G4double dx = (fcn(xumax) - fcn(xumin) ) / nbins;
    G4double binValue = fcn(xumin);
    while ( G4int(edges.size()) <= nbins ) {
      edges.push_back(binValue);
      binValue += dx;
    }
  }  
  else if ( binScheme == kLogBinScheme ) {
    // do not apply fcn 
    G4double dlog 
      = (std::log10(xumax) - std::log10(xumin))/ nbins;
    G4double dx = std::pow(10, dlog);
    G4double binValue = xumin;
    while ( G4int(edges.size()) <= nbins ) {
      edges.push_back(binValue);
      binValue *= dx;
    }
  }
  else if ( binScheme == kUserBinScheme ) {  
    // This should never happen, but let's make sure about it
    // by issuing a warning
    G4ExceptionDescription description;
    description 
      << "    User binning scheme setting was ignored." << G4endl
      << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
    G4Exception("G4Analysis::ComputeEdges",
              "Analysis_W013", JustWarning, description);
  }              
}                                          

//_____________________________________________________________________________
void ComputeEdges(const std::vector<G4double>& edges, 
                  G4double unit, G4Fcn fcn, 
                  std::vector<G4double>& newBins)
{
// Apply function to defined edges
  std::vector<G4double>::const_iterator it;
  for (it = edges.begin(); it != edges.end(); it++ ) {
    newBins.push_back(fcn((*it)/unit));
  }
}
    
}
