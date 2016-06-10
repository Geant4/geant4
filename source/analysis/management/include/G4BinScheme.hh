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
// $Id: G4Fcn.hh 72292 2013-07-15 12:01:43Z ihrivnac $

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4BinScheme_h
#define G4BinScheme_h 1

#include "G4Fcn.hh"
#include "globals.hh"

#include <vector>

// Enumeration for definition available binning schemes

enum class G4BinScheme {
  kLinear,
  kLog,
  kUser
};  

// Utility function
namespace G4Analysis {

G4BinScheme GetBinScheme(const G4String& binSchemeName);

void ComputeEdges(G4int nbins, G4double xmin, G4double xmax, 
                  G4double unit, G4Fcn fcn, G4BinScheme,
                  std::vector<G4double>& edges);

void ComputeEdges(const std::vector<G4double>& edges, 
                  G4double unit, G4Fcn fcn, 
                  std::vector<G4double>& newEdges);
}

#endif  
