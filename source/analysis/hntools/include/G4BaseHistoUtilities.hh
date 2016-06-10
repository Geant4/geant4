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
// $Id$

// Utility functions using tools::histo::base_histo information.
//
// Author: Ivana Hrivnacova, 22/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4BaseHistoUtilities_h
#define G4BaseHistoUtilities_h 1

#include "globals.hh"

#include "tools/histo/base_histo"

namespace G4Analysis
{

using G4ToolsBaseHisto
  = tools::histo::base_histo<double, unsigned int, unsigned int, double, double>;

// Access to data parameters
//  
G4int    GetNbins(const G4ToolsBaseHisto& baseHisto, G4int dimension);
G4double GetMin(const G4ToolsBaseHisto& baseHisto, G4int dimension);
G4double GetMax(const G4ToolsBaseHisto& baseHisto, G4int dimension);
G4double GetWidth(const G4ToolsBaseHisto& baseHisto, G4int dimension,
                  const G4String& hnType);

// Attributes for plotting
//

// Setters
G4bool SetTitle(G4ToolsBaseHisto& baseHisto, const G4String& title);
G4bool SetAxisTitle(G4ToolsBaseHisto& baseHisto, G4int dimension, const G4String& title);

// Accessors
G4String GetTitle(const G4ToolsBaseHisto& baseHisto);
G4String GetAxisTitle(const G4ToolsBaseHisto& baseHisto, G4int dimension, 
                      const G4String& hnType);

}

#endif

