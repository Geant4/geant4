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
    if      ( binSchemeName == "log" )  binScheme = kLogBinScheme;
    // There is no name associated with kUserBinScheme
    else {
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
    
}
