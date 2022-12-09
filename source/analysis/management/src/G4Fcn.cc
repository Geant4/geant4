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

#include "G4Fcn.hh"
#include "G4AnalysisUtilities.hh"

namespace G4Analysis
{

//_____________________________________________________________________________
G4Fcn GetFunction(const G4String& fcnName)
{
  // clang-format off
  G4Fcn fcn = G4FcnIdentity;
  if ( fcnName != "none" ) {
    if      ( fcnName == "log" )  fcn = std::log;
    else if ( fcnName == "log10") fcn = std::log10;
    else if ( fcnName == "exp" )  fcn = std::exp;
    else {
      Warn(
        "\"" + fcnName + "\" function is not supported.\n" +
        "No function will be applied to histogram values.",
        kNamespaceName, "GetFunction");
    }
  }
  return fcn;
  // clang-format on
}

}
