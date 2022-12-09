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
// G4tgbPlaceParameterisation implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgrPlaceParameterisation.hh"
#include "G4tgbPlaceParamCircle.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbRotationMatrixMgr.hh"
#include "G4UIcommand.hh"

// --------------------------------------------------------------------
G4tgbPlaceParameterisation::
G4tgbPlaceParameterisation(G4tgrPlaceParameterisation* tgrParam)
{
  theRotationMatrix =
    G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(
      tgrParam->GetRotMatName());
}

// --------------------------------------------------------------------
G4tgbPlaceParameterisation::~G4tgbPlaceParameterisation()
{
  delete theRotationMatrix;
}

// --------------------------------------------------------------------
void G4tgbPlaceParameterisation::ComputeTransformation(const G4int,
                                                       G4VPhysicalVolume*) const
{
}

// --------------------------------------------------------------------
void G4tgbPlaceParameterisation::
CheckNExtraData(G4tgrPlaceParameterisation* tgrParam, G4int nWcheck,
                WLSIZEtype st, const G4String& methodName)
{
  std::vector<G4double> extraData = tgrParam->GetExtraData();
  G4int ndata                     = (G4int)extraData.size();

  G4String outStr = methodName + " " + tgrParam->GetType() + " ";
  G4bool isOK     = G4tgrUtils::CheckListSize(ndata, nWcheck, st, outStr);

  if(!isOK)
  {
    G4String chartmp = G4UIcommand::ConvertToString(nWcheck);
    outStr += chartmp + G4String(" words");
    G4cerr << outStr;
    G4cerr << " NUMBER OF WORDS " << ndata << G4endl;
    G4Exception("G4tgbPlaceParameterisation::CheckNExtraData", "InvalidData",
                FatalException, "Invalid data size.");
  }
}
