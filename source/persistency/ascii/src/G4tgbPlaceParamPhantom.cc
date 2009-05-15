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
//
// class G4tgbPlaceParamPhantom

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbPlaceParamPhantom.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4UIcommand.hh"

// -------------------------------------------------------------------------
G4tgbPlaceParamPhantom::~G4tgbPlaceParamPhantom()
{
}


// -------------------------------------------------------------------------
G4tgbPlaceParamPhantom::
G4tgbPlaceParamPhantom( G4tgrPlaceParameterisation* tgrParam ) 
{
  //---- Check number of extra data
  std::vector<G4double> extraData = tgrParam->GetExtraData();
  G4int ndata =  tgrParam->GetExtraData().size()+1;
  G4int nWcheck = 6;

  G4String outStr = G4String("PHANTOM") + " " + tgrParam->GetType() + " ";
  G4bool isOK = G4tgrUtils::CheckListSize( ndata, nWcheck, WLSIZE_EQ, outStr );

  if( !isOK )
  { 
    G4String chartmp = G4UIcommand::ConvertToString( nWcheck );
    outStr += chartmp + G4String(" words");
    G4cerr << outStr;
    G4cerr << " NUMBER OF WORDS " << ndata << G4endl;
    G4Exception("G4tgbPlaceParameterisation::CheckNExtraData",
                "InvalidData", FatalException, "Invalid data size.");
  }


  //  G4cout << " tgrParam->GetExtraData() " << extraData.size() << G4endl;
  G4int theNCopies1 = G4UIcommand::ConvertToInt( tgrParam->GetRotMatName() );
  //  G4int theNCopies1 = G4int(extraData[0]);
  G4int theNCopies2 = G4int(extraData[0]);
  G4int theNCopies3 = G4int(extraData[1]);
  theNCopies = theNCopies1 * theNCopies2 * theNCopies3;
  G4double theStep1 = extraData[2];
  G4double theStep2 = extraData[3];
  G4double theStep3 = extraData[4];
  theAxis = kXAxis;
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "G4tgbPlaceParamPhantom: no copies "
           << theNCopies << " = " << theNCopies1
           << " X " << theNCopies2 
           << " X " << theNCopies3 << G4endl
           << " step1 " << theStep1 << G4endl
           << " step2 " << theStep2 << G4endl
           << " step3 " << theStep3 << G4endl;
#endif


  SetVoxelDimensions( theStep1/2., theStep2/2., theStep3/2. );
  SetNoVoxel( theNCopies1, theNCopies2, theNCopies3 );
  SetSkipEqualMaterials(0);

}

