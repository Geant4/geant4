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
#include "DicomBeamDeviceRef.hh"
#include "dcmtk/dcmrt/seq/drtrbs8.h" // DRTReferencedBeamSequenceInRTFractionSchemeModule
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeamDeviceRef::DicomBeamDeviceRef(DRTBeamLimitingDeviceSequenceInRTBeamsModule::Item bldsItem)
{
  OFString fstr;
  Sint32 fint;
  Float64 ffloat;
  OFVector<Float64> fvfloat;

  bldsItem.getRTBeamLimitingDeviceType(fstr);
  G4cout << "   " << " RTBeamLimitingDeviceType " << fstr << G4endl;
  SetType(fstr);
  bldsItem.getSourceToBeamLimitingDeviceDistance(ffloat);
  G4cout << "   " << " SourceToBeamLimitingDeviceDistance " << ffloat << G4endl;
  SetSourceToBeamLimitingDeviceDistance( ffloat ); 
  bldsItem.getNumberOfLeafJawPairs(fint);
  SetNumberOfLeafJawPairs(fint);
  G4cout << "   " << " NumberOfLeafJawPairs " << fint << G4endl;
  bldsItem.getLeafPositionBoundaries(fvfloat);
  if( fint != 1 ) fint++;
  for( int ii = 0; ii < fint; ii++ ) {
    G4cout << "   " << ii << " LeafPositionBoundaries " << fvfloat[ii] << G4endl;
    AddPositionBoundary(fvfloat[ii]);
  }
}

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamDeviceRef::DumpToFile( std::ofstream& fout )
{

  fout << ":P " << theType << "_Z " << theSourceToBeamLimitingDeviceDistance << G4endl;

  /*  if( theType == "MLCX" || theType == "MLCY" ) {
    G4int nLeafs = theNumberOfLeafJawPairs;
    for( G4int jj = 0; jj < nLeafs; jj++ ){
      fout << ":P " << theType << "_" + std::to_string(jj+1) + "_CROSS "
           << thePositionBoundaries[jj] << G4endl;
    }
    }*/
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamDeviceRef::Print( std::ostream&  )
{

}
