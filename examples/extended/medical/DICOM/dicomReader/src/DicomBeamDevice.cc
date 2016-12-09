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
#include "DicomBeamDevice.hh"
#include "dcmtk/dcmrt/seq/drtrbs8.h" // DRTReferencedBeamSequenceInRTFractionSchemeModule
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeamDevice::DicomBeamDevice(DRTBeamLimitingDeviceSequenceInRTBeamsModule::Item bldsItem)
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
DicomBeamDevice::DicomBeamDevice(DRTBeamLimitingDevicePositionSequence::Item bldpsItem)
{
  OFString fstr;
  Float64 ffloat;
  
  bldpsItem.getRTBeamLimitingDeviceType(fstr);
  G4cout << "    " << " BeamLimitingDeviceType " << fstr << G4endl;
  SetType(fstr);
  for(size_t ii = 0;; ii++ ){
    if( bldpsItem.getLeafJawPositions(ffloat,ii) == EC_Normal ){
      G4cout << "    " << ii << " LeafPositionBoundaries " << ffloat << G4endl;
      AddPositionBoundary(ffloat);
    } else {
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamDevice::Print( std::ostream&  )
{

}
