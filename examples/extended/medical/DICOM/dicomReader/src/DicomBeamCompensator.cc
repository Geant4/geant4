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
#include "DicomBeamCompensator.hh"
#include "dcmtk/dcmrt/seq/drtcos.h"
#include "G4UIcommand.hh"

// DOC at https://www.dabsoft.ch/dicom/3/C.8.8.14/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeamCompensator::DicomBeamCompensator(DRTCompensatorSequence::Item bcompItem)
{
  OFString fstr;
  Sint32 fint;
  Float64 ffloat;
  OFVector<Float64> fvfloat;
  OFCondition cond; 
  G4cout << " DicomBeamCompensator::DicomBeamCompensator " << G4endl;
  cond = bcompItem.getCompensatorNumber(fint);
  theCompensatorNumber = fint;
  G4cout << " Number " << fint << G4endl;
    
  cond = bcompItem.getCompensatorColumns(fint);
  theCompensatorColumns = fint;
  cond = bcompItem.getCompensatorRows(fint);
  theCompensatorRows = fint;
  //  first value is the spacing between the center of adjacent rows, and the second value
  // (column spacing) is the spacing between the center of adjacent columns.  
  cond = bcompItem.getCompensatorPixelSpacing(fvfloat);
  theCompensatorPixelSpacing = fvfloat;

  cond = bcompItem.getCompensatorPosition(fvfloat);
  theCompensatorPosition = fvfloat;

  cond = bcompItem.getCompensatorTransmissionData(fvfloat);
  if( cond.good() ) theCompensatorTransmissionData = fvfloat;
  cond = bcompItem.getCompensatorThicknessData(fvfloat);
  if( cond.good() ) theCompensatorThicknessData = fvfloat;

  cond = bcompItem.getCompensatorTrayID(fstr);
  cond = bcompItem.getCompensatorType(fstr);

  cond = bcompItem.getMaterialID(fstr);
  if( cond.good() ) theMaterialID = fstr;
  cond = bcompItem.getSourceToCompensatorDistance(fvfloat);
  if( cond.good() ) theSourceToCompensatorDistance = fvfloat;
  cond = bcompItem.getSourceToCompensatorTrayDistance(ffloat);
  theSourceToCompensatorTrayDistance = ffloat;

  cond = bcompItem.getCompensatorDescription(fstr);
  cond = bcompItem.getCompensatorDivergence(fstr);
  cond = bcompItem.getCompensatorID(fstr);
  cond = bcompItem.getCompensatorMountingPosition(fstr);
  cond = bcompItem.getAccessoryCode(fstr);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamCompensator::Print( std::ostream&  )
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamCompensator::DumpToFile( std::ofstream& fout )
{
  std::string name  = ":P COMP_" +G4UIcommand::ConvertToString(theCompensatorNumber) + "_";
  fout << name << "PixelSpacing_1 " << theCompensatorPixelSpacing[0] << G4endl;
  fout << name << "PixelSpacing_2 " << theCompensatorPixelSpacing[1] << G4endl;

  fout << name << "POSX " << theCompensatorPosition[0] << G4endl;
  fout << name << "POSY " << theCompensatorPosition[1] << G4endl;
  
  fout << name << "SourceToCompensatorTrayDistance " << theSourceToCompensatorTrayDistance <<G4endl;

  for( size_t ii = 0; ii < theSourceToCompensatorDistance.size(); ii++ ) {
    int iCol = ii%theCompensatorColumns;
    int iRow = ii/theCompensatorColumns;
    fout << name << "SourceToCompensatorDistance_" +G4UIcommand::ConvertToString(iRow) + "_" 
      + G4UIcommand::ConvertToString(iCol) << " " << theSourceToCompensatorDistance[ii] << G4endl;
  }

  /*  for( size_t ii = 0; ii < theCompensatorTransmissionData.size(); ii++ ) {
    int iCol = ii%theCompensatorColumns;
    int iRow = ii/theCompensatorColumns;
    fout << name << "Transmission_" +G4UIcommand::ConvertToString(iRow) + "_" 
       +G4UIcommand::ConvertToString(iCol) << " " << theCompensatorTransmissionData[ii] << G4endl;
  }

  for( size_t ii = 0; ii < theCompensatorThicknessData.size(); ii++ ) {
    int iCol = ii%theCompensatorColumns;
    int iRow = ii/theCompensatorColumns;
    fout << name << "Thickness_" +G4UIcommand::ConvertToString(iRow) + "_" 
  +G4UIcommand::ConvertToString(iCol) << " " << theCompensatorThicknessData[ii] << G4endl;
  }
  */
}
