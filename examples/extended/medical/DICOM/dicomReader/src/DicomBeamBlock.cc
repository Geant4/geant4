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
#include "DicomBeamBlock.hh"
#include "dcmtk/dcmrt/seq/drtcos.h"
#include "G4UIcommand.hh"

// DOC at https://www.dabsoft.ch/dicom/3/C.8.8.14/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomBeamBlock::DicomBeamBlock(DRTBlockSequenceInRTBeamsModule::Item bblItem)
{
  OFString fstr;
  Sint32 fint;
  Float64 ffloat;
  OFVector<Float64> fvfloat;
  OFCondition cond; 
  G4cout << " DicomBeamBlock::DicomBeamBlock " << G4endl;

  bblItem.getBlockName(fstr);
  bblItem.getBlockType(fstr);
  bblItem.getBlockNumber(fint);
  theBlockNumber = fint;
  G4cout << " Number " << fint << G4endl;
  bblItem.getBlockTrayID(fstr);
  bblItem.getMaterialID(fstr);
  bblItem.getBlockNumberOfPoints(fint);
  bblItem.getBlockData(fvfloat);
  
  bblItem.getSourceToBlockTrayDistance(ffloat);
  theSourceToBlockTrayDistance = ffloat;

  bblItem.getBlockMountingPosition(fstr);
  theBlockMountingPosition = fstr;
  
  bblItem.getBlockDivergence(fstr);
  bblItem.getBlockThickness(ffloat);
  bblItem.getBlockTransmission(ffloat);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamBlock::Print( std::ostream&  )
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomBeamBlock::DumpToFile( std::ofstream& fout )
{
  std::string name  = ":P BLOCK_" + G4UIcommand::ConvertToString(theBlockNumber) + "_";

  fout << ":PS BLOCK_POS " << theBlockMountingPosition << G4endl;
  
  fout << name << "SourceToBlockTrayDistance " << theSourceToBlockTrayDistance << G4endl;

}
