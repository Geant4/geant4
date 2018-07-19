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
#ifndef DicomBeamBlock__HH
#define DicomBeamBlock__HH

#include "dcmtk/dcmrt/seq/drtbl2.h"

#include <vector>
#include <iostream>

class DicomBeamBlock 
{ 
public:
  DicomBeamBlock(DRTBlockSequenceInRTBeamsModule::Item bcompItem);
  ~DicomBeamBlock(){};

public:

  void Print( std::ostream& out );

  void DumpToFile( std::ofstream& out );

private:
  OFVector<Float64> theBlockData;
  OFString theBlockDivergence;
  OFString theBlockMountingPosition;
  OFString theBlockName;
  Sint32 theBlockNumber;
  Sint32 theBlockNumberOfPoints;
  Float64 theBlockThickness;
  Float64 theBlockTransmission;
  OFString theBlockTrayID;
  OFString theBlockType;
  OFString theMaterialID;
  Float64 theSourceToBlockTrayDistance;

};

#endif  
