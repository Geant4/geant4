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
#ifndef DicomBeamCompensator__HH
#define DicomBeamCompensator__HH

#include "dcmtk/dcmrt/seq/drtcos.h"

#include <vector>
#include <iostream>

class DicomBeamCompensator 
{ 
public:
  DicomBeamCompensator(DRTCompensatorSequence::Item bcompItem);
  ~DicomBeamCompensator(){};

public:

  void Print( std::ostream& out );

  void DumpToFile( std::ofstream& out );

private:
  OFString theAccessoryCode;
  Sint32 theCompensatorColumns;
  OFString theCompensatorDescription;
  OFString theCompensatorDivergence;
  OFString theCompensatorID;
  OFString theCompensatorMountingPosition;
  Sint32 theCompensatorNumber;
  OFVector<Float64> theCompensatorPixelSpacing;
  OFVector<Float64> theCompensatorPosition;
  Sint32 theCompensatorRows;
  OFVector<Float64> theCompensatorThicknessData;
  OFVector<Float64> theCompensatorTransmissionData;
  OFString theCompensatorTrayID;
  OFString theCompensatorType;
  OFString theMaterialID;
  OFVector<Float64> theSourceToCompensatorDistance;
  Float64 theSourceToCompensatorTrayDistance;
};

#endif  
