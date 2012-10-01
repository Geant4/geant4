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
// *                                                                  *
// * Parts of this code which have been  developed by QinetiQ Ltd     *
// * under contract to the European Space Agency (ESA) are the        *
// * intellectual property of ESA. Rights to use, copy, modify and    *
// * redistribute this software for general public use are granted    *
// * in compliance with any licensing, distribution and development   *
// * policy adopted by the Geant4 Collaboration. This code has been   *
// * written by QinetiQ Ltd for the European Space Agency, under ESA  *
// * contract 19770/06/NL/JD (Technology Research Programme).         *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr02/src/G4FullGlaubAADataSet.cc
/// \brief Implementation of the G4FullGlaubAADataSet class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4FullGlaubAADataSet.cc
//
// Version:             0.B
// Date:                02/04/08
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            19770/06/NL/JD
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#ifdef G4_USE_DPMJET


#include "G4FullGlaubAADataSet.hh"

#include "G4DPMJET2_5Interface.hh"

#include <iomanip>
///////////////////////////////////////////////////////////////////////////////
//
// G4FullGlaubAADataSet
//
// Constructor simply resets all variables to zero.
//
G4FullGlaubAADataSet::G4FullGlaubAADataSet()
{
  baseArrayPtrn = &bsiten[0][0];
  baseArrayPtrm = &bsitem[0][0];
  arrayPtrn     = baseArrayPtrn;
  arrayPtrm     = baseArrayPtrm;

  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j++) {
      bsiten[i][j] = 0.0;
      bsitem[i][j] = 0.0;
    }
  }

  wu10 = std::sqrt(10.0);        // WU10=SQRT(10.)

  glauberDataSetType = 0;
}
///////////////////////////////////////////////////////////////////////////////
//
// G4FullGlaubAADataSet
//
// Constructor resets all variables to zero and creates data.
//
G4FullGlaubAADataSet::G4FullGlaubAADataSet(const G4int AP1, const G4int AT1)
{
  baseArrayPtrn = &bsiten[0][0];
  baseArrayPtrm = &bsitem[0][0];
  arrayPtrn     = baseArrayPtrn;
  arrayPtrm     = baseArrayPtrm;

  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j++) {
      bsiten[i][j] = 0.0;
      bsitem[i][j] = 0.0;
    }
  }

  wu10 = std::sqrt(10.0);        // WU10=SQRT(10.)

  glauberDataSetType = 0;

  CreateGlauberData (AP1, AT1);
}
///////////////////////////////////////////////////////////////////////////////
//
// ~G4FullGlaubAADataSet
//
// If you thought the contructor was boring, the destructor is even worse!.
// It doesn't do anything.
//
G4FullGlaubAADataSet::~G4FullGlaubAADataSet()
{}
////////////////////////////////////////////////////////////////////////////////
//
G4bool G4FullGlaubAADataSet::CreateGlauberData (const G4int AP1, const G4int AT1)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"AT G4FullGlaubAADataSet::CreateGlauberData: BEFORE SHMAKI CALL"
           <<G4endl;
  }
#endif
  AP = AP1;
  AT = AT1;

  nucc_.ijproj   = 1;                                  // IJPROJ=1
  collis_.ijprox = 1;                                  // IJPROX=1
  nuccc_.jjproj  = 1;                                  // JJPROJ=1
  G4int jjprox   = 1;                                  // JJPROX=1
  G4int ishc     = 0;                                  // ISHC=0

  G4int ZP1 = stabZ[AP];
  G4int ZT1 = stabZ[AT];

  for (G4int ig = 0; ig < maxig; ig++)
  {
    G4double ppn = std::pow(wu10,ig+2);                // PPN = WU10**(IG+1)
    shmaki_ (&AP, &ZP1, &AT, &ZT1, &rptshm_.rproj, &rptshm_.rtarg, &ppn);
    for (G4int i=0; i<maxArray; i++)
    {
      bsiten[ig][i] = dshm_.bsite[i][1];
    }
    ishc++;                                            // ISHC = ISHC + 1
    if (ishc == 1) {                                   // IF(ISHC.EQ.1)THEN
      rproj = rptshm_.rproj;                           // RPROJJ(MATNUM) = RPROJ
      rtarg = rptshm_.rtarg;                           // RTARGG(MATNUM) = RTARG
      bmax  = dshm_.bmax;                              // BMAXX(MATNUM)  = BMAX
      bstep = dshm_.bstep;                             // BSTEPP(MATNUM) = BSTEP
    }
  }
  nucc_.ijproj   = 13;                                 // IJPROJ=13
  collis_.ijprox = 13;                                 // IJPROX=13
  nuccc_.jjproj  = 13;                                 // JJPROJ=13
  jjprox         = 13;                                 // JJPROX=13

  for (G4int ig = 0; ig < maxig; ig++)
  {
    G4double ppn = std::pow(wu10,ig+2);                // PPN = WU10**(IG+1)
    shmaki_ (&AP, &ZP1, &AT, &ZT1, &rptshm_.rproj, &rptshm_.rtarg, &ppn);
    for (G4int i=0; i<maxArray; i++)
    {
      bsitem[ig][i] = dshm_.bsite[i][1];
    }
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout <<"AT G4FullGlaubAADataSet::CreateGlauberData: AFTER SHMAKI CALL"
           <<G4endl;
  }
#endif

  return true;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double *G4FullGlaubAADataSet::GetArrayPointerN (const G4double ppn)
{
  if (ppn < 1.0E-10) {
    return baseArrayPtrn;
  }
  else {
    G4int ig = G4int(2.0*std::log10(ppn)) - 2;
    if (ig < 0) ig = 0;

    for (G4int j=0; j<maxArray; j++) {
      dtumat_.bsiten[0][ig][j] = bsiten[ig][j];
    }
    dtumat_.ntaxx[0]  = AT;
    dtumat_.nztaxx[0] = ZT;
    dtumat_.nprxx[0]  = AP;
    dtumat_.nzprxx[0] = ZP;
    dtumat_.rprojj[0] = rproj;
    dtumat_.rtagg[0]  = rtarg;
    dtumat_.bstepp[0] = bstep;
    dtumat_.bmaxx[0]  = bmax;
    arrayPtrn = baseArrayPtrn + maxArray*ig; //NOTE IS THIS A VALID OPERATION??
  }
  return arrayPtrn;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double *G4FullGlaubAADataSet::GetArrayPointerM (const G4double ppn)
{
  if (ppn < 1.0E-10) {
    return baseArrayPtrm;
  }
  else {
    G4int ig = G4int(2.0*std::log10(ppn)) - 2;
    if (ig < 0) ig = 0;

    for (G4int j=0; j<maxArray; j++) {
      dtumat_.bsitem[0][ig][j] = bsitem[ig][j];
    }
    dtumat_.ntaxx[0]  = AT;
    dtumat_.nztaxx[0] = ZT;
    dtumat_.nprxx[0]  = AP;
    dtumat_.nzprxx[0] = ZP;
    dtumat_.rprojj[0] = rproj;
    dtumat_.rtagg[0]  = rtarg;
    dtumat_.bstepp[0] = bstep;
    dtumat_.bmaxx[0]  = bmax;
    arrayPtrm = baseArrayPtrm + maxArray*ig; //NOTE IS THIS A VALID OPERATION??
  }
  return arrayPtrm;
}
////////////////////////////////////////////////////////////////////////////////
//
// WriteDataToFile
//
// This appends the Glauber data from the arrays to the output file stream.
// The format is intended to match the standard GLAUBER data file format.
//
std::ofstream & G4FullGlaubAADataSet::WriteDataToFile (std::ofstream &File) const
{
  File.unsetf(std::ios::fixed);
  File.setf(std::ios::scientific|std::ios::right|std::ios::adjustfield);
  
  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j+=5) {
      for (G4int k=j; k<j+5; k++) {
        File <<std::setw(16) <<bsiten[i][k];
      }
      File <<G4endl;
    }
  }

  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j+=5) {
      for (G4int k=j; k<j+5; k++) {
        File <<std::setw(16) <<bsitem[i][k];
      }
      File <<G4endl;
    }
  }

  return File;
}
////////////////////////////////////////////////////////////////////////////////
//
// DumpData
//
// This reads the Glauber data into the arrays from the input file stream.
// The format is intended to match the standard GLAUBER data file format.
//
std::ifstream & G4FullGlaubAADataSet::ReadDataFromFile (std::ifstream &File)
{
  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j+=5) {
      File >>bsiten[i][j]
           >>bsiten[i][j+1]
           >>bsiten[i][j+2]
           >>bsiten[i][j+3]
           >>bsiten[i][j+4];
    }
  }

  for (G4int i=0; i<maxig; i++) {
    for (G4int j=0; j<maxArray; j+=5) {
      File >>bsitem[i][j]
           >>bsitem[i][j+1]
           >>bsitem[i][j+2]
           >>bsitem[i][j+3]
           >>bsitem[i][j+4];
    }
  }

  return File;
}
////////////////////////////////////////////////////////////////////////////////
//
#endif
