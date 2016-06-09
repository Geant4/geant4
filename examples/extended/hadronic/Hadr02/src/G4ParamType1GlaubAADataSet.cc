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
/// \file hadronic/Hadr02/src/G4ParamType1GlaubAADataSet.cc
/// \brief Implementation of the G4ParamType1GlaubAADataSet class
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4ParamType1GlaubAADataSet.cc
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


#include "G4ParamType1GlaubAADataSet.hh"

using namespace std;
#include "G4DPMJET2_5Interface.hh"

#include <iomanip>
///////////////////////////////////////////////////////////////////////////////
//
// G4ParamType1GlaubAADataSet
//
// Constructor simply resets all variables to zero.
//
G4ParamType1GlaubAADataSet::G4ParamType1GlaubAADataSet()
{
  glauberDataSetType = 1;
}
///////////////////////////////////////////////////////////////////////////////
//
// G4ParamType1GlaubAADataSet
//
// Constructor instructing to generate Glauber data based on AP1 and AT1.
//
G4ParamType1GlaubAADataSet::G4ParamType1GlaubAADataSet(const G4int AP1,
  const G4int AT1)
{
  glauberDataSetType = 1;

  CreateGlauberData (AP1,AT1);
}
///////////////////////////////////////////////////////////////////////////////
//
// G4ParamType1GlaubAADataSet
//
// Constructor instructing to generate Glauber data based full glauber set.
//
G4ParamType1GlaubAADataSet::G4ParamType1GlaubAADataSet
  (G4FullGlaubAADataSet *fullGlauberDataSet)
{
  glauberDataSetType = 1;

  CreateGlauberData (fullGlauberDataSet);
  AP    = fullGlauberDataSet->AP;
  AT    = fullGlauberDataSet->AT;
  rproj = fullGlauberDataSet->rproj;
  rtarg = fullGlauberDataSet->rtarg;
  bstep = fullGlauberDataSet->bstep;
  bmax  = fullGlauberDataSet->bmax;
  
}
///////////////////////////////////////////////////////////////////////////////
//
// ~G4ParamType1GlaubAADataSet
//
// If you thought the contructor was boring, the destructor is even worse!.
// It doesn't do anything.
//
G4ParamType1GlaubAADataSet::~G4ParamType1GlaubAADataSet()
{}
////////////////////////////////////////////////////////////////////////////////
//
G4bool G4ParamType1GlaubAADataSet::CreateGlauberData
  (const G4int AP1, const G4int AT1)
{
//
//
// Create an full Glauber data set.
//
  G4FullGlaubAADataSet *fullGlauberDataSet = new G4FullGlaubAADataSet();
  if (fullGlauberDataSet->CreateGlauberData(AP1,AT1)) {
    CreateGlauberData(fullGlauberDataSet);
    AP    = AP1;
    AT    = AT1;
    rproj = fullGlauberDataSet->rproj;
    rtarg = fullGlauberDataSet->rtarg;
    bstep = fullGlauberDataSet->bstep;
    bmax  = fullGlauberDataSet->bmax;
    
    return true;
  }
  else {
    return false;
  }

  delete fullGlauberDataSet;
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool G4ParamType1GlaubAADataSet::CreateGlauberData
  (G4FullGlaubAADataSet *fullGlauberDataSet)
{
  G4double *ptrArrayN = fullGlauberDataSet->GetArrayPointerN ();
  G4double *ptrArrayM = fullGlauberDataSet->GetArrayPointerM ();
//
//
// Loop over the momenta.  Create hit parameters using GetFitParameters member
// function, send a pointer to the relevant part of bsiten and bsitem array in
// fullGlauberData object.  The parameters are automatically inserted into the
// calling parameters array by the calling routine.
//
  for (G4int i=0; i<maxig; i++) {
    GetFitParameters (ptrArrayN, &paramn[i][0]);
    ptrArrayN  += maxArray;
    mun1[i]     = paramn[i][5] / paramn[i][1];
    mun2[i]     = paramn[i][5] / paramn[i][3] - mun1[i];
    G4double c1 = std::log(paramn[i][0]);
    G4double c2 = std::log(paramn[i][2]);
    cn[i]       = std::exp((c2*paramn[i][1] - c1*paramn[i][3]) /
      (paramn[i][3]-paramn[i][1]));
    
    GetFitParameters (ptrArrayM, &paramm[i][0]);
    ptrArrayM  += maxArray;
    mum1[i]     = paramm[i][5] / paramm[i][1];
    mum2[i]     = paramm[i][5] / paramm[i][3] - mun1[i];
    c1          = std::log(paramm[i][0]);
    c2          = std::log(paramm[i][2]);
    cm[i]       = std::exp((c2*paramm[i][1] - c1*paramm[i][3]) /
      (paramm[i][3]-paramm[i][1]));
  }
  
  return true;
}
////////////////////////////////////////////////////////////////////////////////
//
/*G4double *G4ParamType1GlaubAADataSet::GetArrayPointerN (const G4double ppn)
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if (ig > 23) ig = 23; 
  
  for (G4int j=0; j<maxArray; j++) {
    bsiten[j]                = GetInverseValueN(j,ig);
    dtumat_.bsiten[0][ig][j] = bsiten[j];
  }
  dtumat_.ntaxx[0]  = AT;
  dtumat_.nztaxx[0] = ZT;
  dtumat_.nprxx[0]  = AP;
  dtumat_.nzprxx[0] = ZP;
  dtumat_.rprojj[0] = rproj;
  dtumat_.rtagg[0]  = rtarg;
  dtumat_.bstepp[0] = bstep;
  dtumat_.bmaxx[0]  = bmax;

  arrayPtrn = baseArrayPtrn;
  
  return arrayPtrn;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double *G4ParamType1GlaubAADataSet::GetArrayPointerM (const G4double ppn)
{
  G4int ig = 0;
  if (ppn < 1.0E-10) {
    return 0;
  }
  else {
    ig = G4int(2.0*std::log10(ppn)) - 2;
  }
  if (ig > 23) ig = 23; 
  
  for (G4int j=0; j<maxArray; j++) {
    bsitem[j]                = GetInverseValueM(j,ig);
    dtumat_.bsitem[0][ig][j] = bsitem[j];
  }
  dtumat_.ntaxx[0]  = AT;
  dtumat_.nztaxx[0] = ZT;
  dtumat_.nprxx[0]  = AP;
  dtumat_.nzprxx[0] = ZP;
  dtumat_.rprojj[0] = rproj;
  dtumat_.rtagg[0]  = rtarg;
  dtumat_.bstepp[0] = bstep;
  dtumat_.bmaxx[0]  = bmax;

  arrayPtrm = baseArrayPtrm;
  
  return arrayPtrm;
}*/
////////////////////////////////////////////////////////////////////////////////
//
// WriteDataToFile
//
// This appends the Glauber data from the arrays to the output file stream.
// The format is intended to match the standard GLAUBER data file format.
//
std::ofstream & G4ParamType1GlaubAADataSet::WriteDataToFile (std::ofstream &File) const
{
  File.unsetf(std::ios::fixed);
  File.setf(std::ios::scientific|std::ios::right|std::ios::adjustfield);
  
  File <<"              c1"
       <<"              m1"
       <<"              c2"
       <<"              m2"
       <<"          Itrans"
       <<"           gamma"
       <<"              d0"
       <<"              d1"
       <<"              d2"
       <<"              d3"
       <<G4endl;

  for (G4int i=0; i<maxig; i++) {
    File <<std::setw(16) <<paramn[i][0]
         <<std::setw(16) <<paramn[i][1]
         <<std::setw(16) <<paramn[i][2]
         <<std::setw(16) <<paramn[i][3]
         <<std::setw(16) <<paramn[i][4]
         <<std::setw(16) <<paramn[i][5]
         <<std::setw(16) <<paramn[i][6]
         <<std::setw(16) <<paramn[i][7]
         <<std::setw(16) <<paramn[i][8]
         <<std::setw(16) <<paramn[i][9]
         <<G4endl;
  }
  
  for (G4int i=0; i<maxig; i++) {
    File <<std::setw(16) <<paramm[i][0]
         <<std::setw(16) <<paramm[i][1]
         <<std::setw(16) <<paramm[i][2]
         <<std::setw(16) <<paramm[i][3]
         <<std::setw(16) <<paramm[i][4]
         <<std::setw(16) <<paramm[i][5]
         <<std::setw(16) <<paramm[i][6]
         <<std::setw(16) <<paramm[i][7]
         <<std::setw(16) <<paramm[i][8]
         <<std::setw(16) <<paramm[i][9]
         <<G4endl;
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
std::ifstream & G4ParamType1GlaubAADataSet::ReadDataFromFile (std::ifstream &File)
{
  G4String dummy[10];
  File >>dummy[0]
       >>dummy[1]
       >>dummy[2]
       >>dummy[3]
       >>dummy[4]
       >>dummy[5]
       >>dummy[6]
       >>dummy[7]
       >>dummy[8]
       >>dummy[9];

  for (G4int i=0; i<maxig; i++) {
    File >>paramn[i][0]
         >>paramn[i][1]
         >>paramn[i][2]
         >>paramn[i][3]
         >>paramn[i][4]
         >>paramn[i][5]
         >>paramn[i][6]
         >>paramn[i][7]
         >>paramn[i][8]
         >>paramn[i][9];
  }
  
  for (G4int i=0; i<maxig; i++) {
    File >>paramm[i][0]
         >>paramm[i][1]
         >>paramm[i][2]
         >>paramm[i][3]
         >>paramm[i][4]
         >>paramm[i][5]
         >>paramm[i][6]
         >>paramm[i][7]
         >>paramm[i][8]
         >>paramm[i][9];
  }
  
  return File;
}
////////////////////////////////////////////////////////////////////////////////
//
#endif
