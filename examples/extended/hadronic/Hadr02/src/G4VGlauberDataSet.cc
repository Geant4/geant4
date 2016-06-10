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
/// \file hadronic/Hadr02/src/G4VGlauberDataSet.cc
/// \brief Implementation of the G4VGlauberDataSet class
//
// $Id: G4VGlauberDataSet.cc 77519 2013-11-25 10:54:57Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4VGlauberDataSet.cc
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


#include "G4VGlauberDataSet.hh"

#include "G4DPMJET2_5Interface.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//
// G4G4VGlauberDataSet
//
// Constructor simply resets all variables to zero.
//
G4VGlauberDataSet::G4VGlauberDataSet()
{
  rproj = 0.0;
  rtarg = 0.0;
  bstep = 0.0;
  bmax  = 0.0;
  AP    = 0;
  ZP    = 0;
  AT    = 0;
  ZT    = 0;
  
  maxArray = 200;
  maxig    = 24;

  glauberDataSetType = -1;
  
  verboseLevel = 0;
}
///////////////////////////////////////////////////////////////////////////////
//
// ~G4G4VGlauberDataSet
//
// If you thought the contructor was boring, the destructor is even worse!.
// It doesn't do anything.
//
G4VGlauberDataSet::~G4VGlauberDataSet()
{}
////////////////////////////////////////////////////////////////////////////////
//
G4int G4VGlauberDataSet::GetAP () const
{
  return AP;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int G4VGlauberDataSet::GetZP () const
{
  return ZP;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int G4VGlauberDataSet::GetAT () const
{
  return AT;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int G4VGlauberDataSet::GetZT () const
{
  return ZT;
}
////////////////////////////////////////////////////////////////////////////////
//
void G4VGlauberDataSet::SetArrayPointer (const G4int i)
{
  if (i<0 || i>=maxArray) {
    G4cerr <<"WARNING G4G4VGlauberDataSet::SetArrayPointer" <<G4endl;
    G4cerr <<"ATTEMPT TO SET POINTER TO VALUE OUTSIDE [0,"
           <<maxArray-1
           <<"]"
           <<G4endl;
  }
  else {
    arrayPtrn = baseArrayPtrn + i;
    arrayPtrm = baseArrayPtrm + i;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4double * G4VGlauberDataSet::GetArrayPointerN (const G4double ppn)
{
  if (ppn < 1.0E-10) return baseArrayPtrn;
  else               return arrayPtrn;
}
////////////////////////////////////////////////////////////////////////////////
//
G4double * G4VGlauberDataSet::GetArrayPointerM (const G4double ppn)
{
  if (ppn < 1.0E-10) return baseArrayPtrm;
  else               return arrayPtrm;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int G4VGlauberDataSet::GetGlauberDataSetType () const
{
  return glauberDataSetType;
}
////////////////////////////////////////////////////////////////////////////////
//
std::ofstream & G4VGlauberDataSet::WriteDataToFile (std::ofstream &File) const
{
//
//
// Dummy member function;
//
  return File;
}
std::ifstream & G4VGlauberDataSet::ReadDataFromFile (std::ifstream &File)
{
//
//
// Dummy member function;
//
  return File;
}
///////////////////////////////////////////////////////////////////////////////
//
// operator <<
//
// Output file-stream operator.  This is intended to match the standard
// GLAUBER data file format.
//
std::ofstream & operator << (std::ofstream &File, const G4VGlauberDataSet &q)
{
  File.unsetf(std::ios::scientific);
  File.setf(std::ios::fixed|std::ios::right|std::ios::adjustfield);
  File.precision(0);
  File <<std::setw(1) <<q.glauberDataSetType
       <<"NUCLEUS  "
       <<std::setw(10) <<q.AT
       <<std::setw(10) <<q.ZT
       <<std::setw(10) <<q.AP
       <<std::setw(10) <<q.ZP
       <<G4endl;

  File.unsetf(std::ios::fixed);
  File.setf(std::ios::fixed|std::ios::right|std::ios::adjustfield);
  File.precision(5);

  File <<std::setw(10) <<q.bmax
       <<std::setw(10) <<q.bstep
       <<std::setw(10) <<q.rproj
       <<std::setw(10) <<q.rtarg
       <<G4endl;

  File.precision(8);

  return q.WriteDataToFile (File);
}
///////////////////////////////////////////////////////////////////////////////
//
// operator >>
//
// Input file-stream operator.  This is assumed to the format matches the
// standard GLAUBER data file format.
//
std::ifstream & operator >> (std::ifstream &File, G4VGlauberDataSet &q)
{
  G4String dummy;
  File >>dummy
       >>q.AT
       >>q.ZT
       >>q.AP
       >>q.ZP;

  File >>q.bmax
       >>q.bstep
       >>q.rproj
       >>q.rtarg;

  return q.ReadDataFromFile (File);
}
#endif
