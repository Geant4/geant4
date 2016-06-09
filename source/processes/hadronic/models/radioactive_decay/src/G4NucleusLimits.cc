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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4NucleusLimits.cc
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 14 April 2000, F Lei, DERA UK
// 0.b.4 release. No change in this class
//
//	02 November 2011, L Desorgher
//	Remove the maximum limits of A<250 and Z<100
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4NucleusLimits.hh"
////////////////////////////////////////////////////////////////////////////////
//
G4NucleusLimits::G4NucleusLimits ()
  : aMin(1), aMax(250), zMin(0), zMax(100)
  //
  // Default constructor sets the limits to cover all nuclei with Z<100.
  //
{;}
////////////////////////////////////////////////////////////////////////////////
//
G4NucleusLimits::G4NucleusLimits
  (G4int aMin1, G4int aMax1, G4int zMin1, G4int zMax1)
{
  //
  // check that aMin1 is within [1, 250].
  //
  if (aMin1 < 1) {aMin = 1;}
  //else if (aMin1 > 250) {aMin = 240;}
  else {aMin = aMin1;}

  //
  // check that aMax1 is within [1, 250].
  //
  if (aMax1 < 1) {aMax = 1;}
  //else if (aMax1 > 250) {aMax = 250;}
  else {aMax = aMax1;}

  //
  // If aMin > aMax then swap them.
  //
  if (aMin > aMax) {
    G4int t = aMin;
    aMin = aMax;
    aMax = t;
  }
  //
  // check that zMin1 is within [0, 100].
  //
  if (zMin1 < 0) {zMin = 0;}
  //else if (zMin1 > 100) {zMin = 100;}
  else {zMin = zMin1;}

  //
  // check that zMax1 is within [0, 100].
  //
  if (zMax1 < 0) {zMax = 0;}
  //else if (zMax1 > 100) {zMax = 100;}
  else {zMax = zMax1;}

  //
  // If zMin > zMax then swap them.
  //
  if (zMin > zMax) {
    G4int t = zMin;
    zMin = zMax;
    zMax = t;
  }
}


G4NucleusLimits::~G4NucleusLimits ()
{;}


std::ostream &operator << (std::ostream& strm, const G4NucleusLimits& q)
{
  // Definition of the insertion operator << to provide the nucleus limits to
  // ostream
  strm << "Atomic weight: " << q.GetAMin() << "->" << q.GetAMax()
       << "Atomic number: " << q.GetZMin() << "->" << q.GetZMax();
  return strm;
}

