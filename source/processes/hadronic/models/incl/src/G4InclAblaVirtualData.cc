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
// $Id: G4InclAblaVirtualData.cc,v 1.3 2010-06-14 16:10:01 gcosmo Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4InclAblaVirtualData.hh"

G4InclAblaVirtualData::G4InclAblaVirtualData() {}
G4InclAblaVirtualData::~G4InclAblaVirtualData() {}

bool G4InclAblaVirtualData::setAlpha(int A, int Z, double value)
{
  alpha[A][Z] = value;

  return true;
}

bool G4InclAblaVirtualData::setEcnz(int A, int Z, double value)
{
  ecnz[A][Z] = value;

  return true;
}

bool G4InclAblaVirtualData::setVgsld(int A, int Z, double value)
{
  vgsld[A][Z] = value;

  return true;
}

bool G4InclAblaVirtualData::setPace2(int A, int Z, double value)
{
  pace2[A][Z] = value;

  return true;
}

double G4InclAblaVirtualData::getAlpha(int A, int Z)
{
  return alpha[A][Z];
}

double G4InclAblaVirtualData::getEcnz(int A, int Z)
{
  return ecnz[A][Z];
}

double G4InclAblaVirtualData::getVgsld(int A, int Z)
{
  return vgsld[A][Z];
}

double G4InclAblaVirtualData::getPace2(int A, int Z)
{
  return pace2[A][Z];
}

int G4InclAblaVirtualData::getAlphaRows()
{
  return alphaRows;
}

int G4InclAblaVirtualData::getAlphaCols()
{
  return alphaCols;
}
int G4InclAblaVirtualData::getPaceRows()
{
  return paceRows;
}
int G4InclAblaVirtualData::getPaceCols()
{
  return paceCols;
}
