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
// $Id$ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "G4AblaVirtualData.hh"

G4AblaVirtualData::G4AblaVirtualData()
{
	
}

bool G4AblaVirtualData::setAlpha(int A, int Z, double value)
{
  alpha[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setEcnz(int A, int Z, double value)
{
  ecnz[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setVgsld(int A, int Z, double value)
{
  vgsld[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setPace2(int A, int Z, double value)
{
  pace2[A][Z] = value;

  return true;
}

double G4AblaVirtualData::getAlpha(int A, int Z)
{
  return alpha[A][Z];
}

double G4AblaVirtualData::getEcnz(int A, int Z)
{
  return ecnz[A][Z];
}

double G4AblaVirtualData::getVgsld(int A, int Z)
{
  return vgsld[A][Z];
}

double G4AblaVirtualData::getPace2(int A, int Z)
{
  return pace2[A][Z];
}

int G4AblaVirtualData::getAlphaRows()
{
  return alphaRows;
}

int G4AblaVirtualData::getAlphaCols()
{
  return alphaCols;
}
int G4AblaVirtualData::getPaceRows()
{
  return paceRows;
}
int G4AblaVirtualData::getPaceCols()
{
  return paceCols;
}
