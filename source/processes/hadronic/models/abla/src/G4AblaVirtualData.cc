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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, CEA (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//
#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4AblaVirtualData.hh"

#ifdef ABLAXX_IN_GEANT4_MODE
G4AblaVirtualData::G4AblaVirtualData() {}
#else
G4AblaVirtualData::G4AblaVirtualData(G4INCL::Config *) {}
#endif
G4AblaVirtualData::~G4AblaVirtualData() {}

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

bool G4AblaVirtualData::setRms(int A, int Z, double value)
{
  rms[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setMexp(int A, int Z, double value)
{
  mexp[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setMexpID(int A, int Z, int value)
{
  mexpid[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setBeta2(int A, int Z, double value)
{
  beta2[A][Z] = value;

  return true;
}

bool G4AblaVirtualData::setBeta4(int A, int Z, double value)
{
  beta4[A][Z] = value;

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

double G4AblaVirtualData::getRms(int A, int Z)
{
  return rms[A][Z];
}

double G4AblaVirtualData::getMexp(int A, int Z)
{
  return mexp[A][Z];
}

int G4AblaVirtualData::getMexpID(int A, int Z)
{
  return mexpid[A][Z];
}

double G4AblaVirtualData::getBeta2(int A, int Z)
{
  return beta2[A][Z];
}

double G4AblaVirtualData::getBeta4(int A, int Z)
{
  return beta4[A][Z];
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
