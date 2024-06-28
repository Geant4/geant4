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
//
// -------------------------------------------------------------------
//
//      Author:        E.Mendoza
// 
//      Creation date: May 2024
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  NuDEX code (https://doi.org/10.1016/j.nima.2022.167894)
// 


#ifndef NUDEXRANDOM_HH
#define NUDEXRANDOM_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

//COMPILATIONTYPE==1 compile with ROOT
//COMPILATIONTYPE==2 compile with GEANT4

#define COMPILATIONTYPE 2

#if COMPILATIONTYPE == 1
//------------------------------------------------------------
// ROOT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TRandom2.h"
#pragma GCC diagnostic pop
//------------------------------------------------------------
#elif COMPILATIONTYPE == 2
//------------------------------------------------------------
// GEANT4
#include "Randomize.hh"
#include "globals.hh"
#include "G4Exception.hh"
//------------------------------------------------------------
#else
  #error Unsupported COMPILATIONTYPE setting
#endif

void NuDEXException(const char* originOfException,const char* exceptionCode,const char* description);

class G4NuDEXRandom{

public:
  G4NuDEXRandom(unsigned int seed);
  ~G4NuDEXRandom();

public:
  void SetSeed(unsigned int seed);
  unsigned int GetSeed();
  G4double Uniform(G4double Xmin=0,G4double Xmax=1);
  unsigned int Integer(unsigned int IntegerMax);
  G4double Exp(G4double tau);
  G4double Gaus(G4double mean=0,G4double sigma=1);
  G4long Poisson(G4double mean);

private:

#if COMPILATIONTYPE == 1
  TRandom2* theRandom;
#elif COMPILATIONTYPE == 2
  CLHEP::HepJamesRandom* theEngine;
  CLHEP::RandFlat* theRandFlat;
  CLHEP::RandExponential* theRandExponential;
  CLHEP::RandGauss* theRandGauss;
  CLHEP::RandPoisson* theRandPoisson;
#endif
};


#endif




