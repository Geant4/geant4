//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Randomize.hh,v 1.5 2001-07-11 10:00:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

// Additional distributions ...
//
#include <CLHEP/Random/RandGaussQ.h>
#include <CLHEP/Random/RandGaussT.h>
#include <CLHEP/Random/RandPoissonQ.h>
#include <CLHEP/Random/RandPoissonT.h>
#include <CLHEP/Random/RandLandau.h>
#include <CLHEP/Random/RandBit.h>

#define G4RandGauss RandGaussQ
#define G4UniformRand() HepRandom::getTheEngine()->flat()

#endif // randomize_h 
