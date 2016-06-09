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
// $Id: TiaraRandom.cc,v 1.3 2003/06/25 09:13:10 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include "TiaraRandom.hh"

#include "Randomize.hh"

void setRandomSeed(long seed) {
  HepRandom::setTheSeed(seed);
}

void setRandomStatus(const char *filename){
    HepRandom::restoreEngineStatus(filename);
}

void saveRandomStatus(const char *filename){
  HepRandom::saveEngineStatus(filename);
}


