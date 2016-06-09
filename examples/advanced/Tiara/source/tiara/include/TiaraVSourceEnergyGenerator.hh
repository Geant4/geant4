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
// $Id: TiaraVSourceEnergyGenerator.hh,v 1.3 2003/06/25 09:12:55 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraVSourceEnergyGenerator
//

#ifndef TiaraVSourceEnergyGenerator_hh
#define TiaraVSourceEnergyGenerator_hh TiaraVSourceEnergyGenerator_hh

#include "globals.hh"

class TiaraVSourceEnergyGenerator {
public:
  TiaraVSourceEnergyGenerator();
  virtual ~TiaraVSourceEnergyGenerator();

  virtual G4double GetEnergy() = 0;
  virtual TiaraVSourceEnergyGenerator *Clone() const = 0;

};

#endif
