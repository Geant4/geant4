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
// $Id: TiaraVDirectionGenerator.hh,v 1.3 2003/06/25 09:12:54 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
//
// Class TiaraVDirectionGenerator
//

#ifndef TiaraVDirectionGenerator_hh
#define TiaraVDirectionGenerator_hh TiaraVDirectionGenerator_hh

#include "globals.hh"
#include "G4ThreeVector.hh"


class TiaraVDirectionGenerator {
public:
  TiaraVDirectionGenerator();
  virtual ~TiaraVDirectionGenerator();

  virtual G4ThreeVector GetDirection() = 0;
  virtual TiaraVDirectionGenerator *Clone() const = 0;
};

#endif
