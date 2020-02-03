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
// G4Nsplit_Weight
//
// Class description:
// 
// A class (struct) used by importance sampling. It contains the number
// of tracks a mother track should be split into and their weight.
 
// Author: Michael Dressel (CERN), 2002
// ----------------------------------------------------------------------
#ifndef G4NSPLIT_WEIGHT_HH
#define G4NSPLIT_WEIGHT_HH 1

#include "globals.hh"

class G4Nsplit_Weight
{
  public:

    G4int fN = 0;
      // number of tracks a mother track should be split into
      // including the mother track

    G4double fW = 0.0;
      // the weight to be given to the tracks
};

std::ostream& operator<<(std::ostream& out, const G4Nsplit_Weight& nw);

#endif
