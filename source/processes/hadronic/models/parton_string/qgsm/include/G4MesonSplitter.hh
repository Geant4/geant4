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
#ifndef G4MesonSplitter_h
#define G4MesonSplitter_h

// HPW Feb 1999, based on Annihilator prototype.
// Simple class to split a meson, only one trivial method at the moment
// liable for improvement. @@@
// interfaces need change. @@@

#include "globals.hh"

class G4MesonSplitter
{
public:
  G4bool SplitMeson(G4int PDGcode, G4int* aEnd, G4int* bEnd);

private:

};

#endif
