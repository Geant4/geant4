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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DiffractiveHHScatterer.hh,v 1.4 2001/08/01 17:05:44 hpw Exp $

#ifndef G4DiffractiveHHScatterer_h
#define G4DiffractiveHHScatterer_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
class G4DiffractiveExcitation;
class G4LundStringFragmentation;
class G4KineticTrack;
#include "G4KineticTrackVector.hh"

class G4DiffractiveHHScatterer
{
public:

   G4DiffractiveHHScatterer();
   
   G4KineticTrackVector * Scatter(const G4KineticTrack & aTrack, const G4KineticTrack & bTrack);

private:

const G4DiffractiveExcitation * theExcitation;
G4LundStringFragmentation * theStringFragmentation;
};

#endif
