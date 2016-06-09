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
// $Id: G4DiffractiveHHScatterer.hh,v 1.8 2010-06-14 16:26:46 gcosmo Exp $

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
#include "G4FTFParameters.hh"
#include "G4ExcitedString.hh"

class G4DiffractiveHHScatterer
{
public:

   G4DiffractiveHHScatterer();
   virtual ~G4DiffractiveHHScatterer();

//   G4KineticTrackVector * Scatter(const G4KineticTrack & aTrack, const G4KineticTrack & bTrack);

   virtual void CreateStrings() const;
/*
                                        (G4VSplitableHadron * aHadron, 
                                         G4bool isProjectile,
                                         G4ExcitedString * FirstString, 
                                         G4ExcitedString * SecondString,
                                         G4FTFParameters *theParameters) const;
*/
private:

const G4DiffractiveExcitation * theExcitation;
G4LundStringFragmentation * theStringFragmentation;
G4FTFParameters  *theParameters;
};

#endif
