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
// $Id: G4PStep.hh,v 1.3 2002-04-10 13:13:06 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PStep
//
// Class description:
//
// This class serves to address the "cell" a track previously 
// touched and a "cell" a track is currently in. It is used 
// for scoring and importance sampling in the "mass" geometry as well 
// as in a "parallel" geometry. 
// The "cell" information is available with the fPreTouchableKey and the 
// fPostTouchableKey.
// The fCrossBoundary member is set true in case the step
// crosses a boundary in the geometry it this G4PStep 
// refers to.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PStep_hh
#define G4PStep_hh G4PStep_hh

#include "G4PTouchableKey.hh"

class G4PStep
{

public:  // with description

  G4PStep(const G4PTouchableKey &preKey, const G4PTouchableKey &postKey)
   : fPreTouchableKey(preKey), fPostTouchableKey(postKey), 
     fCrossBoundary(false) {}
    // initialise pre and post G4PTouchableKey 

  ~G4PStep(){}

public:  // without description

  G4PTouchableKey fPreTouchableKey;
    // addressing the  "cell" the track previously touched 

  G4PTouchableKey fPostTouchableKey;  
    // addressing the current "cell"

  G4bool fCrossBoundary;
    // true if step crosses boundary of the geometry it refers to
};

#endif
