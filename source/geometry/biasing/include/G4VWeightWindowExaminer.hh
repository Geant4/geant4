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
// $Id: G4VWeightWindowExaminer.hh,v 1.2 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// Class G4VWeightWindowExaminer
//
// Class description:
//
// Interface class to a weight window examiner. The examiner
// has to return the number of tracks and their weight according to
// the initial track weight and the energy of the track.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VWeightWindowExaminer_hh
#define G4VWeightWindowExaminer_hh G4VWeightWindowExaminer_hh

#include "globals.hh"
#include "G4Nsplit_Weight.hh"

class G4VWeightWindowExaminer
{

public:  // with description
  
  // constructor and destructor
  G4VWeightWindowExaminer();
  virtual ~G4VWeightWindowExaminer();
  
  

  virtual G4Nsplit_Weight Examine(G4double w, G4double energy) const = 0; 
    // evaluate the number of tracks and their weight according
    // to the initial track weight and energy

};

#endif
