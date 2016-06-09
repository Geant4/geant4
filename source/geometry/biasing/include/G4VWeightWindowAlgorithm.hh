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
// $Id: G4VWeightWindowAlgorithm.hh,v 1.7 2003/08/19 15:44:57 dressel Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// Class G4VWeightWindowAlgorithm
//
// Class description:
// 
// Interface class for an weight window algorithm. It calculates
// the number of tracks and their weight according to the inital
// track weight and the lower energy bound in the energy-space
// cell.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VWeightWindowAlgorithm_hh
#define G4VWeightWindowAlgorithm_hh G4VWeightWindowAlgorithm_hh 

#include "G4Nsplit_Weight.hh"

class G4VWeightWindowAlgorithm
{

public:  // with description

  G4VWeightWindowAlgorithm();

  virtual ~G4VWeightWindowAlgorithm();

  virtual G4Nsplit_Weight Calculate(G4double init_w,
				    G4double lowerWeightBound) const = 0;
    // calculate number of tracks and their weight according
    // to the initial track weight and the lower energy bound

};

#endif
