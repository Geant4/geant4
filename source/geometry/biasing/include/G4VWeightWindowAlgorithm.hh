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
// $Id: G4VWeightWindowAlgorithm.hh,v 1.5 2002-10-14 12:36:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VWeightWindowAlgorithm
//
// Class description:
// 
// This is a base class for an algorithm used for 
// weight window importance sampling. It calculates
// the new weight of particles and the number of copies
// it should be split into according to the initial weight 
// of the track and the importance of the cell. The window
// is defined by the values of Upper and Lower.
// 
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

  virtual void SetUpperLimit(G4double Upper) = 0;
    // set upper limiting factor for window
  
  virtual void SetLowerLimit(G4double Lower) = 0;
    // set lower limiting facotr for window

  virtual G4Nsplit_Weight Calculate(G4double init_w, 
				    G4double importance) const = 0;
    // calculate the number of tracks and their weight according 
    // to the upper and lower limmiting factors of the window
    // and the initial weight
  
};

#endif
