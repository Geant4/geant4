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
// $Id: G4VWeightWindowAlgorithm.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
